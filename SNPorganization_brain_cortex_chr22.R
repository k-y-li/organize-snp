library(readr)
library(arrow)
library(RSQLite)
library(DBI)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(KEGGREST)
library(dplyr)


# --- I. Configuration & File Paths ---
# Path to your GTEx v10 SIGNIFICANT PAIRS file (PARQUET format is best)
LOCAL_GTEX_EQTL_FILE <- "/Users/Kevin/gtex data/GTEx_Analysis_v10_eQTL_updated/Brain_Cortex.v10.eQTLs.signif_pairs.parquet"

# Path to the large variant lookup table from GTEx Portal
LOOKUP_TABLE_GZ_FILE <- "/Users/Kevin/gtex data/GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz"

# Path for the SQLite database
VARIANT_DB_FILE <- "gtex_variant_lookup.sqlite"

# This script is configured for a single chromosome.
TARGET_CHR_PREFIX_FOR_EQTL_FILTER <- "chr22_"

# P-value filtering configuration (for signif_pairs, pval_nominal is often used)
PVALUE_COLUMN_FOR_FILTERING <- "pval_nominal"  
PVALUE_THRESHOLD <- 5e-8 # Standard GWAS threshold

# The name of the rsID column in the lookup file
RSID_COLUMN_NAME_IN_LOOKUP <- "rs_id_dbSNP155_GRCh38p13"

# Cache file for KEGG pathway data
KEGG_DATA_RDS_FILE <- "kegg_hsa_pathways_via_bioc_v10_hg38.rds"

# Output files 
OUTPUT_CSV_FILE <- "chr22_v10_brain_eqtls_sqlite_final.csv"
OUTPUT_RDS_FILE <- "chr22_v10_brain_eqtls_sqlite_final.rds"


# --- II. Create Variant Lookup Database (One-Time Setup) ---
message("\n--- Part 1: Setting up Variant Lookup Database ---")

if (!file.exists(VARIANT_DB_FILE)) {
  message(paste("Database not found. Creating", VARIANT_DB_FILE, "from", basename(LOOKUP_TABLE_GZ_FILE)))
  message("This is a one-time setup and may take several minutes...")
  
  if (!file.exists(LOOKUP_TABLE_GZ_FILE)) {
    stop(paste("FATAL: Variant lookup file not found at:", LOOKUP_TABLE_GZ_FILE))
  }
  
  con <- dbConnect(RSQLite::SQLite(), VARIANT_DB_FILE)
  
  read_tsv_chunked(
    file = LOOKUP_TABLE_GZ_FILE,
    callback = function(chunk, pos) {
      names(chunk)[names(chunk) == RSID_COLUMN_NAME_IN_LOOKUP] <- "rsID"
      dbWriteTable(con, "variants", chunk, append = TRUE)
    },
    chunk_size = 500000,
    col_types = cols_only(
      variant_id = col_character(),
      !!RSID_COLUMN_NAME_IN_LOOKUP := col_character()
    ),
    show_col_types = FALSE
  )
  
  message("Creating index on 'variant_id' for fast lookups... (This is crucial!)")
  dbExecute(con, "CREATE UNIQUE INDEX idx_variant_id ON variants (variant_id);")
  dbDisconnect(con)
  message("Database setup and indexing complete.")
  
} else {
  message(paste("Found existing variant lookup database:", VARIANT_DB_FILE))
}


# --- III. Load and Filter Parquet Data ---
message("\n--- Part 2: Loading and Filtering GTEx Parquet Data ---")

if (!file.exists(LOCAL_GTEX_EQTL_FILE)) {
  stop(paste("GTEx Parquet file not found at:", LOCAL_GTEX_EQTL_FILE))
}

message(paste("Reading GTEx Parquet file:", LOCAL_GTEX_EQTL_FILE))
gtex_chr22_eqtls <- read_parquet(LOCAL_GTEX_EQTL_FILE) %>%
  filter(
    .data[[PVALUE_COLUMN_FOR_FILTERING]] <= PVALUE_THRESHOLD,
    !is.na(.data[[PVALUE_COLUMN_FOR_FILTERING]]),
    startsWith(variant_id, TARGET_CHR_PREFIX_FOR_EQTL_FILTER)
  )

message(paste(format(nrow(gtex_chr22_eqtls), big.mark=","), "eQTLs found for target chromosome after filtering."))

if (nrow(gtex_chr22_eqtls) == 0) {
  stop("No eQTLs found for the target chromosome after filtering. Consider relaxing filters.")
}


# --- IV. Add rsIDs via Efficient Database Lookup ---
message("\n--- Part 3: Adding rsIDs via SQLite Lookup ---")

variants_to_find <- unique(gtex_chr22_eqtls$variant_id)
message(sprintf("Querying database for rsIDs of %s unique variants...", format(length(variants_to_find), big.mark=",")))

con <- dbConnect(RSQLite::SQLite(), VARIANT_DB_FILE)
query <- dbSendQuery(con, "SELECT variant_id, rsID FROM variants WHERE variant_id = ?")
dbBind(query, list(variants_to_find))
rsid_map <- dbFetch(query)
dbClearResult(query)
dbDisconnect(con)

message(sprintf("Found %s matching rsIDs in the database.", format(nrow(rsid_map), big.mark=",")))

gtex_chr22_with_rsids <- left_join(gtex_chr22_eqtls, rsid_map, by = "variant_id")
rm(gtex_chr22_eqtls, rsid_map); gc()


# --- V. Map eGene Ensembl IDs to Entrez IDs ---
message("\n--- Part 4: Mapping eGene Ensembl IDs to Entrez IDs ---")

unique_egene_ensembl_ids_no_version <- unique(sub("\\..*$", "", gtex_chr22_with_rsids$gene_id))
message(paste("Found", length(unique_egene_ensembl_ids_no_version), "unique eGenes for Entrez mapping."))

entrez_map <- mapIds(org.Hs.eg.db, keys = unique_egene_ensembl_ids_no_version,
                     column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
ensembl_to_entrez_tbl <- tibble(
  ensembl_id_no_version = names(entrez_map[!is.na(entrez_map)]),
  entrez_gene_id = entrez_map[!is.na(entrez_map)]
)
message(paste("Mapped", nrow(ensembl_to_entrez_tbl), "Ensembl IDs to Entrez IDs."))

gtex_eqtls_with_entrez <- gtex_chr22_with_rsids %>%
  mutate(ensembl_id_no_version = sub("\\..*$", "", gene_id)) %>%
  left_join(ensembl_to_entrez_tbl, by = "ensembl_id_no_version") %>%
  filter(!is.na(entrez_gene_id))

message(paste(format(nrow(gtex_eqtls_with_entrez), big.mark=","), "eQTLs remaining for KEGG analysis."))
if (nrow(gtex_eqtls_with_entrez) == 0) stop("No eQTLs remaining after Entrez ID mapping.")
rm(gtex_chr22_with_rsids, ensembl_to_entrez_tbl); gc()


# --- VI. Get KEGG Pathway Data ---
message("\n--- Part 5: Fetching KEGG Pathway Data ---")
kegg_pathways_df <- NULL
if (file.exists(KEGG_DATA_RDS_FILE)) {
  message("Loading cached KEGG pathway data from ", KEGG_DATA_RDS_FILE)
  kegg_pathways_df <- read_rds(KEGG_DATA_RDS_FILE)
} else {
  message("Constructing KEGG pathway-gene relationships using org.Hs.eg.db...")
  kegg_mappings <- tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = "ENTREZID"), 
                                                  columns = c("PATH", "SYMBOL"), keytype = "ENTREZID"), 
                            error = function(e) { message("Error: ", e$message); NULL })
  if (is.null(kegg_mappings)) stop("Could not get KEGG mappings.")
  
  kegg_mappings_filtered <- as_tibble(kegg_mappings) %>% filter(!is.na(PATH))
  
  kegg_pathways_dt <- tibble(
    pathway_id_num = kegg_mappings_filtered$PATH, 
    entrez_gene_id_kegg = as.character(kegg_mappings_filtered$ENTREZID), 
    gene_symbol_from_bioc = kegg_mappings_filtered$SYMBOL
  ) %>% mutate(pathway_id = paste0("hsa", pathway_id_num)) %>% dplyr::select(-pathway_id_num)
  
  message("Fetching all KEGG pathway names for 'hsa'...")
  all_kegg_names_info <- tryCatch(keggList("pathway", "hsa"), error = function(e) { NULL })
  
  pathway_names_dt <- tibble()
  if (!is.null(all_kegg_names_info)) {
    pathway_names_dt <- tibble(
      pathway_id_full_kegg = names(all_kegg_names_info), 
      pathway_name = as.character(all_kegg_names_info)
    ) %>% mutate(pathway_id = sub("path:", "", pathway_id_full_kegg)) %>% dplyr::select(-pathway_id_full_kegg)
  }
  
  kegg_pathways_df <- left_join(kegg_pathways_dt, pathway_names_dt, by = "pathway_id") %>%
    mutate(pathway_name = ifelse(is.na(pathway_name), "Name not available", pathway_name)) %>%
    dplyr::select(pathway_id, pathway_name, entrez_gene_id_kegg, gene_symbol_kegg = gene_symbol_from_bioc)
  
  write_rds(kegg_pathways_df, KEGG_DATA_RDS_FILE)
  message("Saved KEGG pathway data to ", KEGG_DATA_RDS_FILE)
}
message(paste("Loaded", n_distinct(kegg_pathways_df$pathway_id), "KEGG pathways."))


# --- VII. Merge with KEGG and Finalize ---
message("\n--- Part 6: Merging with KEGG Pathways & Finalizing ---")

kegg_pathways_df <- mutate(kegg_pathways_df, entrez_gene_id_kegg = as.character(entrez_gene_id_kegg))
gtex_eqtls_with_entrez <- mutate(gtex_eqtls_with_entrez, entrez_gene_id = as.character(entrez_gene_id))

merged_data <- left_join(
  gtex_eqtls_with_entrez, 
  kegg_pathways_df, 
  by = c("entrez_gene_id" = "entrez_gene_id_kegg")
) %>%
  filter(!is.na(pathway_id))

if (nrow(merged_data) == 0) stop("No overlap found between eGenes and KEGG pathway genes.")

final_output_table <- merged_data %>%
  # Step 1: CREATE the new column with a constant value
  mutate(Tissue = "Brain_Cortex_v10_hg38") %>%
  # Step 2: SELECT and RENAME all desired columns (including the one you just made)
  dplyr::select(
    Pathway_ID = pathway_id,
    Pathway_Name = pathway_name,
    Tissue, # Now you can select the 'Tissue' column because it exists
    SNP_chr_pos_ref_alt = variant_id, 
    SNP_rsID = rsID,
    eQTL_Gene_Ensembl_Versioned = gene_id,
    eQTL_Gene_Ensembl_NoVersion = ensembl_id_no_version,
    eQTL_Gene_Entrez = entrez_gene_id,
    eQTL_Gene_Symbol = gene_symbol_kegg,
    any_of(c("pval_nominal", "slope", "tss_distance")) 
  ) %>%
  rename_with(~paste0("eQTL_", .), any_of(c("pval_nominal", "slope", "tss_distance"))) %>%
  distinct()

message(paste("\nFinal table has", format(nrow(final_output_table), big.mark=","), "unique rows."))


# --- VIII. Save Output ---
message("\n--- Part 7: Saving Output ---")
message(paste("Saving final table to CSV:", OUTPUT_CSV_FILE))
write_csv(final_output_table, OUTPUT_CSV_FILE)

message(paste("Saving final table to RDS:", OUTPUT_RDS_FILE))
write_rds(final_output_table, OUTPUT_RDS_FILE, compress = "gz")

message("\n--- Analysis Complete ---")