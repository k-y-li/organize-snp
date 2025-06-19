library(data.table)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(KEGGREST)
library(readr)
library(arrow)   
library(RSQLite) 
library(DBI)    

# --- I. Configuration & File Paths ---
# Path to GTEx v10 significant pairs file
LOCAL_GTEX_EQTL_FILE <- "/Users/Kevin/gtex data/GTEx_Analysis_v10_eQTL_updated/Brain_Cortex.v10.eQTLs.signif_pairs.parquet"

# Path to the large variant lookup table from GTEx Portal
LOOKUP_TABLE_GZ_FILE <- "/Users/Kevin/gtex data/GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz"

# Path to create/use the SQLite database for fast variant lookups
VARIANT_DB_FILE <- "gtex_variant_lookup.sqlite"

# Chromosomes to analyze (autosomal only)
TARGET_CHROMOSOMES <- paste0("chr", 1:22)

# P-value filtering configuration
PVALUE_COLUMN_FOR_FILTERING <- "pval_nominal"
PVALUE_THRESHOLD <- 5e-8 # A more standard threshold for genome-wide significance

# The name of the rsID column in the lookup file
RSID_COLUMN_NAME_IN_LOOKUP <- "rs_id_dbSNP155_GRCh38p13"

# Cache and Output file paths
KEGG_DATA_RDS_FILE <- "kegg_hsa_pathways_via_bioc_v10_hg38.rds"
OUTPUT_RDS_FILE <- "all_chromosomes_v10_brain_eqtls_kegg_pathways_COMBINED.rds"


# --- Ia. Create Variant Lookup Database (One-Time Setup) ---
message("\n--- Part 1: Setting up Variant Lookup Database ---")

if (!file.exists(VARIANT_DB_FILE)) {
  message(paste("Database not found. Creating", VARIANT_DB_FILE, "from", basename(LOOKUP_TABLE_GZ_FILE)))
  message("This is a one-time setup and may take several minutes...")
  
  if (!file.exists(LOOKUP_TABLE_GZ_FILE)) {
    stop(paste("FATAL: Variant lookup file not found at:", LOOKUP_TABLE_GZ_FILE))
  }
  
  con <- dbConnect(RSQLite::SQLite(), VARIANT_DB_FILE)
  
  read_tsv_chunked(
    LOOKUP_TABLE_GZ_FILE,
    callback = function(chunk, pos) {
      names(chunk)[names(chunk) == RSID_COLUMN_NAME_IN_LOOKUP] <- "rsID"
      dbWriteTable(con, "variants", as.data.frame(chunk), append = TRUE)
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


# --- II. Load and Prepare GTEx Data ---
message("\n--- Part 2: Loading GTEx v10 Significant Pairs data (GRCh38/hg38) ---")

if (!file.exists(LOCAL_GTEX_EQTL_FILE)) {
  stop(paste("GTEx v10 Parquet file not found at:", LOCAL_GTEX_EQTL_FILE))
}
message(paste("Reading GTEx v10 Parquet file:", LOCAL_GTEX_EQTL_FILE))

# Use arrow to read the parquet file and convert to data.table
full_gtex_data <- as.data.table(read_parquet(LOCAL_GTEX_EQTL_FILE))
message(paste("Loaded", format(nrow(full_gtex_data), big.mark=","), "total significant pairs."))

# Apply p-value filter
message(paste("Applying p-value filter:", PVALUE_COLUMN_FOR_FILTERING, "<=", PVALUE_THRESHOLD))
original_rows <- nrow(full_gtex_data)
full_gtex_data <- full_gtex_data[get(PVALUE_COLUMN_FOR_FILTERING) <= PVALUE_THRESHOLD & !is.na(get(PVALUE_COLUMN_FOR_FILTERING))]
message(sprintf("P-value filter: %s -> %s rows", format(original_rows, big.mark=","), format(nrow(full_gtex_data), big.mark=",")))

if (nrow(full_gtex_data) == 0) {
  stop("No eQTLs remaining after filtering. Consider relaxing the p-value threshold.")
}

# --- III. Load KEGG Pathway Data Once ---
message("\n--- Part 3: Loading KEGG Pathway Data ---")
kegg_pathways_df <- NULL
if (file.exists(KEGG_DATA_RDS_FILE)) {
  message("Loading cached KEGG pathway data from ", KEGG_DATA_RDS_FILE)
  kegg_pathways_df <- as.data.table(read_rds(KEGG_DATA_RDS_FILE))
} else {
  message("Constructing KEGG pathway-gene relationships using org.Hs.eg.db...")
  kegg_mappings <- tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = "ENTREZID"), 
                                                  columns = c("PATH", "SYMBOL"), keytype = "ENTREZID"), 
                            error = function(e) { message("Error: ", e$message); NULL })
  if (is.null(kegg_mappings) || nrow(kegg_mappings) == 0) stop("Could not get KEGG mappings from org.Hs.eg.db.")
  kegg_mappings_filtered <- kegg_mappings[!is.na(kegg_mappings$PATH), ]
  if (nrow(kegg_mappings_filtered) == 0) stop("No genes with KEGG pathway annotations in org.Hs.eg.db.")
  
  kegg_pathways_dt <- data.table(pathway_id_num = kegg_mappings_filtered$PATH, 
                                 entrez_gene_id_kegg = as.character(kegg_mappings_filtered$ENTREZID), 
                                 gene_symbol_from_bioc = kegg_mappings_filtered$SYMBOL)
  kegg_pathways_dt[, pathway_id := paste0("hsa", pathway_id_num)][, pathway_id_num := NULL]
  
  message("Fetching all KEGG pathway names for 'hsa'...")
  all_kegg_names_info <- tryCatch(keggList("pathway", "hsa"), 
                                  error = function(e) { 
                                    message("WARNING: Could not fetch KEGG pathway names via API: ", e$message); 
                                    NULL 
                                  })
  pathway_names_dt <- data.table()
  if (!is.null(all_kegg_names_info) && length(all_kegg_names_info) > 0) {
    pathway_names_dt <- data.table(pathway_id_full_kegg = names(all_kegg_names_info), 
                                   pathway_name = as.character(all_kegg_names_info))
    pathway_names_dt[, pathway_id := sub("path:", "", pathway_id_full_kegg)][, pathway_id_full_kegg := NULL]
  }
  
  kegg_pathways_df <- merge(kegg_pathways_dt, pathway_names_dt, by = "pathway_id", all.x = TRUE)
  kegg_pathways_df[is.na(pathway_name), pathway_name := "Name not available"]
  kegg_pathways_df <- kegg_pathways_df[, .(pathway_id, pathway_name, entrez_gene_id_kegg, 
                                           gene_symbol_kegg = gene_symbol_from_bioc)]
  write_rds(kegg_pathways_df, KEGG_DATA_RDS_FILE)
  message("Saved KEGG pathway data to ", KEGG_DATA_RDS_FILE)
}
message(paste("Loaded", uniqueN(kegg_pathways_df$pathway_id), "KEGG pathways."))


# --- IV. Function to Process One Chromosome ---
process_chromosome <- function(chr_name, gtex_data, kegg_data, db_path) {
  message(sprintf("\n--- Processing %s ---", chr_name))
  
  # 1. Filter GTEx data for the current chromosome
  chr_eqtls <- gtex_data[variant_id %like% paste0(chr_name, "_")]
  if (nrow(chr_eqtls) == 0) {
    message(sprintf("No significant eQTLs found for %s post-filtering.", chr_name))
    return(NULL)
  }
  
  # 2. Efficiently look up rsIDs from the SQLite database
  unique_variants <- unique(chr_eqtls$variant_id)
  con <- dbConnect(RSQLite::SQLite(), db_path)
  query <- dbSendQuery(con, "SELECT variant_id, rsID FROM variants WHERE variant_id = ?")
  dbBind(query, list(unique_variants))
  rsid_map <- dbFetch(query)
  dbClearResult(query)
  dbDisconnect(con)
  
  # 3. Merge rsIDs back into the eQTL data
  chr_eqtls <- merge(chr_eqtls, as.data.table(rsid_map), by = "variant_id", all.x = TRUE)
  
  # 4. Map Ensembl Gene IDs to Entrez IDs
  unique_ensembl_ids <- unique(sub("\\..*$", "", chr_eqtls$gene_id))
  entrez_map <- mapIds(org.Hs.eg.db, keys = unique_ensembl_ids, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  ensembl_to_entrez_dt <- data.table(ensembl_id_no_version = names(entrez_map[!is.na(entrez_map)]),
                                     entrez_gene_id = entrez_map[!is.na(entrez_map)])
  if (nrow(ensembl_to_entrez_dt) == 0) {
    message(sprintf("WARNING: No genes could be mapped to Entrez IDs for %s.", chr_name))
    return(NULL)
  }
  
  chr_eqtls[, ensembl_id_no_version := sub("\\..*$", "", gene_id)]
  chr_eqtls_with_entrez <- merge(chr_eqtls, ensembl_to_entrez_dt, by = "ensembl_id_no_version", all.x = TRUE)
  
  # 5. Merge with KEGG pathway data
  chr_eqtls_for_kegg <- chr_eqtls_with_entrez[!is.na(entrez_gene_id)]
  if (nrow(chr_eqtls_for_kegg) == 0) {
    message(sprintf("WARNING: No eQTLs with mappable Entrez IDs for %s.", chr_name))
    return(NULL)
  }
  
  chr_eqtls_for_kegg[, entrez_gene_id := as.character(entrez_gene_id)]
  merged_data <- merge(chr_eqtls_for_kegg, kegg_data, by.x = "entrez_gene_id", by.y = "entrez_gene_id_kegg", allow.cartesian = TRUE)
  if (nrow(merged_data) == 0) {
    message(sprintf("No KEGG pathway matches found for %s.", chr_name))
    return(NULL)
  }
  
  # 6. Create final, clean output table
  chr_output <- merged_data[, .(
    Chromosome = chr_name, Pathway_ID = pathway_id, Pathway_Name = pathway_name,
    Tissue = "Brain_Cortex_v10_hg38", SNP_chr_pos_ref_alt = variant_id, SNP_rsID = rsID,
    eQTL_Gene_Ensembl_Versioned = gene_id, eQTL_Gene_Symbol = gene_symbol_kegg,
    eQTL_Gene_Entrez = entrez_gene_id, eQTL_pval_nominal = pval_nominal,
    eQTL_slope = slope, eQTL_TSS_Distance = tss_distance
  )]
  
  chr_output <- unique(chr_output)
  message(sprintf("Final %s results: %s unique associations.", chr_name, format(nrow(chr_output), big.mark=",")))
  
  return(chr_output)
}


# --- V. Process All Chromosomes ---
message("\n--- Part 4: Processing All Chromosomes ---")
all_chromosome_results <- list()
summary_stats <- data.table()

for (chr in TARGET_CHROMOSOMES) {
  chr_result <- process_chromosome(chr, full_gtex_data, kegg_pathways_df, VARIANT_DB_FILE)
  if (!is.null(chr_result) && nrow(chr_result) > 0) {
    all_chromosome_results[[chr]] <- chr_result
    
    chr_summary <- data.table(
      Chromosome = chr,
      Total_Associations = nrow(chr_result),
      Unique_Pathways = uniqueN(chr_result$Pathway_ID),
      Unique_Genes = uniqueN(chr_result$eQTL_Gene_Entrez),
      Unique_Variants = uniqueN(chr_result$SNP_chr_pos_ref_alt),
      Variants_with_rsID = sum(!is.na(chr_result$SNP_rsID)),
      Percent_with_rsID = round(100 * sum(!is.na(chr_result$SNP_rsID)) / nrow(chr_result), 1)
    )
    summary_stats <- rbind(summary_stats, chr_summary)
  }
  gc()
}

# --- VI. Save Results to Separate CSV Files ---
message("\n--- Part 5: Saving results to separate CSV files ---")

if (length(all_chromosome_results) == 0) {
  warning("No results were generated for any chromosome. Skipping CSV file saving.")
} else {
  # Create a directory for the results if it doesn't exist
  if (!dir.exists("chromosome_results_csv")) {
    dir.create("chromosome_results_csv")
  }
  
  for (chr_name in names(all_chromosome_results)) {
    chr_data <- all_chromosome_results[[chr_name]]
    output_csv_name <- file.path("chromosome_results_csv", sprintf("Brain_Cortex_%s_eQTL_pathways.csv", chr_name))
    
    message(sprintf("Saving data for %s to: %s", chr_name, output_csv_name))
    fwrite(chr_data, output_csv_name)
  }
}


# --- VII. Save Combined RDS File ---
message("\n--- Part 6: Saving Combined RDS File for future R analysis ---")
if (length(all_chromosome_results) > 0) {
  combined_results <- rbindlist(all_chromosome_results, fill = TRUE)
  saveRDS(combined_results, OUTPUT_RDS_FILE)
  message(paste("All chromosomes combined and saved to:", OUTPUT_RDS_FILE))
}


# --- VIII. Final Summary ---
message("\n\n=== ANALYSIS COMPLETE ===")
message(sprintf("Analyzed chromosomes: %s", paste(TARGET_CHROMOSOMES, collapse = ", ")))
message(sprintf("Chromosomes with results: %d", length(all_chromosome_results)))
if (length(all_chromosome_results) > 0) {
  message(sprintf("Total associations across all chromosomes: %s", format(nrow(combined_results), big.mark=",")))
  message(sprintf("Total unique pathways: %d", uniqueN(combined_results$Pathway_ID)))
  message(sprintf("Total unique genes: %d", uniqueN(combined_results$eQTL_Gene_Entrez)))
  message(sprintf("Total unique variants: %d", uniqueN(combined_results$SNP_chr_pos_ref_alt)))
}

message(sprintf("\nOutput files created:"))
message(sprintf("- Separate CSV files in 'chromosome_results_csv/' directory"))
message(sprintf("- Combined RDS file:   %s", OUTPUT_RDS_FILE))

message("\nFinal Summary Statistics:")
print(summary_stats)