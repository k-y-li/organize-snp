# --- 0. Load Libraries ---
# Ensure these are installed:
# install.packages(c("data.table", "dplyr", "readr", "BiocManager"))
# BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db", "KEGGREST"))

library(data.table)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(KEGGREST)
library(readr)

# --- I. Configuration & File Paths ---
# Path to your GTEx v10 eGenes file (GRCh38/hg38)
LOCAL_GTEX_EQTL_FILE <- "/Users/Kevin/gtex data/GTEx_Analysis_v10_eQTL_updated/Brain_Cortex.v10.eGenes.txt.gz"

# CRITICAL: VERIFY this prefix by looking at your file's variant_id column
TARGET_CHR_PREFIX_FOR_EQTL_FILTER <- "chr22_"

# Optional TSS distance filter (GTEx eGenes are already significant)
MAX_TSS_DISTANCE <- 1000000

# P-value filtering configuration
PVALUE_COLUMN_FOR_FILTERING <- "pval_perm"  
PVALUE_THRESHOLD <- 0.05

# rsID column name from your GTEx v10 file
RSID_COLUMN_NAME <- "rs_id_dbSNP155_GRCh38p13"

# Cache file for KEGG pathway data (v10/hg38 compatible)
KEGG_DATA_RDS_FILE <- "kegg_hsa_pathways_via_bioc_v10_hg38.rds"

# Output files 
OUTPUT_CSV_FILE <- "chr22_v10_brain_eqtls_kegg_pathways_final_filtered.csv"
OUTPUT_RDS_FILE <- "chr22_v10_brain_eqtls_kegg_pathways_final_filtered.rds"

# --- II. Process and Filter GTEx v10 eGenes Data ---
message("--- Part 1: Processing GTEx v10 eGenes data (GRCh38/hg38) ---")

if (!file.exists(LOCAL_GTEX_EQTL_FILE)) {
  stop(paste("GTEx v10 eGenes file not found at:", LOCAL_GTEX_EQTL_FILE))
}
message(paste("Found GTEx v10 eGenes file:", LOCAL_GTEX_EQTL_FILE))
message(paste("File size:", round(file.info(LOCAL_GTEX_EQTL_FILE)$size / (1024^2), 2), "MB."))

gtex_final_filtered_eqtls <- NULL

tryCatch({
  message("Reading header of GTEx v10 eGenes file to check columns and variant_id format...")
  header_check_dt <- fread(LOCAL_GTEX_EQTL_FILE, nrows = 5)
  all_columns <- colnames(header_check_dt)
  message("Available columns: ", paste(all_columns, collapse=", "))
  message("First few variant_ids to check format:")
  print(header_check_dt$variant_id)
  
  # Check if rsID column exists
  if (RSID_COLUMN_NAME %in% all_columns) {
    message(paste("Found rsID column:", RSID_COLUMN_NAME))
  } else {
    warning(paste("rsID column", RSID_COLUMN_NAME, "not found. Will set rsIDs to NA."))
    RSID_COLUMN_NAME <- NULL
  }
  
  rm(header_check_dt); gc()
  
  # Select columns to read - include all p-value columns and key data
  cols_to_select <- c("gene_id", "gene_name", "variant_id", "tss_distance", 
                      "pval_nominal", "pval_perm", "pval_beta", "qval", 
                      "slope", "slope_se", "chr", "variant_pos", "ref", "alt")
  
  # Add rsID column if available
  if (!is.null(RSID_COLUMN_NAME)) {
    cols_to_select <- c(cols_to_select, RSID_COLUMN_NAME)
  }
  
  # Only keep columns that actually exist in the file
  cols_to_select <- cols_to_select[cols_to_select %in% all_columns]
  
  message(paste("\nReading selected columns (", paste(cols_to_select, collapse=", "), ") from GTEx v10 eGenes file..."))
  full_data_dt <- fread(LOCAL_GTEX_EQTL_FILE, select = cols_to_select, showProgress = TRUE)
  message(paste("Read", format(nrow(full_data_dt), big.mark=","), "rows (eGenes with lead variants)."))
  
  # Add rsID column if it exists, otherwise set to NA
  if (!is.null(RSID_COLUMN_NAME) && RSID_COLUMN_NAME %in% colnames(full_data_dt)) {
    setnames(full_data_dt, RSID_COLUMN_NAME, "rsID")
    message(paste("Using rsID column:", RSID_COLUMN_NAME))
    rsid_available_count <- sum(!is.na(full_data_dt$rsID))
    message(paste("rsIDs available for", format(rsid_available_count, big.mark=","), "out of", 
                  format(nrow(full_data_dt), big.mark=","), "eQTLs"))
  } else {
    full_data_dt[, rsID := NA_character_]
    message("No rsID column found - setting all rsIDs to NA")
  }
  
  # Apply TSS distance filter if column exists
  if ("tss_distance" %in% colnames(full_data_dt)) {
    message(paste("Original rows before TSS distance filter:", format(nrow(full_data_dt), big.mark=",")))
    full_data_dt <- full_data_dt[abs(tss_distance) <= MAX_TSS_DISTANCE]
    message(paste("Rows after TSS distance filter (abs(tss_distance) <= ", format(MAX_TSS_DISTANCE, big.mark=","), "): ", format(nrow(full_data_dt), big.mark=","), sep=""))
    if (nrow(full_data_dt) == 0) stop("No eQTLs remaining after TSS distance filter.")
  }
  
  # Apply p-value significance filter
  if (PVALUE_COLUMN_FOR_FILTERING %in% colnames(full_data_dt)) {
    message(paste("Original rows before p-value filter:", format(nrow(full_data_dt), big.mark=",")))
    pval_before_filter <- full_data_dt[[PVALUE_COLUMN_FOR_FILTERING]]
    message(sprintf("P-value distribution before filter (%s): min=%.2e, median=%.2e, max=%.2e", 
                    PVALUE_COLUMN_FOR_FILTERING, min(pval_before_filter, na.rm=TRUE), 
                    median(pval_before_filter, na.rm=TRUE), max(pval_before_filter, na.rm=TRUE)))
    
    full_data_dt <- full_data_dt[get(PVALUE_COLUMN_FOR_FILTERING) <= PVALUE_THRESHOLD & !is.na(get(PVALUE_COLUMN_FOR_FILTERING))]
    message(paste("Rows after p-value filter (", PVALUE_COLUMN_FOR_FILTERING, " <= ", PVALUE_THRESHOLD, "): ", format(nrow(full_data_dt), big.mark=","), sep=""))
    
    if (nrow(full_data_dt) == 0) {
      stop(paste("No eQTLs remaining after p-value filter. Consider using a less stringent threshold or different p-value column."))
    }
  } else {
    warning(paste("P-value column", PVALUE_COLUMN_FOR_FILTERING, "not found. Skipping p-value filtering."))
  }
  
  # Filter for chromosome 22
  message(paste("Filtering for chromosome 22 (prefix '", TARGET_CHR_PREFIX_FOR_EQTL_FILTER, "')...", sep=""))
  gtex_final_filtered_eqtls <- full_data_dt[startsWith(variant_id, TARGET_CHR_PREFIX_FOR_EQTL_FILTER)]
  message(paste(format(nrow(gtex_final_filtered_eqtls), big.mark=","), "rows after chromosome filter."))
  rm(full_data_dt); gc()
  
  if (nrow(gtex_final_filtered_eqtls) == 0) {
    stop("No eQTLs found for chromosome 22 after filtering. Check chromosome prefix.")
  }
  
  message("Filtered eQTL data (gtex_final_filtered_eqtls) prepared:")
  print(head(gtex_final_filtered_eqtls))
  
  # Show p-value statistics AFTER filtering
  message("\n=== P-VALUE STATISTICS (AFTER FILTERING) ===")
  pval_cols <- c("pval_nominal", "pval_perm", "pval_beta", "qval")
  pval_cols_present <- pval_cols[pval_cols %in% colnames(gtex_final_filtered_eqtls)]
  
  for (pval_col in pval_cols_present) {
    pval_data <- gtex_final_filtered_eqtls[[pval_col]]
    pval_data <- pval_data[!is.na(pval_data)]
    if (length(pval_data) > 0) {
      message(sprintf("%s: min=%.2e, median=%.2e, max=%.2e, <0.05: %d (%.1f%%)", 
                      pval_col, min(pval_data), median(pval_data), max(pval_data),
                      sum(pval_data < 0.05), 100*sum(pval_data < 0.05)/length(pval_data)))
    }
  }
  
  # Show filtering summary
  message(sprintf("\n=== FILTERING SUMMARY ==="))
  message(sprintf("Filtered by p-value: %s <= %.3f", PVALUE_COLUMN_FOR_FILTERING, PVALUE_THRESHOLD))
  message(sprintf("Final chromosome 22 eQTLs: %s", format(nrow(gtex_final_filtered_eqtls), big.mark=",")))
  if (PVALUE_COLUMN_FOR_FILTERING %in% colnames(gtex_final_filtered_eqtls)) {
    filtered_pvals <- gtex_final_filtered_eqtls[[PVALUE_COLUMN_FOR_FILTERING]]
    message(sprintf("Range of %s in final data: %.2e to %.2e", 
                    PVALUE_COLUMN_FOR_FILTERING, min(filtered_pvals, na.rm=TRUE), max(filtered_pvals, na.rm=TRUE)))
  }
  
}, error = function(e) {
  stop(paste("Error during eQTL data processing: ", e$message))
})

# --- III. Map eGene Ensembl IDs to Entrez IDs ---
message("\n--- Part 2: Mapping eGene Ensembl IDs to Entrez IDs ---")

if (is.null(gtex_final_filtered_eqtls) || nrow(gtex_final_filtered_eqtls) == 0) {
  stop("eQTL data ('gtex_final_filtered_eqtls') is empty. Cannot proceed with ID mapping.")
}

unique_egene_ensembl_ids_no_version <- unique(sub("\\..*$", "", gtex_final_filtered_eqtls$gene_id))
message(paste("Found", length(unique_egene_ensembl_ids_no_version), "unique eGenes for Entrez mapping."))

gtex_eqtls_for_kegg <- NULL

tryCatch({
  entrez_map <- mapIds(org.Hs.eg.db, keys = unique_egene_ensembl_ids_no_version,
                       column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  ensembl_to_entrez_dt <- data.table(ensembl_id_no_version = names(entrez_map[!is.na(entrez_map)]),
                                     entrez_gene_id = entrez_map[!is.na(entrez_map)])
  message(paste("Mapped", nrow(ensembl_to_entrez_dt), "Ensembl IDs to Entrez IDs."))
  if (nrow(ensembl_to_entrez_dt) == 0) stop("No eGenes mapped to Entrez IDs.")
  
  gtex_final_filtered_eqtls[, ensembl_id_no_version := sub("\\..*$", "", gene_id)]
  gtex_eqtls_with_entrez <- merge(gtex_final_filtered_eqtls, ensembl_to_entrez_dt, 
                                  by = "ensembl_id_no_version", all.x = TRUE)
  gtex_eqtls_for_kegg <- gtex_eqtls_with_entrez[!is.na(entrez_gene_id)] # This contains rsID from GTEx file
  message(paste(format(nrow(gtex_eqtls_for_kegg), big.mark=","), "eQTLs remaining for KEGG analysis."))
  if (nrow(gtex_eqtls_for_kegg) == 0) stop("No eQTLs remaining after Entrez ID mapping.")
}, error = function(e) { 
  stop(paste("Error during Ensembl to Entrez ID mapping: ", e$message)) 
})

# --- IV. Get KEGG Pathway Data ---
message("\n--- Part 3: Fetching KEGG Pathway Data ---")

kegg_pathways_df <- NULL
if (file.exists(KEGG_DATA_RDS_FILE)) {
  message("Loading cached KEGG pathway data from ", KEGG_DATA_RDS_FILE)
  kegg_pathways_df <- read_rds(KEGG_DATA_RDS_FILE)
} else {
  message("Constructing KEGG pathway-gene relationships using org.Hs.eg.db...")
  kegg_mappings <- tryCatch(select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = "ENTREZID"), 
                                   columns = c("PATH", "SYMBOL"), keytype = "ENTREZID"), 
                            error = function(e) { 
                              message("Error selecting from org.Hs.eg.db: ", e$message); 
                              NULL 
                            })
  if (is.null(kegg_mappings) || nrow(kegg_mappings) == 0) stop("Could not get KEGG mappings from org.Hs.eg.db.")
  kegg_mappings_filtered <- kegg_mappings[!is.na(kegg_mappings$PATH), ]
  if (nrow(kegg_mappings_filtered) == 0) stop("No genes with KEGG pathway annotations in org.Hs.eg.db.")
  
  kegg_pathways_dt <- data.table(pathway_id_num = kegg_mappings_filtered$PATH, 
                                 entrez_gene_id_kegg = as.character(kegg_mappings_filtered$ENTREZID), 
                                 gene_symbol_from_bioc = kegg_mappings_filtered$SYMBOL)
  kegg_pathways_dt[, pathway_id := paste0("hsa", pathway_id_num)][, pathway_id_num := NULL]
  message(paste("Found", format(nrow(kegg_pathways_dt), big.mark=","), "gene-pathway associations from org.Hs.eg.db."))
  
  message("Fetching all KEGG pathway names for 'hsa'...")
  all_kegg_names_info <- tryCatch(keggList("pathway", "hsa"), 
                                  error = function(e) { 
                                    message("WARNING: Could not fetch all KEGG pathway names via API: ", e$message); 
                                    NULL 
                                  })
  pathway_names_dt <- data.table()
  if (!is.null(all_kegg_names_info) && length(all_kegg_names_info) > 0) {
    pathway_names_dt <- data.table(pathway_id_full_kegg = names(all_kegg_names_info), 
                                   pathway_name = as.character(all_kegg_names_info))
    pathway_names_dt[, pathway_id := sub("path:", "", pathway_id_full_kegg)][, pathway_id_full_kegg := NULL]
    message(paste("Fetched names for", nrow(pathway_names_dt), "KEGG pathways."))
  } else { 
    unique_ids_from_data <- unique(kegg_pathways_dt$pathway_id)
    pathway_names_dt <- data.table(pathway_id = unique_ids_from_data, 
                                   pathway_name = paste("Name API call failed for", unique_ids_from_data))
    message("WARNING: Using placeholder pathway names as API call for all names failed.")
  }
  
  kegg_pathways_df <- merge(kegg_pathways_dt, pathway_names_dt, by = "pathway_id", all.x = TRUE)
  kegg_pathways_df[is.na(pathway_name), pathway_name := paste("Name not in API list for", pathway_id)]
  kegg_pathways_df <- kegg_pathways_df[, .(pathway_id, pathway_name, entrez_gene_id_kegg, 
                                           gene_symbol_kegg = gene_symbol_from_bioc)]
  write_rds(kegg_pathways_df, KEGG_DATA_RDS_FILE)
  message("Saved KEGG pathway data to ", KEGG_DATA_RDS_FILE)
}
if (is.null(kegg_pathways_df) || nrow(kegg_pathways_df) == 0) stop("KEGG pathway data frame is empty.")
message(paste("Loaded/Constructed", uniqueN(kegg_pathways_df$pathway_id), "KEGG pathways."))

# --- V. Merge eQTLs with KEGG Pathways and Prepare Final Output ---
message("\n--- Part 4: Merging eQTLs with KEGG Pathways & Finalizing ---")

if (is.null(gtex_eqtls_for_kegg) || nrow(gtex_eqtls_for_kegg) == 0) stop("eQTL data for KEGG ('gtex_eqtls_for_kegg') is empty.")
gtex_eqtls_for_kegg[, entrez_gene_id := as.character(entrez_gene_id)]
merged_data <- merge(gtex_eqtls_for_kegg, kegg_pathways_df, 
                     by.x = "entrez_gene_id", by.y = "entrez_gene_id_kegg", allow.cartesian = TRUE)
if (nrow(merged_data) == 0) stop("No overlap found between eGenes and KEGG pathway genes.")
message(paste("Found", format(nrow(merged_data), big.mark=","), "eQTL-gene-pathway associations."))

# Create final output table with all p-value columns
final_output_table <- merged_data[, .(
  Pathway_ID = pathway_id,
  Pathway_Name = pathway_name,
  Tissue = "Brain_Cortex_v10_hg38",
  SNP_chr_pos_ref_alt = variant_id, 
  SNP_rsID = rsID,         
  eQTL_Gene_Ensembl_Versioned = gene_id,
  eQTL_Gene_Ensembl_NoVersion = ensembl_id_no_version,
  eQTL_Gene_Entrez = entrez_gene_id,
  eQTL_Gene_Symbol = gene_symbol_kegg,
  eQTL_Gene_Name = if("gene_name" %in% colnames(.SD)) gene_name else NA_character_,
  # Include ALL p-value columns
  eQTL_pval_nominal = if("pval_nominal" %in% colnames(.SD)) pval_nominal else NA_real_,
  eQTL_pval_perm = if("pval_perm" %in% colnames(.SD)) pval_perm else NA_real_,
  eQTL_pval_beta = if("pval_beta" %in% colnames(.SD)) pval_beta else NA_real_,
  eQTL_qval = if("qval" %in% colnames(.SD)) qval else NA_real_,
  # Effect size information
  eQTL_slope = if("slope" %in% colnames(.SD)) slope else NA_real_, 
  eQTL_slope_se = if("slope_se" %in% colnames(.SD)) slope_se else NA_real_,
  eQTL_TSS_Distance = if("tss_distance" %in% colnames(.SD)) tss_distance else NA_integer_,
  # Variant position information
  Variant_Chr = if("chr" %in% colnames(.SD)) chr else NA_character_,
  Variant_Pos = if("variant_pos" %in% colnames(.SD)) variant_pos else NA_integer_,
  Variant_Ref = if("ref" %in% colnames(.SD)) ref else NA_character_,
  Variant_Alt = if("alt" %in% colnames(.SD)) alt else NA_character_
)]

final_output_table <- unique(final_output_table)
message(paste("\nFinal table has", format(nrow(final_output_table), big.mark=","), "unique rows."))

# Show summary of results
message("\n=== FINAL RESULTS SUMMARY ===")
message(sprintf("Filtering applied: %s <= %.3f", PVALUE_COLUMN_FOR_FILTERING, PVALUE_THRESHOLD))
message(sprintf("Unique pathways: %d", uniqueN(final_output_table$Pathway_ID)))
message(sprintf("Unique genes: %d", uniqueN(final_output_table$eQTL_Gene_Entrez)))
message(sprintf("Unique variants: %d", uniqueN(final_output_table$SNP_chr_pos_ref_alt)))
if (sum(!is.na(final_output_table$SNP_rsID)) > 0) {
  message(sprintf("Variants with rsIDs: %d (%.1f%%)", 
                  sum(!is.na(final_output_table$SNP_rsID)),
                  100*sum(!is.na(final_output_table$SNP_rsID))/nrow(final_output_table)))
}

message("First few rows of the final output table (ordered by significance):")
# Order by the p-value column that was used for filtering
if (paste0("eQTL_", PVALUE_COLUMN_FOR_FILTERING) %in% colnames(final_output_table)) {
  order_col <- paste0("eQTL_", PVALUE_COLUMN_FOR_FILTERING)
  print(head(final_output_table[order(get(order_col))], 10))
} else {
  print(head(final_output_table[order(!is.na(SNP_rsID), eQTL_pval_perm)], 10))
}

# --- VI. Save Output ---
message(paste("\nSaving final table to CSV:", OUTPUT_CSV_FILE))
fwrite(final_output_table, OUTPUT_CSV_FILE)
message(paste("Saving final table to RDS:", OUTPUT_RDS_FILE))
write_rds(final_output_table, OUTPUT_RDS_FILE)

message("\n--- Analysis Complete (GTEx v10 workflow with direct rsID from file) ---")