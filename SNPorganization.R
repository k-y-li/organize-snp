# --- 0. Load Libraries ---
# Ensure these are installed:
# install.packages(c("data.table", "dplyr", "readr", "BiocManager"))
# BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db", "KEGGREST")) # SNPlocs not needed for this version

library(data.table)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(KEGGREST)
library(readr)

# --- I. Configuration & File Paths ---
# *** Path to your GTEx "significant pairs" file (ASSUMED to be hg38) ***
LOCAL_GTEX_EQTL_FILE <- "/Users/Kevin/gtex data/GTEx_Analysis_v8_eQTL/Brain_Cortex.v8.signif_variant_gene_pairs.txt.gz"

# *** CRITICAL: VERIFY this prefix by looking at your file's variant_id column ***
# If it's an hg38 file, it's likely "chr22_". If it's an hg19 file (e.g. standard V8 signif pairs), it's "22_".
# For this script to work as an "hg38 workflow", this prefix should match hg38 variant_ids.
TARGET_CHR_PREFIX_IN_FILE <- "chr22_" # <<< ADJUST THIS BASED ON YOUR FILE INSPECTION

# Significance column for reporting (from the significant pairs file)
SIGNIFICANCE_COL_FOR_REPORTING <- "pval_nominal" # Common in GTEx signif_pairs files

# Cache file for KEGG pathway data (ensure it's an hg38-compatible cache if you have one)
KEGG_DATA_RDS_FILE <- "kegg_hsa_pathways_via_bioc_hg38.rds" 

# Output files (reflecting "signif" and hg38, no_rsid)
OUTPUT_CSV_FILE <- "chr22_hg38_signif_eqtls_kegg_pathways_no_rsid.csv"
OUTPUT_RDS_FILE <- "chr22_hg38_signif_eqtls_kegg_pathways_no_rsid.rds"

# --- II. Process and Filter GTEx "significant pairs" eQTL Data ---
message("--- Part 1: Processing GTEx 'significant pairs' eQTL data (ASSUMED hg38) ---")

if (!file.exists(LOCAL_GTEX_EQTL_FILE)) {
  stop(paste("'Significant pairs' file not found at:", LOCAL_GTEX_EQTL_FILE))
}
message(paste("Found local 'significant pairs' file:", LOCAL_GTEX_EQTL_FILE))
message(paste("File size:", round(file.info(LOCAL_GTEX_EQTL_FILE)$size / (1024^2), 2), "MB."))

gtex_final_filtered_eqtls <- NULL 

tryCatch({
  message("Reading header of 'significant pairs' file to check columns and variant_id format...")
  header_check_dt <- fread(LOCAL_GTEX_EQTL_FILE, nrows = 5)
  all_columns <- colnames(header_check_dt)
  message("Available columns: ", paste(all_columns, collapse=", "))
  message("First few variant_ids to check format (e.g., for 'chr' prefix and build suffix like _b38):")
  print(header_check_dt$variant_id)
  
  if (nrow(header_check_dt) > 0 && !any(grepl("_b38|GRCh38", header_check_dt$variant_id, ignore.case = TRUE))) {
    warning(paste("Variant IDs in", basename(LOCAL_GTEX_EQTL_FILE), 
                  "do not clearly indicate hg38/GRCh38 build in the first few entries.",
                  "\nThe rest of the script assumes hg38. Please verify genome build MANUALLY."))
  }
  rm(header_check_dt); gc()
  
  # Define columns to select.
  cols_to_select <- c("variant_id", "gene_id")
  if (SIGNIFICANCE_COL_FOR_REPORTING %in% all_columns) {
    cols_to_select <- c(cols_to_select, SIGNIFICANCE_COL_FOR_REPORTING)
  } else {
    warning(paste("Specified significance column '", SIGNIFICANCE_COL_FOR_REPORTING, "' not found. P-value may be missing in output.", sep=""))
  }
  if ("slope" %in% all_columns) cols_to_select <- c(cols_to_select, "slope")
  if ("tss_distance" %in% all_columns) cols_to_select <- c(cols_to_select, "tss_distance")
  cols_to_select <- unique(cols_to_select)
  
  message(paste("\nReading selected columns (", paste(cols_to_select, collapse=", "), ") from 'significant pairs' file..."))
  full_data_dt <- fread(LOCAL_GTEX_EQTL_FILE, select = cols_to_select, showProgress = TRUE)
  message(paste("Read", format(nrow(full_data_dt), big.mark=","), "rows (all significant by definition)."))
  
  # Optional: Filter by TSS distance if desired (though signif_pairs are already cis)
  # MAX_TSS_DISTANCE <- 1000000
  # if ("tss_distance" %in% colnames(full_data_dt)) {
  #     full_data_dt <- full_data_dt[abs(tss_distance) <= MAX_TSS_DISTANCE]
  #     message(paste("Rows after TSS distance filter: ", format(nrow(full_data_dt), big.mark=",")))
  # }
  
  message(paste("Filtering for chromosome 22 (prefix '", TARGET_CHR_PREFIX_IN_FILE, "')...", sep=""))
  gtex_final_filtered_eqtls <- full_data_dt[startsWith(variant_id, TARGET_CHR_PREFIX_IN_FILE)]
  message(paste(format(nrow(gtex_final_filtered_eqtls), big.mark=","), "rows after chromosome filter."))
  rm(full_data_dt); gc()
  
  # NO ADDITIONAL SIGNIFICANCE FILTERING, as file is already significant pairs.
  
  if (nrow(gtex_final_filtered_eqtls) == 0) {
    stop("No eQTLs found for chromosome 22 after filtering. Check chromosome prefix.")
  }
  message("Filtered eQTL data (gtex_final_filtered_eqtls) prepared:")
  print(head(gtex_final_filtered_eqtls))
  
}, error = function(e) {
  stop(paste("Error during eQTL data processing: ", e$message))
})


# --- III. Map eGene Ensembl IDs to Entrez IDs ---
# This section is the same as your previous hg38 script
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
  gtex_eqtls_for_kegg <- gtex_eqtls_with_entrez[!is.na(entrez_gene_id)]
  message(paste(format(nrow(gtex_eqtls_for_kegg), big.mark=","), "eQTLs remaining for KEGG analysis."))
  if (nrow(gtex_eqtls_for_kegg) == 0) stop("No eQTLs remaining after Entrez ID mapping.")
}, error = function(e) { stop(paste("Error during Ensembl to Entrez ID mapping: ", e$message)) })


# --- IV. Get KEGG Pathway Data (using org.Hs.eg.db and KEGGREST for names) ---
# This section is the same as your previous hg38 script
message("\n--- Part 3: Fetching KEGG Pathway Data ---")
kegg_pathways_df <- NULL
if (file.exists(KEGG_DATA_RDS_FILE)) {
  message("Loading cached KEGG pathway data from ", KEGG_DATA_RDS_FILE)
  kegg_pathways_df <- read_rds(KEGG_DATA_RDS_FILE)
} else {
  message("Constructing KEGG pathway-gene relationships using org.Hs.eg.db...")
  kegg_mappings <- tryCatch(select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = "ENTREZID"), columns = c("PATH", "SYMBOL"), keytype = "ENTREZID"), error = function(e) { message("Error selecting from org.Hs.eg.db: ", e$message); NULL })
  if (is.null(kegg_mappings) || nrow(kegg_mappings) == 0) stop("Could not get KEGG mappings from org.Hs.eg.db.")
  kegg_mappings_filtered <- kegg_mappings[!is.na(kegg_mappings$PATH), ]
  if (nrow(kegg_mappings_filtered) == 0) stop("No genes with KEGG pathway annotations in org.Hs.eg.db.")
  kegg_pathways_dt <- data.table(pathway_id_num = kegg_mappings_filtered$PATH, entrez_gene_id_kegg = as.character(kegg_mappings_filtered$ENTREZID), gene_symbol_from_bioc = kegg_mappings_filtered$SYMBOL)
  kegg_pathways_dt[, pathway_id := paste0("hsa", pathway_id_num)][, pathway_id_num := NULL]
  message(paste("Found", format(nrow(kegg_pathways_dt), big.mark=","), "gene-pathway associations from org.Hs.eg.db."))
  message("Fetching all KEGG pathway names for 'hsa'...")
  all_kegg_names_info <- tryCatch(keggList("pathway", "hsa"), error = function(e) { message("WARNING: Could not fetch all KEGG pathway names via API: ", e$message); NULL })
  pathway_names_dt <- data.table()
  if (!is.null(all_kegg_names_info) && length(all_kegg_names_info) > 0) {
    pathway_names_dt <- data.table(pathway_id_full_kegg = names(all_kegg_names_info), pathway_name = as.character(all_kegg_names_info))
    pathway_names_dt[, pathway_id := sub("path:", "", pathway_id_full_kegg)][, pathway_id_full_kegg := NULL]
    message(paste("Fetched names for", nrow(pathway_names_dt), "KEGG pathways."))
  } else { 
    unique_ids_from_data <- unique(kegg_pathways_dt$pathway_id)
    pathway_names_dt <- data.table(pathway_id = unique_ids_from_data, pathway_name = paste("Name API call failed for", unique_ids_from_data))
    message("WARNING: Using placeholder pathway names as API call for all names failed.")
  }
  kegg_pathways_df <- merge(kegg_pathways_dt, pathway_names_dt, by = "pathway_id", all.x = TRUE)
  kegg_pathways_df[is.na(pathway_name), pathway_name := paste("Name not in API list for", pathway_id)]
  kegg_pathways_df <- kegg_pathways_df[, .(pathway_id, pathway_name, entrez_gene_id_kegg, gene_symbol_kegg = gene_symbol_from_bioc)]
  write_rds(kegg_pathways_df, KEGG_DATA_RDS_FILE)
  message("Saved KEGG pathway data to ", KEGG_DATA_RDS_FILE)
}
if (is.null(kegg_pathways_df) || nrow(kegg_pathways_df) == 0) stop("KEGG pathway data frame is empty.")
message(paste("Loaded/Constructed", uniqueN(kegg_pathways_df$pathway_id), "KEGG pathways."))


# --- V. Merge eQTLs with KEGG Pathways and Prepare Final Output ---
message("\n--- Part 4: Merging eQTLs with KEGG Pathways & Finalizing ---")
if (is.null(gtex_eqtls_for_kegg) || nrow(gtex_eqtls_for_kegg) == 0) stop("eQTL data for KEGG ('gtex_eqtls_for_kegg') is empty.")
gtex_eqtls_for_kegg[, entrez_gene_id := as.character(entrez_gene_id)]
merged_data <- merge(gtex_eqtls_for_kegg, kegg_pathways_df, by.x = "entrez_gene_id", by.y = "entrez_gene_id_kegg", allow.cartesian = TRUE)
if (nrow(merged_data) == 0) stop("No overlap found between eGenes and KEGG pathway genes.")
message(paste("Found", format(nrow(merged_data), big.mark=","), "eQTL-gene-pathway associations."))

# Get the p-value column for reporting.
# Since this is a "significant pairs" file, 'pval_nominal' is usually what's directly available and reported.
# If your file has a different primary significance metric (e.g. qval), adjust SIGNIFICANCE_COL_FOR_REPORTING at the top.
p_value_col_to_report <- SIGNIFICANCE_COL_FOR_REPORTING 
if (!p_value_col_to_report %in% colnames(merged_data)) {
  warning(paste0("Expected p-value column '", p_value_col_to_report, "' not found in merged data. P-value in output may be incorrect or missing."))
  # Try to find any pval column if the specific one is missing
  available_pval_cols <- colnames(merged_data)[grepl("pval|qval", colnames(merged_data), ignore.case=TRUE)]
  if(length(available_pval_cols) > 0) {
    p_value_col_to_report <- available_pval_cols[1]
    message(paste0("Using first available p-value like column for reporting: ", p_value_col_to_report))
  } else {
    p_value_col_to_report <- NA_character_ # Will result in NA column
    message("No p-value like column found for reporting.")
  }
}


final_output_table <- merged_data[, .(
  Pathway_ID = pathway_id,
  Pathway_Name = pathway_name,
  Tissue = "Brain_Cortex_hg38", # ASSUMING hg38 data source
  SNP_chr_pos_ref_alt = variant_id, 
  # SNP_rsID = rsID, # rsID column removed as per request
  eQTL_Gene_Ensembl_Versioned = gene_id,
  eQTL_Gene_Ensembl_NoVersion = ensembl_id_no_version,
  eQTL_Gene_Entrez = entrez_gene_id,
  eQTL_Gene_Symbol = gene_symbol_kegg,
  eQTL_PValue_Reported = if(!is.na(p_value_col_to_report) && p_value_col_to_report %in% colnames(merged_data)) get(p_value_col_to_report) else NA_real_,
  eQTL_Slope = if("slope" %in% colnames(merged_data)) slope else NA_real_, 
  eQTL_TSS_Distance = if("tss_distance" %in% colnames(merged_data)) tss_distance else NA_integer_
)]
if(!is.na(p_value_col_to_report) && p_value_col_to_report %in% colnames(merged_data)) {
  setnames(final_output_table, "eQTL_PValue_Reported", paste0("eQTL_", p_value_col_to_report))
}


final_output_table <- unique(final_output_table)
message(paste("\nFinal table has", format(nrow(final_output_table), big.mark=","), "unique rows."))
message("First few rows of the final output table:")
print(head(final_output_table))


# --- VI. Save Output ---
message(paste("\nSaving final table to CSV:", OUTPUT_CSV_FILE))
fwrite(final_output_table, OUTPUT_CSV_FILE)
message(paste("Saving final table to RDS:", OUTPUT_RDS_FILE))
write_rds(final_output_table, OUTPUT_RDS_FILE)

message("\n--- Analysis Complete (hg38 significant pairs workflow) ---")