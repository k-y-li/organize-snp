# --- 0. Load Libraries ---
# Ensure these are installed:
# install.packages(c("data.table", "dplyr", "readr", "BiocManager"))
# BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db", "KEGGREST", "SNPlocs.Hsapiens.dbSNP155.GRCh38", "GenomicRanges"))

library(data.table)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(KEGGREST)
library(readr)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38) # For hg38 rsID lookup
library(GenomicRanges)                   # Used with SNPlocs

# --- I. Configuration & File Paths ---
# Path to your GTEx "significant pairs" file (ASSUMED to be hg38)
LOCAL_GTEX_EQTL_FILE <- "/Users/Kevin/gtex data/GTEx_Analysis_v8_eQTL/Brain_Cortex.v8.signif_variant_gene_pairs.txt.gz"

# CRITICAL: VERIFY this prefix by looking at your file's variant_id column for hg38 data.
# If variant_ids are like "22_...", change this to "22_".
TARGET_CHR_PREFIX_FOR_EQTL_FILTER <- "chr22_" 

# Optional TSS distance filter (GTEx significant pairs are already cis, usually within 1Mb)
MAX_TSS_DISTANCE <- 1000000 

# Significance column for reporting (from the significant pairs file)
SIGNIFICANCE_COL_FOR_REPORTING <- "pval_nominal" 

# Cache file for KEGG pathway data (hg38 compatible)
KEGG_DATA_RDS_FILE <- "kegg_hsa_pathways_via_bioc_hg38.rds" # Reverted to a clear hg38 name

# Output files (reverting to names used with the "allpairs" hg38 script)
OUTPUT_CSV_FILE <- "chr22_hg38_brain_eqtls_kegg_pathways_final.csv"
OUTPUT_RDS_FILE <- "chr22_hg38_brain_eqtls_kegg_pathways_final.rds"

# --- II. Process and Filter GTEx "significant pairs" eQTL Data (hg38) ---
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
  
  if ("tss_distance" %in% colnames(full_data_dt)) {
    message(paste("Original rows before TSS distance filter:", format(nrow(full_data_dt), big.mark=",")))
    full_data_dt <- full_data_dt[abs(tss_distance) <= MAX_TSS_DISTANCE]
    message(paste("Rows after TSS distance filter (abs(tss_distance) <= ", format(MAX_TSS_DISTANCE, big.mark=","), "): ", format(nrow(full_data_dt), big.mark=","), sep=""))
    if (nrow(full_data_dt) == 0) stop("No eQTLs remaining after TSS distance filter.")
  }
  
  message(paste("Filtering for chromosome 22 (prefix '", TARGET_CHR_PREFIX_IN_FILE, "')...", sep=""))
  gtex_final_filtered_eqtls <- full_data_dt[startsWith(variant_id, TARGET_CHR_PREFIX_IN_FILE)]
  message(paste(format(nrow(gtex_final_filtered_eqtls), big.mark=","), "rows after chromosome filter."))
  rm(full_data_dt); gc()
  
  if (nrow(gtex_final_filtered_eqtls) == 0) {
    stop("No eQTLs found for chromosome 22 after filtering. Check chromosome prefix.")
  }
  message("Filtered eQTL data (gtex_final_filtered_eqtls) prepared:")
  print(head(gtex_final_filtered_eqtls))
  
}, error = function(e) {
  stop(paste("Error during eQTL data processing: ", e$message))
})


# --- III. Map Variant IDs to rsIDs (using hg38 SNPlocs) ---
message("\n--- Part 2: Mapping Variant IDs to rsIDs (hg38/GRCh38) ---")
if (!exists("gtex_final_filtered_eqtls") || is.null(gtex_final_filtered_eqtls) || nrow(gtex_final_filtered_eqtls) == 0) {
  message("gtex_final_filtered_eqtls is empty or does not exist. Skipping rsID mapping.")
  if (exists("gtex_final_filtered_eqtls") && is.data.table(gtex_final_filtered_eqtls)) { # Ensure column exists if table is empty
    gtex_final_filtered_eqtls[, rsID := NA_character_]
  }
} else {
  unique_snps_to_map <- unique(gtex_final_filtered_eqtls[, .(variant_id)])
  message(paste("Found", nrow(unique_snps_to_map), "unique variant_ids to map to rsIDs."))
  
  unique_snps_to_map[, c("chr_col_raw", "pos_str", "ref", "alt", "build_str") := tstrsplit(variant_id, "_", fixed=TRUE)]
  unique_snps_to_map[, pos := as.integer(pos_str)]
  
  snpdb_obj <- SNPlocs.Hsapiens.dbSNP155.GRCh38 
  current_snplocs_style <- seqlevelsStyle(snpdb_obj)[1]
  message("SNPlocs object seqlevelsStyle: ", current_snplocs_style)
  
  unique_snps_to_map[, chr_for_snplocs := chr_col_raw] 
  if (current_snplocs_style == "UCSC") { 
    if (!all(startsWith(unique_snps_to_map$chr_col_raw, "chr"))) {
      message("Adjusting chromosome names from eQTL file to 'chrX' (UCSC) style for SNPlocs...")
      unique_snps_to_map[ !startsWith(chr_col_raw, "chr") & chr_col_raw %in% as.character(c(1:22,"X","Y","M","MT")), 
                          chr_for_snplocs := paste0("chr", chr_col_raw) ]
    }
  } else if (current_snplocs_style == "NCBI") { 
    if (any(startsWith(unique_snps_to_map$chr_col_raw, "chr"))) { 
      message("Adjusting chromosome names from eQTL file from 'chrX' to 'X' (NCBI) style for SNPlocs...")
      unique_snps_to_map[, chr_for_snplocs := sub("^chr", "", chr_col_raw)]
    }
  } else {
    message(paste0("Warning: SNPlocs style is '", current_snplocs_style, "'. Data chr style: '", 
                   substr(unique_snps_to_map$chr_col_raw[1],1,3),"'. Check adjustment."))
  }
  
  message("Sample of parsed SNP components for rsID mapping (using 'chr_for_snplocs'):")
  print(head(unique_snps_to_map[, .(variant_id, chr_col_raw, chr_for_snplocs, pos)]))
  
  target_chromosomes_for_snplocs <- unique(unique_snps_to_map$chr_for_snplocs)
  message("Attempting to fetch SNP locations from SNPlocs for chromosomes: ", paste(target_chromosomes_for_snplocs, collapse=", "))
  
  all_snplocs_for_target_chroms_list <- list()
  for(chrom_name_to_query in target_chromosomes_for_snplocs){
    if(chrom_name_to_query %in% seqlevels(snpdb_obj)){
      message(paste("  Querying SNPlocs for chromosome:", chrom_name_to_query))
      snps_on_chrom_gr <- snpsBySeqname(snpdb_obj, chrom_name_to_query) # Removed as.GRanges=TRUE
      if(length(snps_on_chrom_gr) > 0){
        all_snplocs_for_target_chroms_list[[chrom_name_to_query]] <- data.table(
          chr_for_snplocs_db = chrom_name_to_query, 
          pos_snploc_db = start(snps_on_chrom_gr),
          rsID = mcols(snps_on_chrom_gr)$RefSNP_id
        )
        message(paste("    Retrieved", format(length(snps_on_chrom_gr), big.mark=","), "SNPs from SNPlocs for chromosome", chrom_name_to_query))
      } else { message(paste("    No SNPs found in SNPlocs for chromosome:", chrom_name_to_query)) }
    } else {
      message(paste("    Chromosome '", chrom_name_to_query, "' not found in SNPlocs seqlevels. Skipping."))
    }
  }
  
  if(length(all_snplocs_for_target_chroms_list) > 0) {
    all_snplocs_dt <- rbindlist(all_snplocs_for_target_chroms_list)
    setnames(all_snplocs_dt, c("chr_for_snplocs_db", "pos_snploc_db"), c("chr_for_snplocs", "pos")) 
    all_snplocs_dt_unique_rsid <- all_snplocs_dt[!is.na(rsID), .(rsID = first(rsID)), by = .(chr_for_snplocs, pos)]
    
    setkey(unique_snps_to_map, chr_for_snplocs, pos)
    setkey(all_snplocs_dt_unique_rsid, chr_for_snplocs, pos)
    unique_snps_with_rsids <- merge(unique_snps_to_map, all_snplocs_dt_unique_rsid, all.x = TRUE)
    
    message(paste(sum(!is.na(unique_snps_with_rsids$rsID)), "out of", nrow(unique_snps_with_rsids), "unique SNPs mapped to an rsID."))
    
    gtex_final_filtered_eqtls <- merge(
      gtex_final_filtered_eqtls,
      unique_snps_with_rsids[, .(variant_id, rsID)], 
      by = "variant_id", all.x = TRUE
    )
    message(paste("Mapped rsIDs to", format(sum(!is.na(gtex_final_filtered_eqtls$rsID)), big.mark=","), "eQTL entries."))
    if(sum(!is.na(gtex_final_filtered_eqtls$rsID)) > 0) print(head(gtex_final_filtered_eqtls[!is.na(rsID)]))
  } else {
    message("Could not retrieve any SNP locations from SNPlocs. rsID column will be all NA.")
    gtex_final_filtered_eqtls[, rsID := NA_character_]
  }
}


# --- IV. Map eGene Ensembl IDs to Entrez IDs ---
message("\n--- Part 3: Mapping eGene Ensembl IDs to Entrez IDs ---")
# (This section remains the same as your last working version)
# ... (Paste your working Ensembl to Entrez mapping code here, ensuring input is gtex_final_filtered_eqtls)
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
  gtex_eqtls_for_kegg <- gtex_eqtls_with_entrez[!is.na(entrez_gene_id)] # This now contains rsID if mapped
  message(paste(format(nrow(gtex_eqtls_for_kegg), big.mark=","), "eQTLs remaining for KEGG analysis."))
  if (nrow(gtex_eqtls_for_kegg) == 0) stop("No eQTLs remaining after Entrez ID mapping.")
}, error = function(e) { stop(paste("Error during Ensembl to Entrez ID mapping: ", e$message)) })


# --- V. Get KEGG Pathway Data (using org.Hs.eg.db and KEGGREST for names) ---
message("\n--- Part 4: Fetching KEGG Pathway Data ---")
# (This section remains the same as your last working version, using KEGG_DATA_RDS_FILE)
# ... (Paste your working KEGG pathway fetching code here) ...
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


# --- VI. Merge eQTLs with KEGG Pathways and Prepare Final Output ---
message("\n--- Part 5: Merging eQTLs with KEGG Pathways & Finalizing ---")
if (is.null(gtex_eqtls_for_kegg) || nrow(gtex_eqtls_for_kegg) == 0) stop("eQTL data for KEGG ('gtex_eqtls_for_kegg') is empty.")
gtex_eqtls_for_kegg[, entrez_gene_id := as.character(entrez_gene_id)]
merged_data <- merge(gtex_eqtls_for_kegg, kegg_pathways_df, by.x = "entrez_gene_id", by.y = "entrez_gene_id_kegg", allow.cartesian = TRUE)
if (nrow(merged_data) == 0) stop("No overlap found between eGenes and KEGG pathway genes.")
message(paste("Found", format(nrow(merged_data), big.mark=","), "eQTL-gene-pathway associations."))

p_value_col_to_report <- SIGNIFICANCE_COL_FOR_REPORTING 
if (!p_value_col_to_report %in% colnames(merged_data)) {
  warning(paste0("Expected p-value column '", p_value_col_to_report, "' not found in merged data. P-value in output may be incorrect."))
  available_pval_cols <- colnames(merged_data)[grepl("pval|qval", colnames(merged_data), ignore.case=TRUE)]
  if(length(available_pval_cols) > 0) p_value_col_to_report <- available_pval_cols[1] else p_value_col_to_report <- NA_character_
}

final_output_table <- merged_data[, .(
  Pathway_ID = pathway_id,
  Pathway_Name = pathway_name,
  Tissue = "Brain_Cortex_hg38", # ASSUMING hg38 data source
  SNP_chr_pos_ref_alt = variant_id, 
  SNP_rsID = rsID, # rsID column is now included         
  eQTL_Gene_Ensembl_Versioned = gene_id,
  eQTL_Gene_Ensembl_NoVersion = ensembl_id_no_version,
  eQTL_Gene_Entrez = entrez_gene_id,
  eQTL_Gene_Symbol = gene_symbol_kegg,
  eQTL_PValue_Reported = if(!is.na(p_value_col_to_report) && p_value_col_to_report %in% colnames(.SD)) get(p_value_col_to_report) else NA_real_,
  eQTL_Slope = if("slope" %in% colnames(.SD)) slope else NA_real_, 
  eQTL_TSS_Distance = if("tss_distance" %in% colnames(.SD)) tss_distance else NA_integer_ 
)] # .SD refers to Subset of Data.table, good for dynamic column access within data.table 'j' expression.
if(!is.na(p_value_col_to_report) && p_value_col_to_report %in% colnames(final_output_table)) { # Check if column exists before renaming
  setnames(final_output_table, "eQTL_PValue_Reported", paste0("eQTL_", p_value_col_to_report))
} else if ("eQTL_PValue_Reported" %in% colnames(final_output_table) && is.na(p_value_col_to_report)){
  setnames(final_output_table, "eQTL_PValue_Reported", "eQTL_PValue_Unavailable")
}


final_output_table <- unique(final_output_table)
message(paste("\nFinal table has", format(nrow(final_output_table), big.mark=","), "unique rows."))
message("First few rows of the final output table (check for SNP_rsID):")
print(head(final_output_table[order(!is.na(SNP_rsID), SNP_rsID)], 20))

# --- VII. Save Output ---
message(paste("\nSaving final table to CSV:", OUTPUT_CSV_FILE))
fwrite(final_output_table, OUTPUT_CSV_FILE)
message(paste("Saving final table to RDS:", OUTPUT_RDS_FILE))
write_rds(final_output_table, OUTPUT_RDS_FILE)

message("\n--- Analysis Complete (hg38 significant pairs workflow with rsID mapping) ---")