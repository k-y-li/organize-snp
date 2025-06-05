# --- 0. Load Libraries ---
# Ensure these are installed:
# install.packages(c("data.table", "dplyr", "readr", "BiocManager", "openxlsx"))
# BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db", "KEGGREST"))

library(data.table)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(KEGGREST)
library(readr)
library(openxlsx)

# --- I. Configuration & File Paths ---
# Path to your GTEx v10 eGenes file (GRCh38/hg38)
LOCAL_GTEX_EQTL_FILE <- "/Users/Kevin/gtex data/GTEx_Analysis_v10_eQTL_updated/Brain_Cortex.v10.eGenes.txt.gz"

# Chromosomes to analyze (autosomal only, excluding X, Y, M/MT)
TARGET_CHROMOSOMES <- paste0("chr", 1:22)

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
OUTPUT_EXCEL_FILE <- "all_chromosomes_v10_brain_eqtls_kegg_pathways.xlsx"
OUTPUT_RDS_FILE <- "all_chromosomes_v10_brain_eqtls_kegg_pathways.rds"

# --- II. Load and Prepare GTEx Data ---
message("--- Part 1: Loading GTEx v10 eGenes data (GRCh38/hg38) ---")

if (!file.exists(LOCAL_GTEX_EQTL_FILE)) {
  stop(paste("GTEx v10 eGenes file not found at:", LOCAL_GTEX_EQTL_FILE))
}
message(paste("Found GTEx v10 eGenes file:", LOCAL_GTEX_EQTL_FILE))
message(paste("File size:", round(file.info(LOCAL_GTEX_EQTL_FILE)$size / (1024^2), 2), "MB."))

# Load the full dataset once
message("Reading full GTEx v10 eGenes file...")
tryCatch({
  # Check columns first
  header_check_dt <- fread(LOCAL_GTEX_EQTL_FILE, nrows = 5)
  all_columns <- colnames(header_check_dt)
  message("Available columns: ", paste(all_columns, collapse=", "))
  
  # Check if rsID column exists
  if (RSID_COLUMN_NAME %in% all_columns) {
    message(paste("Found rsID column:", RSID_COLUMN_NAME))
  } else {
    warning(paste("rsID column", RSID_COLUMN_NAME, "not found. Will set rsIDs to NA."))
    RSID_COLUMN_NAME <- NULL
  }
  
  rm(header_check_dt); gc()
  
  # Select columns to read
  cols_to_select <- c("gene_id", "gene_name", "variant_id", "tss_distance", 
                      "pval_nominal", "pval_perm", "pval_beta", "qval", 
                      "slope", "slope_se", "chr", "variant_pos", "ref", "alt")
  
  if (!is.null(RSID_COLUMN_NAME)) {
    cols_to_select <- c(cols_to_select, RSID_COLUMN_NAME)
  }
  
  cols_to_select <- cols_to_select[cols_to_select %in% all_columns]
  
  # Read full dataset
  message("Loading complete dataset...")
  full_gtex_data <- fread(LOCAL_GTEX_EQTL_FILE, select = cols_to_select, showProgress = TRUE)
  message(paste("Loaded", format(nrow(full_gtex_data), big.mark=","), "total eQTL associations."))
  
  # Rename rsID column if it exists
  if (!is.null(RSID_COLUMN_NAME) && RSID_COLUMN_NAME %in% colnames(full_gtex_data)) {
    setnames(full_gtex_data, RSID_COLUMN_NAME, "rsID")
  } else {
    full_gtex_data[, rsID := NA_character_]
  }
  
  # Apply filters that apply to all chromosomes
  if ("tss_distance" %in% colnames(full_gtex_data)) {
    message("Applying TSS distance filter...")
    original_rows <- nrow(full_gtex_data)
    full_gtex_data <- full_gtex_data[abs(tss_distance) <= MAX_TSS_DISTANCE]
    message(sprintf("TSS filter: %s -> %s rows", format(original_rows, big.mark=","), 
                    format(nrow(full_gtex_data), big.mark=",")))
  }
  
  if (PVALUE_COLUMN_FOR_FILTERING %in% colnames(full_gtex_data)) {
    message("Applying p-value significance filter...")
    original_rows <- nrow(full_gtex_data)
    full_gtex_data <- full_gtex_data[get(PVALUE_COLUMN_FOR_FILTERING) <= PVALUE_THRESHOLD & 
                                       !is.na(get(PVALUE_COLUMN_FOR_FILTERING))]
    message(sprintf("P-value filter (%s <= %.3f): %s -> %s rows", 
                    PVALUE_COLUMN_FOR_FILTERING, PVALUE_THRESHOLD,
                    format(original_rows, big.mark=","), format(nrow(full_gtex_data), big.mark=",")))
  }
  
  if (nrow(full_gtex_data) == 0) {
    stop("No eQTLs remaining after filtering. Consider relaxing filters.")
  }
  
}, error = function(e) {
  stop(paste("Error loading GTEx data: ", e$message))
})

# --- III. Load KEGG Pathway Data Once ---
message("\n--- Part 2: Loading KEGG Pathway Data ---")

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
message(paste("Loaded", uniqueN(kegg_pathways_df$pathway_id), "KEGG pathways."))

# --- IV. Function to Process One Chromosome ---
process_chromosome <- function(chr_name, gtex_data, kegg_data) {
  message(sprintf("\n--- Processing %s ---", chr_name))
  
  # Filter for this chromosome
  chr_prefix <- paste0(chr_name, "_")
  chr_eqtls <- gtex_data[startsWith(variant_id, chr_prefix)]
  
  if (nrow(chr_eqtls) == 0) {
    message(sprintf("No eQTLs found for %s", chr_name))
    return(NULL)
  }
  
  message(sprintf("Found %s eQTLs for %s", format(nrow(chr_eqtls), big.mark=","), chr_name))
  
  # Map Ensembl IDs to Entrez IDs
  unique_ensembl_ids <- unique(sub("\\..*$", "", chr_eqtls$gene_id))
  message(sprintf("Mapping %d unique genes to Entrez IDs...", length(unique_ensembl_ids)))
  
  entrez_map <- mapIds(org.Hs.eg.db, keys = unique_ensembl_ids,
                       column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  ensembl_to_entrez_dt <- data.table(ensembl_id_no_version = names(entrez_map[!is.na(entrez_map)]),
                                     entrez_gene_id = entrez_map[!is.na(entrez_map)])
  
  if (nrow(ensembl_to_entrez_dt) == 0) {
    message(sprintf("No genes mapped to Entrez IDs for %s", chr_name))
    return(NULL)
  }
  
  message(sprintf("Mapped %d genes to Entrez IDs", nrow(ensembl_to_entrez_dt)))
  
  # Add Entrez IDs to eQTL data
  chr_eqtls[, ensembl_id_no_version := sub("\\..*$", "", gene_id)]
  chr_eqtls_with_entrez <- merge(chr_eqtls, ensembl_to_entrez_dt, 
                                 by = "ensembl_id_no_version", all.x = TRUE)
  chr_eqtls_for_kegg <- chr_eqtls_with_entrez[!is.na(entrez_gene_id)]
  
  if (nrow(chr_eqtls_for_kegg) == 0) {
    message(sprintf("No eQTLs remaining after Entrez mapping for %s", chr_name))
    return(NULL)
  }
  
  # Merge with KEGG pathways
  chr_eqtls_for_kegg[, entrez_gene_id := as.character(entrez_gene_id)]
  merged_data <- merge(chr_eqtls_for_kegg, kegg_data, 
                       by.x = "entrez_gene_id", by.y = "entrez_gene_id_kegg", allow.cartesian = TRUE)
  
  if (nrow(merged_data) == 0) {
    message(sprintf("No KEGG pathway matches for %s", chr_name))
    return(NULL)
  }
  
  message(sprintf("Found %s eQTL-pathway associations for %s", 
                  format(nrow(merged_data), big.mark=","), chr_name))
  
  # Create final output for this chromosome
  chr_output <- merged_data[, .(
    Chromosome = chr_name,
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
  
  chr_output <- unique(chr_output)
  
  # Summary statistics
  message(sprintf("Final %s results: %s unique associations, %d pathways, %d genes, %d variants", 
                  chr_name, format(nrow(chr_output), big.mark=","), 
                  uniqueN(chr_output$Pathway_ID),
                  uniqueN(chr_output$eQTL_Gene_Entrez),
                  uniqueN(chr_output$SNP_chr_pos_ref_alt)))
  
  return(chr_output)
}

# --- V. Process All Chromosomes ---
message("\n--- Part 3: Processing All Chromosomes ---")

all_chromosome_results <- list()
summary_stats <- data.table()

for (chr in TARGET_CHROMOSOMES) {
  chr_result <- process_chromosome(chr, full_gtex_data, kegg_pathways_df)
  
  if (!is.null(chr_result)) {
    all_chromosome_results[[chr]] <- chr_result
    
    # Collect summary statistics
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
  
  # Clean up memory periodically
  if (length(all_chromosome_results) %% 5 == 0) {
    gc()
  }
}

# --- VI. Save Results to Excel with Multiple Tabs ---
message("\n--- Part 4: Creating Excel Output ---")

if (length(all_chromosome_results) == 0) {
  stop("No chromosome results to save!")
}

# Create Excel workbook
wb <- createWorkbook()

# Add summary sheet first
message("Creating Summary sheet...")
addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_stats, startRow = 1, startCol = 1, rowNames = FALSE)

# Format summary sheet
summary_style <- createStyle(fontSize = 11, fontColour = "black", halign = "center", 
                             fgFill = "grey", border = "TopBottomLeftRight")
addStyle(wb, "Summary", summary_style, rows = 1, cols = 1:ncol(summary_stats), gridExpand = TRUE)
setColWidths(wb, "Summary", cols = 1:ncol(summary_stats), widths = "auto")

# Add each chromosome as a separate sheet
for (chr in names(all_chromosome_results)) {
  chr_data <- all_chromosome_results[[chr]]
  
  # Create worksheet name (Excel has 31 character limit)
  sheet_name <- gsub("chr", "Chr", chr)
  message(sprintf("Adding sheet: %s (%s associations)", sheet_name, format(nrow(chr_data), big.mark=",")))
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, chr_data, startRow = 1, startCol = 1, rowNames = FALSE)
  
  # Format headers
  header_style <- createStyle(fontSize = 10, fontColour = "black", halign = "center", 
                              fgFill = "grey", border = "TopBottomLeftRight", textDecoration = "bold")
  addStyle(wb, sheet_name, header_style, rows = 1, cols = 1:ncol(chr_data), gridExpand = TRUE)
  
  # Auto-size columns
  setColWidths(wb, sheet_name, cols = 1:ncol(chr_data), widths = "auto")
}

# Save Excel file
message(paste("Saving Excel file to:", OUTPUT_EXCEL_FILE))
saveWorkbook(wb, OUTPUT_EXCEL_FILE, overwrite = TRUE)

# --- VII. Save Combined RDS File ---
message("Saving combined RDS file...")
combined_results <- rbindlist(all_chromosome_results, fill = TRUE)
combined_output <- list(
  data = combined_results,
  summary = summary_stats,
  parameters = list(
    pvalue_filter = PVALUE_COLUMN_FOR_FILTERING,
    pvalue_threshold = PVALUE_THRESHOLD,
    max_tss_distance = MAX_TSS_DISTANCE,
    chromosomes_analyzed = TARGET_CHROMOSOMES,
    total_associations = nrow(combined_results)
  )
)
saveRDS(combined_output, OUTPUT_RDS_FILE)

# --- VIII. Final Summary ---
message("\n=== ANALYSIS COMPLETE ===")
message(sprintf("Analyzed chromosomes: %s", paste(TARGET_CHROMOSOMES, collapse = ", ")))
message(sprintf("Chromosomes with results: %d", length(all_chromosome_results)))
message(sprintf("Total associations across all chromosomes: %s", format(nrow(combined_results), big.mark=",")))
message(sprintf("Total unique pathways: %d", uniqueN(combined_results$Pathway_ID)))
message(sprintf("Total unique genes: %d", uniqueN(combined_results$eQTL_Gene_Entrez)))
message(sprintf("Total unique variants: %d", uniqueN(combined_results$SNP_chr_pos_ref_alt)))

message(sprintf("\nOutput files created:"))
message(sprintf("- Excel file: %s", OUTPUT_EXCEL_FILE))
message(sprintf("- RDS file: %s", OUTPUT_RDS_FILE))

message("\nExcel file contains:")
message("- Summary tab: Overview statistics for all chromosomes")
message(sprintf("- %d chromosome tabs: Individual results for each chromosome", length(all_chromosome_results)))

print(summary_stats)