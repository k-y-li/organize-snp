library(data.table)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(KEGGREST)
library(arrow)
library(RSQLite)
library(DBI)

# --- I. Global Configuration & Paths ---
# Directory containing all the GTEx v10 parquet files
GTEX_DATA_DIR <- "/Users/Kevin/gtex data/GTEx_Analysis_v10_eQTL_updated/"

# Path to the large variant lookup table from GTEx Portal
LOOKUP_TABLE_GZ_FILE <- "/Users/Kevin/gtex data/GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz"

# Path to create/use the SQLite database for fast variant lookups
VARIANT_DB_FILE <- "gtex_variant_lookup.sqlite"

# Main directory to store all output results
OUTPUT_DIR <- "/Users/Kevin/organize-snp/gtex_brain_pathway_results"

# Chromosomes to analyze (autosomal only)
TARGET_CHROMOSOMES <- paste0("chr", 1:22)

# P-value filtering configuration
PVALUE_COLUMN_FOR_FILTERING <- "pval_nominal"
PVALUE_THRESHOLD <- 5e-8 # A more standard threshold for genome-wide significance

# The name of the rsID column in the lookup file
RSID_COLUMN_NAME_IN_LOOKUP <- "rs_id_dbSNP155_GRCh38p13"

# Cache file for KEGG data
KEGG_DATA_RDS_FILE <- "kegg_hsa_pathways_via_bioc_v10_hg38.rds"

# --- II. One-Time Setup: Variant Lookup Database ---
message("\n--- Part 1: Setting up Variant Lookup Database (if needed) ---")

if (!file.exists(VARIANT_DB_FILE)) {
  message(paste("Database not found. Creating", VARIANT_DB_FILE, "from", basename(LOOKUP_TABLE_GZ_FILE)))
  message("This is a one-time setup and may take several minutes...")
  
  if (!file.exists(LOOKUP_TABLE_GZ_FILE)) {
    stop(paste("FATAL: Variant lookup file not found at:", LOOKUP_TABLE_GZ_FILE))
  }
  
  # First, examine the file structure
  message("Examining file structure...")
  header_info <- fread(LOOKUP_TABLE_GZ_FILE, nrows = 0)
  message("Available columns: ", paste(colnames(header_info), collapse = ", "))
  
  # Verify required columns exist
  required_cols <- c("variant_id", "rs_id_dbSNP155_GRCh38p13")
  missing_cols <- setdiff(required_cols, colnames(header_info))
  if (length(missing_cols) > 0) {
    stop(paste("FATAL: Required columns not found:", paste(missing_cols, collapse = ", ")))
  }
  con <- dbConnect(RSQLite::SQLite(), VARIANT_DB_FILE)
  
  # Set up chunked reading parameters
  chunk_size <- 1000000  # 1M rows per chunk
  total_rows_processed <- 0
  chunk_number <- 0
  
  message("Starting chunked data processing...")
  
  # Process file in chunks
  repeat {
    chunk_number <- chunk_number + 1
    message(sprintf("Processing chunk %d (starting from row %s)...", 
                    chunk_number, format(total_rows_processed + 1, big.mark=",")))
    chunk_data <- fread(
      file = LOOKUP_TABLE_GZ_FILE,
      skip = total_rows_processed,
      nrows = chunk_size,
      select = c("variant_id", "rs_id_dbSNP155_GRCh38p13"),
      colClasses = list(character = c("variant_id", "rs_id_dbSNP155_GRCh38p13")),
      showProgress = FALSE
    )
    if (nrow(chunk_data) == 0) {
      message("Reached end of file")
      break
    }
    
    message(sprintf("Read %s rows in chunk %d", format(nrow(chunk_data), big.mark=","), chunk_number))
    setnames(chunk_data, "rs_id_dbSNP155_GRCh38p13", "rsID")
    chunk_data <- chunk_data[!is.na(variant_id) & variant_id != ""]
    
    if (nrow(chunk_data) > 0) {
      dbWriteTable(con, "variants", chunk_data, append = TRUE)
      message(sprintf("Wrote %s rows to database", format(nrow(chunk_data), big.mark=",")))
    }
    
    total_rows_processed <- total_rows_processed + nrow(chunk_data)
    if (nrow(chunk_data) < chunk_size) {
      message("Reached end of file (chunk smaller than expected)")
      break
    }
    
    # Progress update every 5 chunks
    if (chunk_number %% 5 == 0) {
      message(sprintf("Progress: %s total rows processed so far", 
                      format(total_rows_processed, big.mark=",")))
    }
  }

  final_count <- dbGetQuery(con, "SELECT COUNT(*) as count FROM variants")$count
  message(sprintf("Database loading complete. Total rows in database: %s", 
                  format(final_count, big.mark=",")))
  
  if (final_count == 0) {
    dbDisconnect(con)
    file.remove(VARIANT_DB_FILE)
    stop("FATAL: No data was loaded into the database")
  }
  
  message("Creating index on 'variant_id' for fast lookups...")
  start_time <- Sys.time()
  dbExecute(con, "CREATE UNIQUE INDEX IF NOT EXISTS idx_variant_id ON variants (variant_id)")
  end_time <- Sys.time()
  message(sprintf("Index created in %.1f seconds", as.numeric(end_time - start_time, units = "secs")))
  
  # Show some sample data
  sample_data <- dbGetQuery(con, "SELECT * FROM variants LIMIT 5")
  message("Sample data from database:")
  print(sample_data)
  
  dbDisconnect(con)
  message("Database setup complete!")
  
} else {
  message(paste("Found existing variant lookup database:", VARIANT_DB_FILE))
  
  con <- dbConnect(RSQLite::SQLite(), VARIANT_DB_FILE)
  row_count <- dbGetQuery(con, "SELECT COUNT(*) as count FROM variants")$count
  message(sprintf("Existing database contains %s variant records", format(row_count, big.mark=",")))
  dbDisconnect(con)
}


# --- III. One-Time Setup: Load KEGG Pathway Data ---
message("\n--- Part 2: Loading KEGG Pathway Data (if needed) ---")
kegg_pathways_df <- NULL
if (file.exists(KEGG_DATA_RDS_FILE)) {
  message("Loading cached KEGG pathway data from ", KEGG_DATA_RDS_FILE)
  kegg_pathways_df <- as.data.table(readRDS(KEGG_DATA_RDS_FILE))
} else {
  message("Constructing KEGG pathway-gene relationships using org.Hs.eg.db...")
  kegg_mappings <- tryCatch(AnnotationDbi::select(org.Hs.eg.db, 
                                                  keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
                                                  columns = c("PATH", "SYMBOL"), 
                                                  keytype = "ENTREZID"),
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
                                  error = function(e) { message("WARNING: Could not fetch KEGG pathway names via API: ", e$message); NULL })
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
  saveRDS(kegg_pathways_df, KEGG_DATA_RDS_FILE)
  message("Saved KEGG pathway data to ", KEGG_DATA_RDS_FILE)
}
message(paste("Loaded", uniqueN(kegg_pathways_df$pathway_id), "KEGG pathways."))


# --- IV. Function to Process a Single Chromosome ---
process_chromosome <- function(chr_name, gtex_data, kegg_data, db_path, tissue_name_label) {
  message(sprintf("--- Processing %s for %s ---", chr_name, tissue_name_label))
  chr_eqtls <- gtex_data[variant_id %like% paste0(chr_name, "_")]
  if (nrow(chr_eqtls) == 0) {
    message(sprintf("No significant eQTLs found for %s post-filtering.", chr_name))
    return(NULL)
  }
  unique_variants <- unique(chr_eqtls$variant_id)
  con <- dbConnect(RSQLite::SQLite(), db_path)
  rsid_map <- dbGetQuery(con, "SELECT variant_id, rsID FROM variants WHERE variant_id IN (?)",
                         params = list(unique_variants))
  dbDisconnect(con)
  chr_eqtls <- merge(chr_eqtls, as.data.table(rsid_map), by = "variant_id", all.x = TRUE)
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
  chr_output <- merged_data[, .(
    Chromosome = chr_name, Pathway_ID = pathway_id, Pathway_Name = pathway_name,
    Tissue = tissue_name_label, SNP_chr_pos_ref_alt = variant_id, SNP_rsID = rsID,
    eQTL_Gene_Ensembl_Versioned = gene_id, eQTL_Gene_Symbol = gene_symbol_kegg,
    eQTL_Gene_Entrez = entrez_gene_id, eQTL_pval_nominal = pval_nominal,
    eQTL_slope = slope, eQTL_TSS_Distance = tss_distance
  )]
  chr_output <- unique(chr_output)
  message(sprintf("Final %s results: %s unique associations.", chr_name, format(nrow(chr_output), big.mark=",")))
  return(chr_output)
}


# --- V. Main Processing Loop for All Brain Tissues ---
message(sprintf("\n\n=== STARTING ANALYSIS OF ALL BRAIN TISSUES in %s ===", GTEX_DATA_DIR))

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

brain_files_to_process <- list.files(
  path = GTEX_DATA_DIR,
  pattern = "^Brain.*\\.signif_pairs\\.parquet$",
  full.names = TRUE
)

if (length(brain_files_to_process) == 0) {
  stop("No brain tissue parquet files found in the specified directory.")
}

message(sprintf("Found %d brain tissue files to process.", length(brain_files_to_process)))
print(basename(brain_files_to_process))

overall_summary <- data.table()

for (gtex_file_path in brain_files_to_process) {
  
  tissue_name <- sub("\\.v10\\.eQTLs\\.signif_pairs\\.parquet$", "", basename(gtex_file_path))
  message(sprintf("\n\n<<<<< Processing Tissue: %s >>>>>", tissue_name))
  
  tryCatch({
    message(paste("Reading:", basename(gtex_file_path)))
    full_gtex_data <- as.data.table(read_parquet(gtex_file_path))
    message(paste("Loaded", format(nrow(full_gtex_data), big.mark=","), "total significant pairs."))
    
    message(paste("Applying p-value filter:", PVALUE_COLUMN_FOR_FILTERING, "<=", PVALUE_THRESHOLD))
    original_rows <- nrow(full_gtex_data)
    full_gtex_data <- full_gtex_data[get(PVALUE_COLUMN_FOR_FILTERING) <= PVALUE_THRESHOLD & !is.na(get(PVALUE_COLUMN_FOR_FILTERING))]
    message(sprintf("P-value filter: %s -> %s rows", format(original_rows, big.mark=","), format(nrow(full_gtex_data), big.mark=",")))
    
    if (nrow(full_gtex_data) == 0) {
      warning("No eQTLs remaining after filtering. Skipping to next tissue.", immediate. = TRUE)
      next 
    }
    
    all_chromosome_results <- data.table()
    summary_stats <- data.table()
    
    for (chr in TARGET_CHROMOSOMES) {
      chr_result <- process_chromosome(chr, full_gtex_data, kegg_pathways_df, VARIANT_DB_FILE, tissue_name)
      if (!is.null(chr_result) && nrow(chr_result) > 0) {
        all_chromosome_results <- rbind(all_chromosome_results, chr_result, fill = TRUE)
        chr_summary <- data.table(
          Chromosome = chr, Total_Associations = nrow(chr_result),
          Unique_Pathways = uniqueN(chr_result$Pathway_ID),
          Unique_Genes = uniqueN(chr_result$eQTL_Gene_Entrez),
          Unique_Variants = uniqueN(chr_result$SNP_chr_pos_ref_alt)
        )
        summary_stats <- rbind(summary_stats, chr_summary)
      }
      gc()
    }
    
    # --- C. Save Results for the current tissue ---
    if (nrow(all_chromosome_results) == 0) {
      warning("No results were generated for any chromosome for this tissue.", immediate. = TRUE)
    } else {
      
      # --- C1. NEW: Save separate CSV files per chromosome ---
      message(sprintf("\nSaving per-chromosome CSV files for %s...", tissue_name))
      tissue_csv_dir <- file.path(OUTPUT_DIR, tissue_name)
      if (!dir.exists(tissue_csv_dir)) {
        dir.create(tissue_csv_dir, recursive = TRUE)
      }
      
      for (chr_name in unique(all_chromosome_results$Chromosome)) {
        chr_data <- all_chromosome_results[Chromosome == chr_name]
        output_csv_path <- file.path(tissue_csv_dir, sprintf("%s_%s_eQTL_pathways.csv", tissue_name, chr_name))
        fwrite(chr_data, output_csv_path)
      }
      message(sprintf("CSV files saved in: %s/", tissue_csv_dir))
      
      # --- C2. Save combined RDS file for this tissue ---
      output_rds_file <- file.path(OUTPUT_DIR, paste0(tissue_name, "_eQTL_pathways_COMBINED.rds"))
      saveRDS(all_chromosome_results, output_rds_file)
      message(paste("Combined RDS for this tissue saved to:", output_rds_file))
      
      # --- C3. Update overall summary statistics ---
      tissue_summary <- data.table(
        Tissue = tissue_name,
        Total_Associations = nrow(all_chromosome_results),
        Unique_Pathways = uniqueN(all_chromosome_results$Pathway_ID),
        Unique_Genes = uniqueN(all_chromosome_results$eQTL_Gene_Entrez),
        Unique_Variants = uniqueN(all_chromosome_results$SNP_chr_pos_ref_alt)
      )
      overall_summary <- rbind(overall_summary, tissue_summary)
      
      message(sprintf("\nSummary for %s:", tissue_name))
      print(summary_stats)
    }
    
  }, error = function(e) {
    message(sprintf("\n\n !!! ERROR processing %s: %s. Skipping to next file. !!! \n\n", 
                    tissue_name, e$message))
    error_summary <- data.table(Tissue = tissue_name, Total_Associations = "FAILED",
                                Unique_Pathways = NA, Unique_Genes = NA, Unique_Variants = NA)
    overall_summary <- rbind(overall_summary, error_summary, fill = TRUE)
  })
}


# --- VI. Final Overall Summary ---
message("=========== ANALYSIS COMPLETE ===========")

if (nrow(overall_summary) > 0) {
  message("\nOverall Summary Across All Processed Tissues:")
  print(overall_summary)
  
  summary_csv_path <- file.path(OUTPUT_DIR, "00_overall_brain_tissue_summary.csv")
  fwrite(overall_summary, summary_csv_path)
  
  message(sprintf("\nAll results saved in the '%s' directory.", OUTPUT_DIR))
  message("This includes:")
  message(" - A summary file: 00_overall_brain_tissue_summary.csv")
  message(" - A combined .rds file for each tissue (e.g., Brain_Cortex_eQTL_pathways_COMBINED.rds)")
  message(" - A sub-directory for each tissue containing per-chromosome .csv files.")
  
} else {
  message("No tissues were successfully processed.")
}