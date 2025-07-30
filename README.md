# Overview

**Single nucleotide polymorphisms (SNPs)** are sequence differences that affect a single nucleotide and are considered to be the most common genetic variation type in the human genome. SNPs are highly abundant throughout the human genome, occurring at a rate of approximately one every 1,000 base pairs. While most SNPs are located in the non-coding regions of a DNA sequence (AKA the exon) and therefore have no noticeable effect, some SNPs can have functions that lead to protein structure differences or gene expression regulation. These SNPs are referred to as **eQTLs**, or **expression quantitative trait loci**, and are formally defined as the locations in a species genome where a genetic variation is associated with differing gene expression levels (Nica & Dermitzakis, 2013). Specifically, a lot of research focuses on the purpose of *cis-eQTLs*, which are located close to and regulate one or several nearby genes (approximately 1 million bp from a transcription start site) (Pritchard, 2023).

eQTLs are essential for genome-wide association studies (GWAS) because a SNP being a cis-eQTL for a gene can provide a hypothesis that such a mutation can influence trait risk by changing that gene's expression level (Uffelmann et al., 2021). However, while we can use GWAS to associate SNPs with traits and diseases, it might be difficult to determine the mechanisms behind these associations without further analysis. One way researchers do so is through **partitioning heritability**, in which *heritability* is the proportion of disease and trait variation explained by genetic factors. Instead of treating all the variation as one group, we can "partition" SNPs and eQTLs into categories that represent common biological functions and pathways (Finucane et al., 2015). Once we fully associate our eQTLs through partitioning, we can then perform statistical analysis to determine if a pathway contributes to a trait's heritability more than expected by chance.

In short, this project aims to **perform a brain tissue eQTL-pathway analysis**. We will do this by **partitioning SNPs into groups based on the biological pathways they influence** and **creating a database that potentially associates eQTLs with heavily-researched traits and diseases**. Our specific dataset focuses on brain eQTLs downloaded from version 10 of the public *Adult Genotype-Tissue Expression (GTEx)* project and associates them with *KEGG (Kyoto Encyclopedia of Genes and Genomes)* biological pathways to facilitate pathway-based heritability partitioning analyses. All genomic coordinates and variant annotations in this analysis are based on the *GRCh38/hg38 human genome reference assembly (patch 13)*, ensuring compatibility with current genomic databases and downstream GWAS applications. These resources are pooled together to map genetic variants that could significantly affect gene expression in brain tissue and connect these genes to known biological pathways, creating new datasets that could potentially help researchers understand the genetic architecture behind human neurological traits and psychiatric disorders.

# Prerequisites

This dataset was created using R version 4.5.0, running in R Studio version 2024.12.1+563.

## R Packages Used

### Data Processing

```{r load-data-packages, eval=FALSE}
library(data.table)  # Fast data manipulation and filtering of large eQTL datasets
library(dplyr)       # Data wrangling and transformation
library(arrow)       # Reading GTEx parquet files
```

### Database Management

```{r load-db-packages, eval=FALSE}
library(RSQLite)     # Creating and querying variant lookup database for rsID mapping
library(DBI)         # Database interface
```

### Genomic Annotations

```{r load-genomic-packages, eval=FALSE}
library(AnnotationDbi)  # Gene annotation database interface
library(org.Hs.eg.db)   # Human gene annotations (Ensembl to Entrez ID mapping)
library(KEGGREST)       # Downloading KEGG pathway data and annotation
```

## Input Data Files

### GTEx v10 Files

-   Brain tissue eQTL files in `.signif_pairs.parquet` format
    -   Available in [Open Access Datasets \> QTL \> GTEx Analysis V10](https://www.gtexportal.org/home/downloads/adult-gtex/reference)
-   Variant lookup table: `GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz`
    -   Available in [Open Access Datasets \> References \> GTEx Analysis V10](https://www.gtexportal.org/home/downloads/adult-gtex/reference)

### Generated During Analysis

-   SQLite database for variant-to-rsID mapping (created automatically)
-   Cached KEGG pathway data (downloaded automatically on first run)

# Configuration

Our analysis requires two main input files from the GTEx v10 dataset (GTEx Consortium, 2020). The main source of data we will obtain our SNPs from will consist of brain tissue eQTL files in parquet format, which follow the naming convention `Brain_[TissueName].v10.eQTLs.signif_pairs.parquet` and contain pre-filtered significant eQTL associations. These parquet files include essential columns such as:

-   `variant_id` (genomic coordinates in `chr_pos_ref_alt` format)
-   `gene_id` (Ensembl gene IDs with version numbers)
-   `pval_nominal` (nominal p-values for statistical filtering)
-   `slope` (effect sizes)
-   `tss_distance` (distance from transcription start site)

The second required input file is the GTEx variant lookup table named `GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz`. Since the parquet files do not have reference SNP cluster IDs or *rsIDs* (which is the main naming convention used to identify SNPs with), this file serves to map between GTEx variant identified and dbSNP version 151 rsIDs. This compressed text file contains columns including:

-   `variant_id` (matching the format in parquet files)
-   `rs_id_dbSNP155_GRCh38p13` (corresponding dbSNP rsIDs)

These columns enable the annotation of genetic variants with standardized database identifiers essential for downstream analysis and integration with other genomic datasets.

# Methodology

## Variant Lookup Database Creation

Our analysis begins with creating an SQLite database for variant-to-rsID mapping, which addresses the challenge of needing to repeatedly refer to a multi-gigabyte lookup table during dataset creation. The system first checks for an existing database and only performs the setup when necessary:

```{r variant-db-setup, eval=FALSE}
if (!file.exists(VARIANT_DB_FILE)) {
  message(paste("Database not found. Creating", VARIANT_DB_FILE, 
                "from", basename(LOOKUP_TABLE_GZ_FILE)))
con <- dbConnect(RSQLite::SQLite(), VARIANT_DB_FILE)
```

The code then attempts to read the entire compressed lookup table into memory, but can fall back to chunked processing should memory constraints be encountered:

```{r variant-db-read, eval=FALSE}
tryCatch({
  lookup_data <- fread(
    LOOKUP_TABLE_GZ_FILE,
    select = c("variant_id", RSID_COLUMN_NAME_IN_LOOKUP),
    colClasses = c("character", "character")
  )
  setnames(lookup_data, RSID_COLUMN_NAME_IN_LOOKUP, "rsID")
  dbWriteTable(con, "variants", lookup_data, append = TRUE)
}, error = function(e) {
  # Chunked reading fallback for memory-constrained systems
  chunk_size <- 500000
  # ... chunked processing code ...
})
```

While creating the SQLite database, a unique index is created on the `variant_id` column of the lookup table to assist with rapid lookups during chromosome processing. This entire process results in a database that gives us file scans that would take mere minutes instead of the hours it might take with common packages like `data.table`.

## KEGG Pathway Data Caching

Since our analysis requires using the KEGG API, it might be useful to cache the data so we don't have to constantly depend on external connections should we need to run our code again. Like the previous section, our pipeline first checks for existing cached data and loads it directly if available:

```{r kegg-cache-check, eval=FALSE}
kegg_pathways_df <- NULL
if (file.exists(KEGG_DATA_RDS_FILE)) {
  message("Loading cached KEGG pathway data from ", KEGG_DATA_RDS_FILE)
  kegg_pathways_df <- as.data.table(read_rds(KEGG_DATA_RDS_FILE))
}
```

If the code cannot find cached data, it proceeds to use the `org.Hs.eg.db` Bioconductor annotation package to construct KEGG pathway-gene relationships while attempting to perform error handling to manage potential database connectivity issues:

```{r kegg-construct, eval=FALSE}
message("Constructing KEGG pathway-gene relationships using org.Hs.eg.db...")
kegg_mappings <- tryCatch(AnnotationDbi::select(org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
  columns = c("PATH", "SYMBOL"),
  keytype = "ENTREZID"),
  error = function(e) { message("Error: ", e$message); NULL })
```

The pathway data construction itself involves filtering for genes with valid KEGG annotations and standardizing pathway identifiers to ensure compatibility with KEGG databases:

```{r kegg-filter, eval=FALSE}
kegg_mappings_filtered <- kegg_mappings[!is.na(kegg_mappings$PATH), ]
if (nrow(kegg_mappings_filtered) == 0) {
  stop("No genes with KEGG pathway annotations in org.Hs.eg.db.")
}

kegg_pathways_dt <- data.table(
  pathway_id_num = kegg_mappings_filtered$PATH,
  entrez_gene_id_kegg = as.character(kegg_mappings_filtered$ENTREZID),
  gene_symbol_from_bioc = kegg_mappings_filtered$SYMBOL
)

kegg_pathways_dt[, pathway_id := paste0("hsa", pathway_id_num)][, pathway_id_num := NULL]
```

Human-readable pathway names are retrieved from the KEGG API with fallback handling for network connectivity issues, ensuring the analysis can proceed even when external services are unavailable. Eventually, our completed KEGG pathway dataset is then cached as an RDS file for future analyses, eliminating the need for repeated API calls and ensuring consistent annotations across analysis runs.

## Core Processing Pipeline

The heart of the analysis centers around the `process_chromosome()` function, which systematically transforms raw eQTL data into pathway-annotated associations for individual chromosomes. This function handles the complex data integration challenges inherent in multi-database genomic analysis:

```{r process-chromosome-function, eval=FALSE}
process_chromosome <- function(chr_name, gtex_data, kegg_data, db_path, tissue_name_label) {
  message(sprintf("--- Processing %s for %s ---", chr_name, tissue_name_label))
  
  chr_eqtls <- gtex_data[variant_id %like% paste0(chr_name, "_")]
  if (nrow(chr_eqtls) == 0) {
    message(sprintf("No significant eQTLs found for %s post-filtering.", chr_name))
    return(NULL)
  }
```

The function begins by filtering the tissue-specific eQTL data to the target chromosome using pattern matching on variant coordinates, immediately returning `NULL` for chromosomes without significant associations to prevent downstream processing errors. Variant annotation leverages the pre-built SQLite database for rapid rsID mapping:

```{r variant-annotation, eval=FALSE}
  unique_variants <- unique(chr_eqtls$variant_id)
  con <- dbConnect(RSQLite::SQLite(), db_path)
  rsid_map <- dbGetQuery(con, "SELECT variant_id, rsID FROM variants WHERE variant_id IN (?)",
    params = list(unique_variants))
  dbDisconnect(con)
  chr_eqtls <- merge(chr_eqtls, as.data.table(rsid_map), by = "variant_id", all.x = TRUE)
```

Gene identifier conversion represents a critical step where Ensembl gene IDs from GTEx are mapped to Entrez gene IDs required for KEGG pathway assignment:

```{r gene-conversion, eval=FALSE}
  unique_ensembl_ids <- unique(sub("\\..*$", "", chr_eqtls$gene_id))
  entrez_map <- mapIds(org.Hs.eg.db, keys = unique_ensembl_ids, column = "ENTREZID",
    keytype = "ENSEMBL", multiVals = "first")
  ensembl_to_entrez_dt <- data.table(
    ensembl_id_no_version = names(entrez_map[!is.na(entrez_map)]),
    entrez_gene_id = entrez_map[!is.na(entrez_map)]
  )
```

The main processing loop performs tissue-level analysis with error handling to ensure robustness across varying data quality and completeness:

```{r main-processing-loop, eval=FALSE}
for (gtex_file_path in brain_files_to_process) {
  tissue_name <- sub("\\.v10\\.eQTLs\\.signif_pairs\\.parquet$", "", basename(gtex_file_path))
  message(sprintf("\n\n<<<<< Processing Tissue: %s >>>>>", tissue_name))
  
  tryCatch({
    full_gtex_data <- as.data.table(read_parquet(gtex_file_path))
    message(paste("Loaded", format(nrow(full_gtex_data), big.mark=","), "total significant pairs."))
    
    original_rows <- nrow(full_gtex_data)
    full_gtex_data <- full_gtex_data[get(PVALUE_COLUMN_FOR_FILTERING) <= PVALUE_THRESHOLD &
      !is.na(get(PVALUE_COLUMN_FOR_FILTERING))]
```

Each tissue undergoes systematic chromosome-by-chromosome processing to ensure the entire human genome is being covered:

```{r chromosome-processing, eval=FALSE}
    for (chr in TARGET_CHROMOSOMES) {
      chr_result <- process_chromosome(chr, full_gtex_data, kegg_pathways_df, VARIANT_DB_FILE, tissue_name)
      if (!is.null(chr_result) && nrow(chr_result) > 0) {
        all_chromosome_results[[chr]] <- chr_result
        chr_summary <- data.table(
          Chromosome = chr, 
          Total_Associations = nrow(chr_result),
          Unique_Pathways = uniqueN(chr_result$Pathway_ID),
          Unique_Genes = uniqueN(chr_result$eQTL_Gene_Entrez),
          Unique_Variants = uniqueN(chr_result$SNP_chr_pos_ref_alt)
        )
        summary_stats <- rbind(summary_stats, chr_summary)
      }
      gc()
    }
```

Output generation follows a structured approach with each tissue has one CSV file for each chromosome, as well as an RDS file where all chromosomes are combined for more comprehensive analysis:

```{r output-generation, eval=FALSE}
    tissue_csv_dir <- file.path(OUTPUT_DIR, tissue_name)
    for (chr_name in names(all_chromosome_results)) {
      chr_data <- all_chromosome_results[[chr_name]]
      output_csv_path <- file.path(tissue_csv_dir, 
        sprintf("%s_%s_eQTL_pathways.csv", tissue_name, chr_name))
      fwrite(chr_data, output_csv_path)
    }
    
    combined_results <- rbindlist(all_chromosome_results, fill = TRUE)
    output_rds_file <- file.path(OUTPUT_DIR, paste0(tissue_name, "_eQTL_pathways_COMBINED.rds"))
    saveRDS(combined_results, output_rds_file)
```

The pipeline concludes by generating summary statistics across all tissues, providing insights into analysis completeness and data coverage:

```{r summary-generation, eval=FALSE}
if (length(overall_summary) > 0) {
  final_summary_table <- rbindlist(overall_summary, fill = TRUE)
  message("\nOverall Summary Across All Processed Tissues:")
  print(final_summary_table)
  summary_csv_path <- file.path(OUTPUT_DIR, "00_overall_brain_tissue_summary.csv")
  fwrite(final_summary_table, summary_csv_path)
}
```

# Output

The main output we attain from our code are the per-chromosome CSV files (`[Tissue]_[Chromosome]_eQTL_pathways.csv`) that contain the eQTL-pathway associations for a brain tissue and chromosome combination, with each row representing a variant-gene-pathway association. Each row consists of the same column structure:

-   **genomic coordinates** (`SNP_chr_pos_ref_alt`)
-   **dbSNP identifiers** (`SNP_rsID`)
-   **gene annotations** (`eQTL_Gene_Ensembl_Versioned`, `eQTL_Gene_Symbol`, `eQTL_Gene_Entrez`)
-   **pathway information** (`Pathway_ID`, `Pathway_Name`)
-   **eQTL statistics** (`eQTL_pval_nominal`, `eQTL_slope`, `eQTL_TSS_Distance`)

We also include compressed RDS data files that contain the entire eQTL-pathway association dataset for each brain tissue across all processed chromosomes (`[Tissue]_eQTL_pathways_COMBINED.rds`) while maintaining the same column structure as the per-chromosome CSV files. The RDS format allows for researchers to quickly load large datasets into R environments, allowing for tissue-level analysis and statistical modeling.

The master summary file (`00_overall_brain_tissue_summary.csv`) provides analysis completeness information across all processed brain tissues, giving counts for:

-   total associations identified
-   unique pathways represented
-   unique genes analyzed
-   unique variants processed

## Data Interpretation Guidelines

Each row consists of a statistically significant eQTL (p-value \< 5e-8) where a genetic variant influences gene expression, and that gene participates in a KEGG pathway. The `eQTL_slope` values indicate effect direction and magnitude, meaning that:

-   a **positive value** shows that an eQTL increases gene expression relative to the wild-type sequence
-   a **negative value** suggests decreased expression

The magnitude itself is measured in standard deviations of normalized expression levels, allowing researchers to focus on eQTLs with larger slope values if they choose to do so (Mohammadi et al., 2017).

`eQTL_TSS_Distance` provides information about the regulatory relationship between variant and gene:

-   **Negative distances** indicate variants located upstream of the transcription start site (typically in promoter regions)
-   **Positive distances** represent downstream locations (often in gene bodies or 3' regions)

Variants closer to the TSS are more likely to affect core promoter function, while more distant variants may represent enhancer or silencer elements. The values we obtain from the `eQTL_TSS_Distance` column can be useful for interpreting the plausibility of eQTL associations when compared with other characteristics, like the aforementioned `eQTL_slope` column or transcription factor binding site predictions (GTEx Consortium, 2017).

# References

-   Finucane, H. K., Bulik-Sullivan, B., Gusev, A., Trynka, G., et al. (2015). *Partitioning heritability by functional annotation using genome-wide association summary statistics.* Nature Genetics, 47(11), 1228-1235.
-   GTEx Consortium (2017). *Genetic effects on gene expression across human tissues*. Nature, 550(7675), 204-213.
-   GTEx Consortium. (2020). *The GTEx Consortium atlas of genetic regulatory effects across human tissues.* Science, 369(6509), 1318-1330.
-   Mohammadi, P., Castel, S. E., Brown, A. A., & Lappalainen, T. (2017). *Quantifying the regulatory effect size of cis-acting genetic variation using allelic fold change.* Genome Research, 27(11), 1872-1884.
-   Nica, A. C., & Dermitzakis, E. T. (2013). *Expression quantitative trait loci: present and future.* Philosophical Transactions of the Royal Society B: Biological Sciences, 368(1620), 20120362.
-   Pritchard, J. (2023). *An owner's guide to the human genome: An introduction to human population genetics, variation and disease.* Pritchard Lab, Stanford University.
-   Uffelmann, E., Huang, Q. Q., Munung, N. S., de Vries, J., et al. (2021). *Genome-wide association studies.* Nature Reviews Methods Primers, 1(1).

## R Package References

-   Barrett, T., Dowle, M., Srinivasan, A., Gorecki, J., Chirico, M., Hocking, T., Schwendinger, B., & Krylov, I. (2025). *data.table: Extension of data.frame [R package]*. Version 1.17.4.
-   Carlson, M. (2025). *org.Hs.eg.db: Genome wide annotation for Human [R package]*. Version 3.21.0.
-   Müller, K., Wickham, H., James, D. A., & Falcon, S. (2025). *RSQLite: SQLite Interface for R [R package]*. Version 2.4.1.
-   Pagès, H., Carlson, M., Falcon, S., & Li, N. (2025). *AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor [R package]*. Version 1.70.0.
-   R Special Interest Group on Databases, Wickham, H., & Müller, K. (2024). *DBI: R Database Interface [R package]*. Version 1.2.3.
-   Richardson, N., Cook, I., Crane, N., Dunnington, D., François, R., Keane, J., Moldovan-Grünfeld, D., Ooms, J., Wujciak-Jens, J., & Apache Arrow. (2025). *arrow: Integration to 'Apache' 'Arrow' [R package]*. Version 20.0.0.2.
-   Tenenbaum, D., & Maintainer, B. (2025). *KEGGREST: Client-side REST access to the Kyoto Encyclopedia of Genes and Genomes (KEGG) [R package]*. Version 1.48.0.
