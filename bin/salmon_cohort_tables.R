#!/usr/bin/env Rscript
#
# Author: Marina Reixachs Sole [M.Reixachssole@prinsesmaximacentrum.nl]
#
# Input requirements: =========================================================
# 1. quant_paths: text file with paths to all salmon output quant.sf files
# 2. gtf_file: reference or custom gtf for TX2GENE
# 3. prefix: prefix for the output files
#
# Output: =====================================================================
# 1. transcript_counts.txt: file containing transcript counts for all samples in quant_paths
# 2. transcript_tpms.txt: file containing transcript TPM quantifications for all samples in quant_paths
# 3. gene_counts.txt: file containing aggregated gene counts for all samples in quant_paths
# 4. gene_tpms.txt: file containing aggregated gene TPM quantifications for all samples in quant_paths


# LIBRARIES-----------------
suppressPackageStartupMessages({
  library(tximport)
  library(readr)
  library(dplyr)
  library(stringr)
})



# READ ARGS-----------------
args <- commandArgs(trailingOnly = TRUE)

quant_paths = args[1]
gtf_file = args[2]
prefix = args[3]

# FUNCTIONS -----------------

#Function to generate tx2gene
get_tx2gene <- function(gtf_data) {
  # Read the GTF file
  colnames(gtf_data) <- paste("field", c(1:ncol(gtf_data)), sep = "_")
  # Filter transcript entries and extract transcript and gene IDs
  tx2gene <- gtf_data %>%
    filter(field_3 == "transcript") %>%
    mutate(
      transcript_id = str_match(field_9, "transcript_id \"(.*?)\"")[, 2],
      gene_id = str_match(field_9, "gene_id \"(.*?)\"")[, 2],
      gene_name = str_match(field_9, "gene_name \"(.*?)\"")[, 2]
    ) %>%
    select(transcript_id, gene_id, gene_name)
  return(tx2gene)
}

# Function to write properly formatted tables
write_tsv <- function(data, file) {
  write.table(data, file = file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
}


# CODE -----------------

#Read gtf_file
gtf_data <- read_delim(gtf_file, delim = "\t", col_names = FALSE, comment = "#")

# Get tx2gene object
tx2gene <- get_tx2gene(gtf_data)

# Get the list of quant.sf files
#quant_files <- get_quant_files(base_dir)

quant_files <- readLines(quant_paths)

# Check if there are any quant.sf files found
if (length(quant_files) == 0) {
  stop("No quant.sf files found in the specified directory.")
}

# Create a named vector for files
names(quant_files) <- basename(dirname(quant_files))
print(names(quant_files))

# Import the quantification files using tximport
txi <- tximport(files = quant_files, type = "salmon", txOut = TRUE)


# Summarise gene counts by gene_id
txi.id <- summarizeToGene(txi, tx2gene[,c(1:2)], ignoreTxVersion = TRUE)

# Summarise gene counts by gene_name
txi.name <- summarizeToGene(txi, tx2gene[,c(1,3)], ignoreTxVersion = TRUE)


# Define output file names
counts_file <- file.path(paste0(prefix, "_transcript_counts.tsv"))
tpms_file <- file.path(paste0(prefix, "_transcript_tpms.tsv"))

gene_ID_counts_file <- file.path(paste0(prefix, "_gene_ID_counts.tsv"))
gene_ID_tpms_file <- file.path(paste0(prefix, "_gene_ID_tpms.tsv"))

gene_name_counts_file <- file.path(paste0(prefix, "_gene_name_counts.tsv"))
gene_name_tpms_file <- file.path(paste0(prefix, "_gene_name_tpms.tsv"))

# Write counts and TPMs to separate TSV files

write_tsv(txi$counts, counts_file)
write_tsv(txi$tpms, tpms_file)

cat("Successfully wrote transcript-level counts to:", counts_file, "\n")
cat("Successfully wrote transcript-level TPMs to:", tpms_file, "\n")

write_tsv(txi.id$counts, gene_ID_counts_file)
write_tsv(txi.id$abundance, gene_ID_tpms_file)

cat("Successfully wrote gene-level counts (by gene ID) to:", gene_ID_counts_file, "\n")
cat("Successfully wrote gene-level TPMs (by gene ID) to:", gene_ID_tpms_file, "\n")

write_tsv(txi.name$counts, gene_name_counts_file)
write_tsv(txi.name$abundance, gene_name_tpms_file)

cat("Successfully wrote gene-level counts (by gene name) to:", gene_name_counts_file, "\n")
cat("Successfully wrote gene-level TPMs (by gene name) to:", gene_name_tpms_file, "\n")
