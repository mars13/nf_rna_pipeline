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
  library(rtracklayer)
  library(tximport)
  library(readr)
  library(dplyr)
  library(stringr)
})

# READ ARGS-----------------
args <- commandArgs(trailingOnly = TRUE)

quant_paths <- args[1]
gtf_file <- args[2]
prefix <- args[3]

# FUNCTIONS -----------------

# Function to generate tx2gene
get_tx2gene <- function(gtf_data) {
  # Filter transcript entries
  transcript_entries <- as.data.frame(gtf_data) %>%
    filter(type == "transcript")

  if (nrow(transcript_entries) == 0) {
    stop("No transcript entries found in the GTF file. Ensure the GTF file contains transcripts.")
  }

  has_gene_id <- "gene_id" %in% colnames(transcript_entries) || !all(is.na(transcript_entries$gene_id))
  has_gene_name <- "gene_name" %in% colnames(transcript_entries) || !all(is.na(transcript_entries$gene_name))


  # Error handling for missing fields
  if (!has_gene_id && !has_gene_name) {
    warning("No 'gene_id' or 'gene_name' found in the GTF file. Gene expression will not be computed.")

    tx2gene <- NULL
  } else if (!has_gene_id && has_gene_name) {
    warning("No 'gene_id' found in the GTF file.")

    tx2gene <- transcript_entries %>%
      dplyr::select(transcript_id, gene_name)
  } else if (!has_gene_name && has_gene_id) {
    warning("No 'gene_name' found in the GTF file.")

    tx2gene <- transcript_entries %>%
      dplyr::select(transcript_id, gene_id)
  } else {
    tx2gene <- transcript_entries %>%
      dplyr::select(transcript_id, gene_id, gene_name)
  }

  # Remove rows with missing transcript IDs
  tx2gene <- tx2gene[!is.na(tx2gene$transcript_id), ]

  return(tx2gene)
}

# Function to write properly formatted tables
write_tsv <- function(counts_data, file) {
  if (!is.null(counts_data)) {
    write.table(counts_data,
                file = file,
                sep = "\t",
                row.names = TRUE,
                col.names = TRUE,
                quote = FALSE)
    message(paste0("Successfully wrote file: ", file))
  } else {
    warning(paste0("Unable to write output file: ", file))
  }
}


# CODE -----------------

# Import gtf_file
gtf_data <- rtracklayer::import(gtf_file)

# Get tx2gene object
tx2gene <- get_tx2gene(gtf_data)

# Get the list of quant.sf files
quant_files <- readLines(quant_paths)

# Check if there are any quant.sf files found
if (length(quant_files) == 0) {
  stop("No quant.sf files found in the specified directory.")
}

# Create a named vector for files
names(quant_files) <- basename(dirname(quant_files))

# Import the quantification files using tximport
txi <- tximport(files = quant_files, type = "salmon", txOut = TRUE)

# Summarise gene counts by gene_id

# Check if gene_id is present
if ("gene_id" %in% colnames(tx2gene)) {
  # Check for identical transcript_id and gene_id
  if (all(tx2gene$transcript_id == tx2gene$gene_id, na.rm = TRUE)) {
    warning("All 'transcript_id' values are identical to 'gene_id'. Gene ID tables will not be generated.")
    txi.id <- NULL
  } else {
    txi.id <- summarizeToGene(txi, tx2gene[, c("transcript_id", "gene_id")], ignoreTxVersion = TRUE)
  }
} else {
  txi.id <- NULL
}

# Summarise gene counts by gene_name

# Check if gene_name is present
if ("gene_name" %in% colnames(tx2gene)) {
  # Check for identical transcript_id and gene_name
  if (all(tx2gene$transcript_id == tx2gene$gene_name, na.rm = TRUE)) {
    warning("All 'transcript_id' values are identical to 'gene_name'. Gene name tables will not be generated.")
    txi.name <- NULL
  } else {
    txi.name <- summarizeToGene(txi,tx2gene[, c("transcript_id", "gene_id")], ignoreTxVersion = TRUE)
  }
} else {
  txi.name <- NULL
}

# Define output file names
counts_file <- file.path(paste0(prefix, "_transcript_counts.tsv"))
tpms_file <- file.path(paste0(prefix, "_transcript_tpms.tsv"))

gene_ID_counts_file <- file.path(paste0(prefix, "_gene_ID_counts.tsv"))
gene_ID_tpms_file <- file.path(paste0(prefix, "_gene_ID_tpms.tsv"))

gene_name_counts_file <- file.path(paste0(prefix, "_gene_name_counts.tsv"))
gene_name_tpms_file <- file.path(paste0(prefix, "_gene_name_tpms.tsv"))

# Write counts and TPMs to separate TSV files
write_tsv(txi$counts, counts_file)
write_tsv(txi$abundance, tpms_file)

write_tsv(txi.id$counts, gene_ID_counts_file)
write_tsv(txi.id$abundance, gene_ID_tpms_file)

write_tsv(txi.name$counts, gene_name_counts_file)
write_tsv(txi.name$abundance, gene_name_tpms_file)