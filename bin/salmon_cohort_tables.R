#!/usr/bin/env Rscript

# Author: Marina Reixachs Sole [M.Reixachssole@prinsesmaximacentrum.nl]

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

# LOGGING FUNCTION -----------------
log_file <- file.path(paste0(prefix, "_log.txt"))
log_message <- function(message, level = "WARNING") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- paste0("[", level, "] ", timestamp, " - ", message)
  write(log_entry, file = log_file, append = TRUE)
  message(log_entry) # Also print to console for visibility
}

# FUNCTIONS -----------------

# Function to generate tx2gene
get_tx2gene <- function(gtf_data) {
  transcript_entries <- as.data.frame(gtf_data) %>%
    filter(type == "transcript")

  if (nrow(transcript_entries) == 0) {
    stop("No transcript entries found in the GTF file. Ensure the GTF file contains transcripts.")
  }

  has_gene_id <- "gene_id" %in% colnames(transcript_entries) && !all(is.na(transcript_entries$gene_id))
  has_gene_name <- "gene_name" %in% colnames(transcript_entries) && !all(is.na(transcript_entries$gene_name))

  # Error handling for missing fields
  if (!has_gene_id && !has_gene_name) {
    log_message("No 'gene_id' or 'gene_name' found in the GTF file. Gene expression will not be computed.")
    tx2gene <- NULL
  } else if (!has_gene_id && has_gene_name) {
    log_message("No 'gene_id' found in the GTF file.")
    tx2gene <- transcript_entries %>% dplyr::select(transcript_id, gene_name)
  } else if (!has_gene_name && has_gene_id) {
    log_message("No 'gene_name' found in the GTF file.")
    tx2gene <- transcript_entries %>% dplyr::select(transcript_id, gene_id)
  } else {
    tx2gene <- transcript_entries %>% dplyr::select(transcript_id, gene_id, gene_name)
  }

  # Count missing values for gene_id and gene_name
  if (has_gene_id) {
    missing_gene_ids <- sum(is.na(tx2gene$gene_id))
    if (missing_gene_ids > 0) {
      log_message(paste(missing_gene_ids, "transcripts have missing 'gene_id' values."))
    }
  }

  if (has_gene_name) {
    missing_gene_names <- sum(is.na(tx2gene$gene_name))
    if (missing_gene_names > 0) {
      log_message(paste(missing_gene_names, "transcripts have missing 'gene_name' values."))
    }
  }

  # Remove rows with missing transcript IDs
  tx2gene <- tx2gene[!is.na(tx2gene$transcript_id), ]

  return(tx2gene)
}

# Function to write properly formatted tables
write_tsv <- function(data, file) {
  if (!is.null(data)) {
    write.table(data,
      file = file,
      sep = "\t",
      row.names = TRUE,
      col.names = TRUE,
      quote = FALSE)
    log_message(paste("Successfully wrote file:", file), level = "INFO")
  } else {
    log_message(paste("Unable to write output file:", file))
  }
}

# CODE -----------------

# Create or clear log file
write("", file = log_file)

# Import gtf_file
gtf_data <- rtracklayer::import(gtf_file)

# Get tx2gene object
tx2gene <- get_tx2gene(gtf_data)

# Save tx2gene table to file
tx2gene_file <- file.path(paste0(prefix, "_tx2gene.tsv"))
write_tsv(tx2gene, tx2gene_file)

# Get the list of quant.sf files
quant_files <- readLines(quant_paths)

if (length(quant_files) == 0) {
  stop("No quant.sf files found in the specified directory.")
}

# Create a named vector for files
names(quant_files) <- basename(dirname(quant_files))

# Import the quantification files using tximport
txi <- tximport(files = quant_files, type = "salmon", txOut = TRUE)

# Summarise gene counts by gene_id
if ("gene_id" %in% colnames(tx2gene)) {
  if (all(tx2gene$transcript_id == tx2gene$gene_id, na.rm = TRUE)) {
    log_message("All 'transcript_id' values are identical to 'gene_id'. Gene ID tables will not be generated.")
    txi.id <- NULL
  } else {
    txi.id <- summarizeToGene(txi, tx2gene[, c("transcript_id", "gene_id")], ignoreTxVersion = TRUE)
  }
} else {
  txi.id <- NULL
}

# Summarise gene counts by gene_name
if ("gene_name" %in% colnames(tx2gene)) {
  if (all(tx2gene$transcript_id == tx2gene$gene_name, na.rm = TRUE)) {
    log_message("All 'transcript_id' values are identical to 'gene_name'. Gene name tables will not be generated.")
    txi.name <- NULL
  } else {
    txi.name <- summarizeToGene(txi, tx2gene[, c("transcript_id", "gene_name")], ignoreTxVersion = TRUE)
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

# Write outputs
write_tsv(txi$counts, counts_file)
write_tsv(txi$abundance, tpms_file)
write_tsv(txi.id$counts, gene_ID_counts_file)
write_tsv(txi.id$abundance, gene_ID_tpms_file)
write_tsv(txi.name$counts, gene_name_counts_file)
write_tsv(txi.name$abundance, gene_name_tpms_file)