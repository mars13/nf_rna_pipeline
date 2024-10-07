#!/usr/bin/env Rscript

# Input requirements: =========================================================
# 1. gtf_ref_loc: GTF file with ensembl annotations, currently using hg38 v 102
# 2. gtf_refseq_basename: GTF with refseq annotations, used to annotate XR and NR overlap
# 3. gtf_novel_file: generated using gffcompare to merge all stringtie output GTF files
# 4. tracking_file: generated in same way as gtf_novel_file, contains mappings between XLOC / TCONS IDs and reference ENSG and ENST IDs
# 5. min_occurrence: minimum number of samples in which a transcript must be found to be included in the final GTF
# 6. min_tpm: minimum expression value in which a transcript must be found to be included in the final GTF
# 7. outfile: output prefix
# 8. scripts_dir: directory containing filter_annotate_functions.R
#
# Output: =====================================================================
# 1. outfile.gtf: filtered custom annotated GTF
# 2. outfile.log: logging file with version info and filtering steps
# 1. outfile.tsv: expression matrix of novel transcripts

# Track versions:
version <- "1.0"
date <- "07-10-2024"

# Load required packages
message(paste(Sys.time(), "Loading required libraries ..."), sep = "\t")
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# Global arguments
args <- commandArgs(trailingOnly = TRUE)

gtf_ref_loc <- args[1]
gtf_refseq_basename <- args[2]
gtf_novel_file <- args[3]
tracking_file <- args[4]
min_occurrence <- args[5]
min_tpm <- args[6]
outfile <- args[7]
scripts_dir <- args[8]

min_occurrence <- as.numeric(min_occurrence)
min_tpm <- as.numeric(min_tpm)

# Source additional required functions
functions_file <- paste0(scripts_dir, "/filter_annotate_functions.R")
source(functions_file)


filterGTF <- function(gtf_refseq_basename, gtf_ref_path, min_occurrence, novel_gtf_path, tracking_file, output_prefix) {
  # Load reference GTF file for comparison and annotation purposes
  gtf_reference <- rtracklayer::import.gff(gtf_ref_path, colnames = c(
    "type", "source", "gene_id", "gene_name", "gene_biotype", "transcript_id", "transcript_name"
  ))

  # Import the novel GTF file for filtering and fixing
  novel_gtf <- rtracklayer::import(novel_gtf_path)

  # Force chromosome style to Ensembl
  seqlevelsStyle(novel_gtf) <- "Ensembl"
  seqlevelsStyle(gtf_reference) <- "Ensembl"

  ## Remove transcripts located on scaffolds
  novel_gtf <- novel_gtf[seqnames(novel_gtf) %in% c(1:22, "X", "Y"), ]
  ## Remove unstranded transcripts
  novel_gtf <- novel_gtf[strand(novel_gtf) != "*", ]

  # Convert GTF to a data frame for further processing
  gtf_ref_df <- as.data.frame(gtf_reference)
  novel_gtf_df <- data.frame(novel_gtf)

  # TRANSCRIPTOME FILTERING ---------------------------

  # FILTER: Discard unwanted transcript classes
  unwanted_transctripts <- c("=", "c", "j", "m", "n", "e", "r", "s")
  transcripts_discard <- novel_gtf_df[which(novel_gtf_df$type == "transcript" &
    novel_gtf_df$class_code %in% unwanted_transctripts), ]
  novel_gtf_df <- novel_gtf_df[!novel_gtf_df$transcript_id %in% transcripts_discard$transcript_id, ]

  # LOGGING: Starting transcript number
  start_tx_num <- length(unique(novel_gtf_df$transcript_id))

  # Load tracking file and merge information with novel GTF data frame
  tracking <- read_tracking_file(tracking_file)

  novel_gtf_df <- novel_gtf_df %>%
    left_join(tracking, by = c("transcript_id" = "TCONS"))

  # Carry over information from transcript rows to exon rows
  setDT(novel_gtf_df)

  # Now using data.table's by-reference assignment to fill in NA values with the corresponding values from the 'transcript' rows
  novel_gtf_df[, c("gene_name", "oId", "cmp_ref", "class_code", "cmp_ref_gene", "ref_gene_id", "num_samples") := { # nolint
    # Locate the 'transcript' row
    w <- which(type == "transcript")
    # If no 'transcript' row is found, do nothing
    if (length(w) == 0) {
      return(.SD)
    }
    # Carry the values from the 'transcript' row to 'exon' rows
    lapply(.SD, function(x) {
      x[is.na(x)] <- x[w]
      x
    })
  }, by = transcript_id, .SDcols = c("gene_name", "oId", "cmp_ref", "class_code", "cmp_ref_gene", "ref_gene_id", "num_samples")]

  # Update the gene_id column using vectorized operations
  novel_gtf_df$gene_id <- ifelse(
    novel_gtf_df$ref_gene_id != "-" & is.na(novel_gtf_df$cmp_ref_gene),
    novel_gtf_df$ref_gene_id,
    novel_gtf_df$gene_id
  )

  novel_gtf_df$gene_name <- ifelse(
    is.na(novel_gtf_df$gene_name),
    novel_gtf_df$gene_id,
    novel_gtf_df$gene_name
  )

  # FILTER: Remove mono-exonic transcripts
  mono_exonic_novel <- count_mono_exonics(gtf = novel_gtf_df)
  novel_gtf_df <- novel_gtf_df[!(novel_gtf_df$transcript_id %in% mono_exonic_novel$transcript_id), ]

  # LOGGING: Number of non mono-exonic transcripts
  filt_monoexonic <- length(unique(mono_exonic_novel$transcript_id))


  # FILTER: Keep transcripts based on minimum occurrence and expression threshold
  occurence_mask <- filter_tx_occurrence(novel_gtf_df, min_occurence = 2, min_tpm = 1)
  transcripts_keep <- novel_gtf_df[occurence_mask, ]$transcript_id

  novel_gtf_df <- novel_gtf_df[occurence_mask, ]

  # LOGGING: Store transcripts failing occurence threshold
  filt_occurrence <- start_tx_num - filt_monoexonic - length(unique(transcripts_keep))


  # Add biotype to custom annotation based on reference ID
  novel_gtf_df$gene_biotype <- gtf_ref_df$gene_biotype[match(novel_gtf_df$gene_id, gtf_ref_df$gene_id)]
  novel_gtf_df[which(is.na(novel_gtf_df$gene_biotype)), "gene_biotype"] <- "stringtie"


  # FILTER: Transcripts overlapping reference transcripts and same strand I class
  ## Check ref tx overlap and same strand I class transcripts
  ref_overlap_txs <- suppressWarnings(check_tx_overlap(gtf = novel_gtf_df, gtf_reference = gtf_reference))
  ## TODO check why is such slow process
  same_strand_i_txs <- suppressWarnings(filter_i_class(gtf_df = novel_gtf_df, reference_granges = gtf_reference))

  novel_gtf_df <- subset(
    novel_gtf_df,
    !(novel_gtf_df$transcript_id %in% same_strand_i_txs) &
      !(novel_gtf_df$transcript_id %in% ref_overlap_txs)
  )
  # LOGGING: Get final number of transcripts to be added
  final_tx <- length(unique(novel_gtf_df$transcript_id))

  # WARNINGS: Check for genes with transcripts on multiple strands
  novel_gtf_multiple_strands <- novel_gtf_df %>%
    group_by(gene_id) %>%
    summarize(nstrands = length(unique(strand)))

  # Count the number of genes with transcripts on multiple strands
  genes_with_multiple_strands <- sum(novel_gtf_multiple_strands$nstrands > 1)

  # Display a message based on the count
  if (genes_with_multiple_strands > 0) {
    message(paste(genes_with_multiple_strands, "genes have transcripts on multiple strands."))
  } else {
    message("No genes with transcripts on multiple strands found.")
  }

  # WARNINGS: Flagging RefSeq transcripts
  novel_gtf_df <- annotate_overlap(gtf = novel_gtf_df, gtf_refseq_basename = gtf_refseq_basename, x_name = "xr")
  novel_gtf_df <- annotate_overlap(gtf = novel_gtf_df, gtf_refseq_basename = gtf_refseq_basename, x_name = "nr")
  novel_gtf_df <- data.frame(novel_gtf_df)

  xr_overlap <- nrow(subset(novel_gtf_df, type == "transcript" & xr_overlap == "xr_hit"))
  nr_overlap <- nrow(subset(novel_gtf_df, type == "transcript" & nr_overlap == "nr_hit"))

  # Output info
  ## TODO can this tsv have the same headers as the sample IDs?
  expression_columns <- colnames(novel_gtf_df)[grepl("TPM_q", colnames(novel_gtf_df))]
  output_info <- unique(novel_gtf_df[, c("transcript_id", "gene_id", expression_columns)])

  # Output gtf
  output_gtf <- novel_gtf_df[, c(
    "seqnames",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "gene_id",
    "transcript_id",
    "gene_name",
    "gene_biotype",
    "oId",
    "cmp_ref",
    "class_code",
    "tss_id",
    "num_samples",
    "exon_number",
    "cmp_ref_gene",
    "contained_in",
    "xloc",
    "ref_gene_id",
    "ref_transcript_id",
    "xr_overlap",
    "nr_overlap"
  )]

  # DEFINE OUTPUT FILES ---------------------------
  output_gtf_path <- paste(output_prefix, "gtf", sep = ".")
  output_info_path <- paste(output_prefix, "tsv", sep = ".")
  output_log_path <- paste(output_prefix, "log", sep = ".")

  # WRITE LOG FILE ---------------------------

  ## Initialise log file
  cat("#GTF FILTERING", version, "(last updated", date, ")", "\n", file = output_log_path)
  cat("#gtf_ref_path=", gtf_ref_path, "\n", file = output_log_path, append = T)
  cat("#gtf_refseq_basename=", gtf_refseq_basename, "\n", file = output_log_path, append = T)
  cat("#novel_gtf_path=", novel_gtf_path, "\n", file = output_log_path, append = T)
  cat("#tracking_file=", tracking_file, "\n", file = output_log_path, append = T)
  cat("#min_occurrence=", min_occurrence, "\n", file = output_log_path, append = T)
  cat("#min_tpm=", min_tpm, "\n", file = output_log_path, append = T)

  ## Add filtering steps
  cat("Novel transcripts - stranded, not in scaffolds:", "\t", start_tx_num, "\n", file = output_log_path, append = T)
  cat("Filtered monoexonic:", "\t", filt_monoexonic, "\n", file = output_log_path, append = T)
  cat("Filtered min", min_tpm, "TPM in", min_occurrence, "samples:", "\t", filt_occurrence, "\n", file = output_log_path, append = T)
  cat("Filtered reference overlap:", "\t", length(unique(ref_overlap_txs)), "\n", file = output_log_path, append = T)
  cat("Filtered class I same strand overlap:", "\t", length(unique(same_strand_i_txs)), "\n", file = output_log_path, append = T)
  cat("__________________________________________", file = output_log_path, append = T)
  cat("TOTAL NOVEL TRANSCRIPTS ADDED TO REFERENCE:", "\t", final_tx, "\n", file = output_log_path, append = T)
  cat("__________________________________________", file = output_log_path, append = T)
  ## Add extra warnings
  cat("Info - Genes with transcripts on multiple strands:", "\t", genes_with_multiple_strands, "\n", file = output_log_path, append = T)
  cat("Info - RefSeq xr overlap:", "\t", xr_overlap, "\n", file = output_log_path, append = T)
  cat("Info - RefSeq nr overlap:", "\t", nr_overlap, "\n", file = output_log_path, append = T)

  # WRITE FINAL GTF AND INFO TABLE ---------------------------

  ## Turn df into GRanges
  gtf_novel_GR_out <- GenomicRanges::makeGRangesFromDataFrame(output_gtf, keep.extra.columns = T)

  gtf_novel_merged <- c(gtf_reference, gtf_novel_GR_out)

  message(paste(Sys.time(), "Exporting custom gtf ... "), sep = "\t")
  export(object = gtf_novel_merged, con = output_gtf_path, format = "gtf", version = "2")
  lines <- readLines(output_gtf_path)
  modified_lines <- gsub("; ID.*", "", lines)
  writeLines(modified_lines, output_gtf_path)
  write.table(output_info, file = output_info_path, sep = "\t", row.names = F, quote = F)
  message(paste(Sys.time(), "Exporting custom gtf ... Done!"), sep = "\t")
}

# EXECUTE FILTERING  ========================================================
filterGTF(
  gtf_refseq_basename = gtf_refseq_basename,
  gtf_ref_path = gtf_ref_loc,
  min_occurrence = min_occurrence,
  novel_gtf_path = gtf_novel_file,
  tracking_file = tracking_file,
  output_prefix = outfile
)
