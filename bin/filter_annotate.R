#!/usr/bin/env Rscript

# Input requirements: =========================================================
# 1. gtf_ref_loc: GTF file with ensembl annotations, currently using hg38 v 102
# 2. gtf_refseq_basename: GTF with refseq annotations, used to annotate XR and NR overlap
# 3. gtf_novel_file: generated using gffcompare to merge all stringtie output GTF files
# 4. tracking_file: generated in same way as gtf_novel_file, contains mappings between XLOC / TCONS IDs and reference ENSG and ENST IDs
# 5. min_occurrence: minimum number of samples in which a transcript must be found to be included in the final GTF
# 6. outfile: output file name
#
# Output: =====================================================================
# 1. outfile: custom annotated GTF

# Track versions:
version = "0.1"
date = "10-06-2024"

# Load required packages
message(paste(Sys.time(), "Loading required libraries ..."), sep = "\t")
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
})

# Global arguments
args <- commandArgs(trailingOnly = TRUE)

gtf_ref_loc = args[1]
gtf_refseq_basename = args[2]
gtf_novel_file = args[3]
tracking_file = args[4]
min_occurrence = args[5]
outfile = args[6]
scripts_dir = args[7]

min_occurrence <- as.numeric(min_occurrence)

# Source additional required functions
functions_file = paste0(scripts_dir, "/filter_annotate_functions.R")
source(functions_file)

fixGTF <- function(gtf_refseq_basename, gtf_ref_path, min_occurrence, novel_gtf_path, tracking_file, output_prefix) {
  output_gtf_path <- paste(output_prefix, "gtf", sep = ".")
  output_log_path <-  paste(output_prefix, "log", sep = ".")

  # Load reference GTF file for comparison and annotation purposes
  gtf_reference <- rtracklayer::import.gff(gtf_ref_path, colnames = c(
    "type", "source", "gene_id", "gene_name", "gene_biotype", "transcript_id", "transcript_name"
  ))
  gtf_ref_df <- as.data.frame(gtf_reference)

  # Import the novel GTF file for filtering and fixing
  novel_gtf <- import(novel_gtf_path)

  # Filtering steps to ensure transcript quality and relevance
  ## Remove transcripts located on scaffolds
  novel_gtf <- novel_gtf[seqnames(novel_gtf) %in% c(1:22, "X", "Y"), ]
  ## Remove unstranded transcripts
  novel_gtf <- novel_gtf[strand(novel_gtf) != "*", ]

  # Convert GTF to a data frame for further processing
  novel_gtf_df <- data.frame(novel_gtf)

  # Discard unwanted transcript classes
  transcripts_discard <- novel_gtf_df[which(novel_gtf_df$type == "transcript" &
                                              novel_gtf_df$class_code %in% c("=", "c", "j", "m", "n", "e", "r", "s")), ]
  novel_gtf_df <- novel_gtf_df[which(!(novel_gtf_df$transcript_id %in% transcripts_discard$transcript_id)), ]
  start_tx_num = length(unique(novel_gtf_df$transcript_id))

  # Load tracking file and merge information with novel GTF data frame
  tracking <- data.table::fread(tracking_file, header = FALSE, select = c(1:4)) %>%
    tidyr::separate(V3, into = c("ref_gene_id_annotated", "ref_transcript_id_annotated"), sep = "\\|") 
  colnames(tracking) <- c("TCONS", "xloc", "ref_gene_id", "ref_transcript_id", "class_code")

  novel_gtf_df <- novel_gtf_df %>%
    left_join(tracking[, c(1:4)], by = c("transcript_id" = "TCONS")) 

  # Carry over information from transcript rows to exon rows
  setDT(novel_gtf_df)
  # Now using data.table's by-reference assignment to fill in NA values with the corresponding values from the 'transcript' rows
  novel_gtf_df[, c('gene_name', 'oId', 'cmp_ref', 'class_code', 'cmp_ref_gene', 'ref_gene_id') := {
    # Locate the 'transcript' row
    w <- which(type == 'transcript')
    # If no 'transcript' row is found, do nothing
    if (length(w) == 0) return(.SD)
    # Carry the values from the 'transcript' row to 'exon' rows
    lapply(.SD, function(x) { x[is.na(x)] <- x[w]; x })
  }, by = transcript_id, .SDcols = c('gene_name', 'oId', 'cmp_ref', 'class_code', 'cmp_ref_gene', 'ref_gene_id')]
  
  # Update the gene_id column using vectorized operations
  novel_gtf_df$gene_id <- ifelse(
    novel_gtf_df$ref_gene_id != "-" & is.na(novel_gtf_df$cmp_ref_gene), 
    novel_gtf_df$ref_gene_id, 
    novel_gtf_df$gene_id)
  
  novel_gtf_df$gene_name <- ifelse(
    is.na(novel_gtf_df$gene_name),
    novel_gtf_df$gene_id, 
    novel_gtf_df$gene_name
  )
  
  # Remove mono-exonic transcripts  
  mono_exonic_novel <- count_mono_exonics(gtf = novel_gtf_df)
  filt_monoexonic <- length(unique(mono_exonic_novel$transcript_id))
  
  novel_gtf_df <- novel_gtf_df[!(novel_gtf_df$transcript_id %in% mono_exonic_novel$transcript_id), ]
  
  
  # Remove transcripts based on the minimum occurrence threshold
  novel_gtf_df$num_samples <- as.numeric(novel_gtf_df$num_samples)
  transcripts_keep <- subset(novel_gtf_df, type == "transcript" & num_samples >= min_occurrence)$transcript_id
  filt_occurrence <- length(unique(novel_gtf_df$transcript_id)) - length(transcripts_keep)
  
  novel_gtf_df <- novel_gtf_df[which(novel_gtf_df$transcript_id %in% transcripts_keep),]
  
  
  # Check for genes with transcripts on multiple strands
  novel_gtf_multiple_strands <- novel_gtf_df %>%
    group_by(gene_id) %>%
    summarize(nstrands = length(unique(strand)))
  
  # Count the number of genes with transcripts on multiple strands
  genes_with_multiple_strands <- sum(novel_gtf_multiple_strands$nstrands > 1)
  
  # Display a message based on the count
  #TODO should move after filtering?
  if (genes_with_multiple_strands > 0) {
    message(paste(genes_with_multiple_strands, "genes have transcripts on multiple strands."))
  } else {
    message("No genes with transcripts on multiple strands found.")
  }
  
  ## Add biotype to custom annotation based on reference ID 
  novel_gtf_df$gene_biotype <- gtf_ref_df$gene_biotype[match(novel_gtf_df$gene_id, gtf_ref_df$gene_id)]
  novel_gtf_df[which(is.na(novel_gtf_df$gene_biotype)), "gene_biotype"] <- "stringtie"
  
  ## Check ref tx overlap and same strand I class transcripts
  ref_overlap_txs <- suppressWarnings(check_tx_overlap(gtf = novel_gtf_df, gtf_reference = gtf_reference))
  same_strand_i_txs <- suppressWarnings(filter_i_class(gtf_df = novel_gtf_df, reference_granges = gtf_reference))
  
  # Filter transcripts overlapping reference transcripts and same strand I class
  novel_gtf_df <- subset(novel_gtf_df,
                               !(novel_gtf_df$transcript_id %in% same_strand_i_txs) & 
                                 !(novel_gtf_df$transcript_id %in% ref_overlap_txs))
  
  final_tx <- length(unique(novel_gtf_df$transcript_id))
  
  ## Flagging RefSeq transcripts
  novel_gtf_df <- annotate_overlap(gtf = novel_gtf_df, gtf_refseq_basename = gtf_refseq_basename, x_name = "xr")
  novel_gtf_df <- annotate_overlap(gtf = novel_gtf_df, gtf_refseq_basename = gtf_refseq_basename, x_name = "nr")
  novel_gtf_df <- data.frame(novel_gtf_df)
  
  xr_overlap <- nrow(subset(novel_gtf_df, type == "transcript" & xr_overlap == "xr_hit")) 
  nr_overlap <- nrow(subset(novel_gtf_df, type == "transcript" & nr_overlap == "nr_hit")) 
  
  novel_gtf_df <- novel_gtf_df[, c("seqnames", "source", "type", "start", "end", "score", "strand", "phase", "gene_id", "transcript_id", "gene_name", "gene_biotype", "oId", "cmp_ref", "class_code", "tss_id", "num_samples", "exon_number", "cmp_ref_gene", "contained_in", "xloc", "ref_gene_id", "ref_transcript_id", "xr_overlap", "nr_overlap")]
  
  
  # WRITE LOG FILE
  
  ## Initialise log file
  cat("#GTF FILTERING", version, "(last updated", date, ")", "\n", file = output_log_path)
  cat("#gtf_ref_path=", gtf_ref_path,"\n", file = output_log_path, append = T)
  cat("#gtf_refseq_basename=", gtf_refseq_basename,"\n", file = output_log_path, append = T)
  cat("#novel_gtf_path=", novel_gtf_path,"\n", file = output_log_path, append = T)
  cat("#tracking_file=", tracking_file,"\n", file = output_log_path, append = T)
  cat("#min_occurrence=", min_occurrence,"\n", file = output_log_path, append = T)
  
  ## Add filtering steps
  cat("Novel transcripts - stranded, not in scaffolds:", "\t",start_tx_num, "\n", file = output_log_path, append = T)
  cat("Filtered monoexonic:", "\t", filt_monoexonic, "\n", file = output_log_path, append = T)
  cat("Filtered min occurence:", "\t", filt_occurrence, "\n", file = output_log_path, append = T)
  cat("Filtered reference overlap:", "\t", length(unique(ref_overlap_txs)), "\n", file = output_log_path, append = T)
  cat("Filtered class I same strand overlap:", "\t", length(unique(same_strand_i_txs)), "\n", file = output_log_path, append = T)
  
  cat("TOTAL NOVEL TRANSCRIPTS ADDED TO REFERENCE:", "\t", final_tx, "\n", file = output_log_path, append = T)
  
  ## Add extra warnings
  cat("INFO - Genes with transcripts on multiple strands:", "\t", genes_with_multiple_strands, "\n", file = output_log_path, append = T)
  cat("INFO - RefSeq xr overlap:", "\t",  xr_overlap, "\n", file = output_log_path, append = T)
  cat("INFO - RefSeq nr overlap:", "\t",  nr_overlap, "\n", file = output_log_path, append = T)
  
  # WRITE FINAL GTF
  ## Turn df into GRanges
  gtf_novel_GR_out <- GenomicRanges::makeGRangesFromDataFrame(novel_gtf_df, keep.extra.columns = T)
  
  gtf_novel_merged <- c(gtf_reference, gtf_novel_GR_out)

  message(paste(Sys.time(), "Exporting custom gtf ... "), sep = "\t")
  export(object = gtf_novel_merged, con = output_gtf_path, format = "gtf", version = "2")
  lines <- readLines(output_gtf_path)
  modified_lines <- gsub("; ID.*", "", lines)
  writeLines(modified_lines, output_gtf_path)
  message(paste(Sys.time(), "Exporting custom gtf ... Done!"), sep = "\t")
}

fixGTF(
  gtf_refseq_basename = gtf_refseq_basename,
  gtf_ref_path = gtf_ref_loc,
  min_occurrence = min_occurrence,
  novel_gtf_path = gtf_novel_file,
  tracking_file = tracking_file,
  output_prefix = outfile
)
