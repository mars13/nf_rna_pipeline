# GTF filtering functions

count_mono_exonics <- function(gtf) {
  # Input:
  # gtf_novel = Granges object containing novel transcripts
  #
  # Finds transcripts with a single exon and removes those from the object
  # Remove transcripts that are found by StringTie with a single exon

  # Count exons per transcript
  # Subtract 1 to account for transcript entry
  transcript_counts <- table(gtf$transcript_id) - 1
  # Monoexonic has fewer than 2 entries
  monoexonic_tx_ids <- names(transcript_counts[transcript_counts < 2])

  # Filter for monoexonic transcripts from StringTie and novel genes
  monoexonic_transcripts <- subset(
    gtf,
    transcript_id %in% monoexonic_tx_ids &
      source == "StringTie" &
      ref_gene_id == "-"
  )
  return(monoexonic_transcripts)
}

###########################################################################

filter_tx_occurrence <-
  function(gtf_df_tracking,
           min_occurrence = 1,
           min_tpm = 1) {
    # Input:
    # gtf_df_tracking = data.frame with output from GFFcompare
    #                   with tracking information for TPM values
    # min_occurrence = integer (default = 1) how many samples should share a transfrag
    # min_tpm = integer (default = 1) how much coverage each sample should have to count as covered
    #
    # Subsets transcripts on minimal occurrence and expression parameter.

    tpm_cols <- colnames(gtf_df_tracking)[grepl("TPM_q", colnames(gtf_df_tracking))]

    # Convert TPM columns to a matrix for faster computation
    tpm_matrix <- as.matrix(gtf_df_tracking[, ..tpm_cols])

    # Replace NA's and convert to numeric
    tpm_matrix[is.na(tpm_matrix)] <- 0
    tpm_matrix <- apply(tpm_matrix, 2, as.numeric)

    # Count the number of samples with TPM >= 1 for each row
    tpm_count <- rowSums(tpm_matrix >= min_tpm)

    # Filter the rows where at least 2 samples have TPM >= 1
    to_keep <- tpm_count >= min_occurrence

    return(to_keep)
  }

###########################################################################

check_tx_overlap <- function(gtf, gtf_reference) {
  # gtf_novel_processed = input from previous function
  # gtf_reference = ensembl reference GTF Grange
  #
  # This function looks for exonic overlap between query (novel tx exons)
  # and subject (ref tx exons). Novel transcripts that overlap with multiple
  # reference genes are flagged for removal.

  # These genes should have no reference overlap
  gtf_no_overlap <- gtf[gtf$class_code %in% c("u", "x", "i", "y"), ]
  no_overlap_tx <- gtf_no_overlap$transcript_id

  # Create Granges with reference transcript exons using Ensembl canonical transcript
  gene_tx_refs <-
    gtf_reference[which(gsub(".*-", "", gtf_reference$transcript_name) == "201"), ]$transcript_id
  gene_tx_refs <- unique(gene_tx_refs)

  gene_ref_gtf <-
    gtf_reference[which(gtf_reference$type == "exon" &
      gtf_reference$transcript_id %in% gene_tx_refs), ]

  # Grab novel exons
  gtf_novel_exons <-
    GenomicRanges::makeGRangesFromDataFrame(gtf[which(gtf$type == "exon"), ],
      keep.extra.columns = T
    )

  gtf_no_overlap_exons <-
    GenomicRanges::makeGRangesFromDataFrame(
      gtf[which(gtf$type == "exon" &
        gtf$transcript_id %in% no_overlap_tx), ],
      keep.extra.columns = T
    )

  # Calculate overlapping exons
  no_overlap <-
    GenomicRanges::findOverlaps(
      query = gtf_no_overlap_exons,
      subject = gene_ref_gtf,
      type = "any"
    )
  if (!(length(no_overlap) == 0)) {
    message("WARNING: Overlap in unannotated genes detected")
  }

  overlap <- GenomicRanges::findOverlaps(
    query = gtf_novel_exons,
    subject = gene_ref_gtf,
    type = "any"
  )
  overlap_df <- data.frame(
    gtf_novel_exons[
      from(overlap),
      c(
        "transcript_id",
        "type",
        "class_code",
        "cmp_ref",
        "ref_gene_id"
      )
    ],
    gene_ref_gtf[to(overlap), c("gene_id", "gene_name")]
  )

  # Calculate overlapping ref genes
  split_by_tx <- split(overlap_df, overlap_df$transcript_id)

  check <- bind_rows(lapply(split_by_tx, function(x) {
    # Count per transcript the number of unique gene IDs with overlapping
    # exons. Label transcripts based on this.

    single_check <- ifelse(length(unique(x$gene_id)) > 1, "multi_gene", "single_gene")
    df <- data.frame(transcript_id = unique(x$transcript_id), single_check)

    return(df)
  }))

  check_no_overlap <- unique(subset(
    check,
    check$single_check == "multi_gene"
  )$transcript_id)

  return(check_no_overlap)
}

###########################################################################

filter_i_class <- function(gtf_df, reference_granges) {
  # Input:
  # reference_granges = GRanges object of the reference GTF
  # gtf_novel_df = data frame or data table containing new transcripts
  #
  # Searches for same-strand overlap with known genes for i class
  # transcripts.
  tx_i <- subset(gtf_df, type == "transcript" & class_code == "i")$transcript_id
  # gtf_df[which(gtf_df$type == "transcript" &
  #                gtf_df$class_code == "i")]$transcript_id
  gtf_i <- subset(gtf_df, transcript_id %in% tx_i)
  # Convert the entire gtf_i data frame to a GRanges object once
  gtf_i_granges <- makeGRangesFromDataFrame(gtf_i, keep.extra.columns = TRUE)

  # Function to check overlaps
  check_overlap <- function(tx_id) {
    # Subset GRanges object for the current transcript_id
    x <- gtf_i_granges[gtf_i_granges$transcript_id == tx_id]

    # Find overlaps
    overlap <- findOverlaps(query = x, subject = reference_granges, type = "any")

    # Check if there is any overlap
    if (length(overlap) > 0) {
      return(tx_id)
    }
  }

  # Get unique transcript IDs
  unique_tx_ids <- unique(gtf_i$transcript_id)

  # Apply the function to each transcript_id and unlist results
  i_no_pass <- unlist(lapply(unique_tx_ids, check_overlap))

  # Return the results
  return(i_no_pass)
}

###########################################################################

annotate_overlap <- function(gtf, gtf_refseq_basename, x_name) {
  # Input:
  # gtf_novel = Previous Granges output
  # gtf_refseq_basename = name of custom refseq gtf that holds XR and NR transcripts respectively
  # x_name = character vector of the name of the new column and the name of
  #          the RefSeq GFF that contains the transcripts
  #
  # Annotates transcripts found in gtf_novel with any overlap in the X
  # granges object. Currently used to annotate XR_### and NR_### transcripts
  # as these might not be applicable for neo-antigen detection.

  x <-
    rtracklayer::import(paste(gtf_refseq_basename, x_name, "gff", sep = "."))

  # Extract unique seqnames
  seqnames_list <- unique(seqnames(x))
  valid_chromosomes <- c(1:22, "X", "Y")

  # Check if seqnames are Ensembl (numbers or X, Y)
  if (any(valid_chromosomes %in% seqnames_list)) {
    x <- x[seqnames(x) %in% valid_chromosomes]
    # Check if seqnames are UCSC (chr1, chrX, etc.)
  } else if (any(grepl("^chr[0-9XY]+$", seqnames_list))) {
    seqlevels(x) <- gsub("^chr", "", seqlevels(x))
    x <- x[seqnames(x) %in% valid_chromosomes]
    # Check if seqnames are RefSeq (NC_ prefix)
  } else if (any(grepl("^NC_", seqnames_list))) {
    refseq_to_ensembl <- c(
      "NC_000001.11" = "1",
      "NC_000002.12" = "2",
      "NC_000003.12" = "3",
      "NC_000004.12" = "4",
      "NC_000005.10" = "5",
      "NC_000006.12" = "6",
      "NC_000007.14" = "7",
      "NC_000008.11" = "8",
      "NC_000009.12" = "9",
      "NC_000010.11" = "10",
      "NC_000011.10" = "11",
      "NC_000012.12" = "12",
      "NC_000013.11" = "13",
      "NC_000014.9"  = "14",
      "NC_000015.10" = "15",
      "NC_000016.10" = "16",
      "NC_000017.11" = "17",
      "NC_000018.10" = "18",
      "NC_000019.10" = "19",
      "NC_000020.11" = "20",
      "NC_000021.9"  = "21",
      "NC_000022.11" = "22",
      "NC_000023.11" = "X",
      "NC_000024.10" = "Y"
    )
    # Convert RefSeq seqnames to Ensembl equivalents
    x <- x[seqnames(x) %in% names(refseq_to_ensembl)]
    x <- keepSeqlevels(x, value = names(refseq_to_ensembl))
    seqlevels(x) <- refseq_to_ensembl[seqlevels(x)]
  } else {
    stop("Unknown seqname format. Please check the input.")
  }

  x <- subset(x, x$type == "exon")

  gtf_novel_gr <-
    GenomicRanges::makeGRangesFromDataFrame(gtf, keep.extra.columns = T)

  gtf_novel_exons <-
    subset(gtf_novel_gr, gtf_novel_gr$type == "exon")

  x_overlap <-
    GenomicRanges::findOverlaps(
      query = x,
      subject = gtf_novel_exons,
      type = "any"
    )

  novel_x_hits <-
    as.data.frame(unique(gtf_novel_exons[subjectHits(x_overlap)]))
  novel_x_hits <- unique(novel_x_hits[, "transcript_id"])
  elementMetadata(gtf_novel_gr)[[paste0(x_name, "_overlap")]] <-
    ifelse(gtf_novel_gr$transcript_id %in% novel_x_hits,
      paste0(x_name, "_hit"),
      "none"
    )
  return(gtf_novel_gr)
}

###########################################################################

rename_stringtie_transcripts <- function(gtf_novel_df) {
  # Input:
  # gtf_novel_df = data frame or data table containing new transcripts
  #
  # Renames genes without a reference gene in a similar fashion as
  # we name our ORFs using chromosome, start and stop position. Linked
  # transcripts will have the same gene name

  gtf_by_gene <- split(gtf_novel_df, gtf_novel_df$gene_id)

  gtf_renamed <-
    suppressWarnings(bind_rows(lapply(gtf_by_gene, function(x) {
      if (is.na(x$ref_gene_id)[1]) {
        chr <- unique(x$seqnames)
        start <- min(x$start)
        end <- max(x$end)
        x$gene_name <- paste0(chr, ":", start, "-", end)
      }

      return(x)
    })))
}

###########################################################################

read_tracking_file <- function(tracking_file) {
  # Read the tracking file
  df <- read.table(tracking_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

  # Extract the first 4 columns
  df_cleaned <- df[, 1:3] %>%
    tidyr::separate(V3,
      into = c("ref_gene_id_annotated", "ref_transcript_id_annotated"),
      sep = "\\|"
    )
  colnames(df_cleaned) <- c("TCONS", "xloc", "ref_gene_id", "ref_transcript_id")

  # Extract TPMs from the remaining columns
  df_tpm <- df[, 5:ncol(df)] %>%
    apply(1, function(row) {
      # Split each column by "|", extract the third-to-last value (the TPM value), or return NA if not found
      sapply(row, function(x) {
        if (x == "-" || is.na(x)) {
          return(NA)
        } # Handle missing values
        parts <- unlist(strsplit(x, "\\|"))
        if (length(parts) >= 3) {
          return(parts[length(parts) - 2])
        } else {
          return(NA)
        } # Handle cases with fewer than 3 fields
      })
    }) %>%
    t() %>% # Transpose the result to match the original data structure
    as.data.frame()

  # Add sample names (q1, q2, ...) as column names for the TPMs
  colnames(df_tpm) <- paste0("TPM_q", 1:ncol(df_tpm))

  # Combine the first 4 columns and the extracted TPMs into a final data frame
  final_df <- cbind(df_cleaned, df_tpm)

  # View the result
  return(final_df)
}
