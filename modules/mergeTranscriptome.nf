// Define process for GTF merging
process mergeGTF {
    label "assembly"
    publishDir "${outdir}/gffcompare", mode: 'copy'

    input:
    path gtf_list         // List of previously created gtf files to be merged
    val masked_fasta      // Path to input reference fasta file
    val reference_gtf     // Path to input reference gtf file
    val output_basename   // Val containing the id/name given to the output files
    val outdir            // Path to output directory

    output:
    path "${output_basename}/${output_basename}*"
    path "${output_basename}/${output_basename}.combined.gtf", emit: merged_gtf
    path "${output_basename}/${output_basename}.tracking", emit: tracking

    when:
    gtf_list.exists()

    script:
    """
    mkdir -p ${output_basename}
    gffcompare \
        -V \
        -r ${reference_gtf} \
        -s ${masked_fasta} \
        -o "${output_basename}/${output_basename}" \
        -i "${gtf_list}"
    """
}


// Define process for transcript filtering and annotation
process filterAnnotate {
    label "assembly"
    publishDir "${outdir}/customannotation/", mode: 'copy'

    input:
    val reference_gtf   // Path to the input reference gtf file
    val refseq_gtf      // Path to input refseq gtf file
    path gtf_novel      // Path to the merged gtf file
    path gtf_tracking   // Path to the tracking file created by the merge step
    val min_occurrence  // Val contatining the minimum occurence of transcripts for filtering
    val min_tpm         // Val containing the minium tpm of transcripts for filtering
    val output_basename // Val containing output basename
    val scripts_dir     // Path location of input R scripts
    val outdir          // Path to output directory

    output:
    path "${output_basename}_novel_filtered.gtf", emit: gtf
    path "${output_basename}_novel_filtered.log"
    path "${output_basename}_novel_filtered.tsv"

    script:
    def refseq_arg = refseq_gtf ? "\"${refseq_gtf}\"" : ""

    """
    filter_annotate.R \
    "${reference_gtf}" \
    "${gtf_novel}" \
    "${gtf_tracking}" \
    "${min_occurrence}" \
    "${min_tpm}" \
    "${output_basename}_novel_filtered" \
    "${scripts_dir}" \
    ${refseq_arg}
    """
}


// TODO FIX Create custom annotation for RiboseQC and ORFquant OR do it in Ribo-seq pipeline
// TODO: discuss this function
/* process customAnotation {
    clusterOptions '--mem=10G --cpus-per-task=2 --gres=tmpspace:50G --time=24:00:00'
    containerOptions '/hpc:/hpc",${TMPDIR}:${TMPDIR} --env "LC_CTYPE=en_US.UTF-8'
    publishDir "${outdir}/customannotation/", mode: 'copy'

    input:
    val twobit
    val merged_gtf
    val output_basename
    val outdir          
    
    output:
    path "custom_annotation/"

    script:
    """
    Rscript orfquant_custom_annotation.R \
    ${params.twobit} \
    ${merged_gtf}\
    ${params.output_basename}/" \
    ${params.output_basename} \
    ${package_install_loc}
    """
} */

// Creates a fasta file of the transcript sequence using the reference fasta file and the transcriptome gtf
process transcriptome_fasta {
    label "gffread"
    publishDir "${outdir}/stringtie", mode: 'copy'

    input:
    val merged_filtered_gtf // Merged and filtered transcriptome file
    val masked_fasta        // Path to input reference fasta file
    val outdir              // Path to output directory

    output:
    file "stringtie_transcriptome.fa"

    script:
    """
    gffread -w stringtie_transcriptome.fa -g ${masked_fasta} ${merged_filtered_gtf}
    """
}