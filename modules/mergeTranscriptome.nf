// Define process for GTF merging
process mergeGTF {
    label "assembly"
    publishDir "${outdir}/gffcompare", mode: 'copy'

    input:
    path gtf_list
    val masked_fasta
    val output_basename
    val outdir

    output:
    path "${output_basename}/${output_basename}*"

    when:
    gtf_list.exists()

    script:
    """
    mkdir -p ${output_basename}
    gffcompare \
        -V \
        -r "${reference_gtf}" \
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
    val reference_gtf
    val refseq_gtf
    path gtf_novel
    path gtf_tracking
    val min_occurrence
    val output_basename
    val scripts_dir
    val outdir

    output:
    path "**${output_basename}*"

    script:
    """
    #mkdir -p plots
    filter_annotate.R \
    "${reference_gtf}" \
    "${refseq_gtf}" \
    "${gtf_novel}" \
    "${gtf_tracking}" \
    "${min_occurrence}" \
    "${output_basename}_novel_filtered" \
    "${scripts_dir}"
    """
}

// TODO FIX Create custom annotation for RiboseQC and ORFquant OR do it in Ribo-seq pipeline
process customAnotation {
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
}

