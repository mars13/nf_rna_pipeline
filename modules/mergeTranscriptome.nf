// Define process for GTF merging
process mergeGTF {
    label "assembly"

    input:
        path(gtf_list)

    output:
        publishDir "${params.outdir}/gffcompare", mode: 'copy'
        path("${params.output_basename}/${params.output_basename}*")

    when:
        gtf_list.exists()

    script:
    """
    mkdir -p ${params.output_basename}
    gffcompare \
        -V \
        -r "${params.reference_gtf}" \
        -s ${params.masked_fasta} \
        -o "${params.output_basename}/${params.output_basename}" \
        -i "${gtf_list}"
    """
}

// Define process for transcript filtering and annotation
process filterAnnotate {
    label "assembly"

    input:
        val(reference_gtf)
        val(refseq_gtf)
        path(gtf_novel)
        path(gtf_tracking)
        val(min_occurrence)
        val(output_basename)
        val(scripts_dir)


    output:
        publishDir "${params.outdir}/customannotation/", mode: 'copy'
        path("**${output_basename}*")


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

    input:
        merged_gtf
    output:
        publishDir "${params.outdir}/customannotation/", mode: 'copy'
        path("custom_annotation/")


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

