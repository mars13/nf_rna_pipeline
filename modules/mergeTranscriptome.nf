// Define process for GTF merging
process mergeGTF {
    label "assembly"
    cpus 2 //change to 4
    time '24h' 
    memory '10 GB' //change 48

    input: 
        path(gtf_list)

    output:
        publishDir "${params.outdir}/gffcompare", mode: 'copy'
        path("${params.merged_gtf_basename}/${params.merged_gtf_basename}*")

    when:
        gtf_list.exists()

    script:
    """
    mkdir -p ${params.merged_gtf_basename}
    gffcompare \
        -V \
        -r "${params.reference_gtf}" \
        -s ${params.masked_fasta} \
        -o "${params.merged_gtf_basename}/${params.merged_gtf_basename}" \
        -i "${gtf_list}"
    """
}

// Define process for transcript filtering and annotation
process filterAnnotate {
    label "assembly"
    cpus 2
    time '24h'
    memory '10 GB' //change to 24

    input:
        path(gtf_novel)
        path(gtf_tracking)
        val(min_occurrence)


    output:
        publishDir "${params.outdir}/customannotation/", mode: 'copy'
        path("**${params.merged_gtf_basename}*")


    script:
    """
    mkdir -p plots
    filter_annotate.R \
    "${params.reference_gtf}" \
    "${params.refseq_gtf}" \
    "${gtf_novel}" \
    "${gtf_tracking}" \
    "${min_occurrence}" \
    "${params.merged_gtf_basename}"_novel_filtered.gtf

    """
}

// TODO FIX Create custom annotation for RiboseQC and ORFquant
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
    ${params.merged_gtf_basename}/" \
    ${params.merged_gtf_basename} \
    ${package_install_loc}
    """
}

