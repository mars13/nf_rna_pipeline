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

// Create custom annotation for RiboseQC and ORFquant
process customAnotation {
    label "assembly"
    cpus 2
    time '24h' 
    memory '10 GB'
    clusterOptions '--gres=tmpspace:G'


    input: 
        tuple val(sample_id), path(reads)
        val(paired_end)
    
    output:
        publishDir "${params.outdir}", mode: 'copy'
        tuple val("${sample_id}"), path("trimgalore/${sample_id}/*")


    script:
    """
    """
}

