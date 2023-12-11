// Define process for GTF merging
process mergeGTF {
    label "assembly"
    cpus 2 //change to 4
    time '24h' 
    memory '10 GB' //change 48

    input: 
        path(gtf_list)
    
    output:
        publishDir "${params.outDir}/gffcompare", mode: 'copy'
        path("${params.mergedGTFbasename}/${params.mergedGTFbasename}*")

    when:
        gtf_list.exists()

    script:
    """
    mkdir -p ${params.mergedGTFbasename}
    gffcompare \
        -V \
        -r "${params.referenceGTF}" \
        -s ${params.maskedFasta} \
        -o "${params.mergedGTFbasename}/${params.mergedGTFbasename}" \
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
        publishDir "${params.outDir}/customannotation/", mode: 'copy'
        path("**${params.mergedGTFbasename}*")


    script:
    """
    mkdir -p plots
    filter_annotate.R \
    "${params.referenceGTF}" \
    "${params.refseqGTF}" \
    "${gtf_novel}" \
    "${gtf_tracking}" \
    "${min_occurrence}" \
    "${params.mergedGTFbasename}"_novel_filtered.gtf

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
        val(pairedEnd)
    
    output:
        publishDir "${params.outDir}", mode: 'copy'
        tuple val("${sample_id}"), path("trimgalore/${sample_id}/*")


    script:
    """
    """
}

