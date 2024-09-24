// Define process for checking strandedness
process checkStrand {
    label "qc"
    publishDir "${params.outdir}/check_strandedness", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val(paired_end)

    output:
    env(strand_info, emit: strand)
    path("**") // used to copy all process output to the output dir

    script:
    if (paired_end == true){
        """
        check_strandedness \
            -g ${params.reference_gtf} \
            -n 1000000 \
            -r1 ${reads[0]} \
            -r2 ${reads[1]} \
            -k "${params.kallisto_index}" >> "${sample_id}.txt"
        strandedness=\$(tail -n 1 ${sample_id}.txt | awk 'NF>1{print \$NF}')
        strand_info=\$(printf "%s\t%s\n" "${sample_id}" "\$strandedness")
        """
    } else {
        println "checkStrand does not currently support single end reads"
        //TODO: the github page says it does support this 
        """
        check_strandedness \
            -g ${params.reference_gtf} \
            -n 1000000 \
            -r1 ${reads[0]} \
            -k "${params.kallisto_index}" >> "${sample_id}.txt"
        strandedness=\$(tail -n 1 ${sample_id}.txt | awk 'NF>1{print \$NF}')
        strand_info=\$(printf "%s\t%s\n" "${sample_id}" "\$strandedness")
        """
    }
}