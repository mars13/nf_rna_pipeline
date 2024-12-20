// Define process for checking strandedness
process checkStrand {
    label "qc"
    publishDir "${outdir}/check_strandedness", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) // Tuple containing sample id and read paths
    val paired_end                    // Bool, true if data is paired end
    val kallisto_index                // Path to kallisto index dir
    val reference_gtf                 // Path to input reference gtf file
    val outdir                        // Path to output directory

    output:
    env strand_info, emit: strand
    path "${sample_id}.txt", emit: strand_file

    script:
    if (paired_end == true){
        """
        check_strandedness \
            -g ${params.reference_gtf} \
            -n 1000000 \
            -r1 ${reads[0]} \
            -r2 ${reads[1]} \
            -k "${kallisto_index}" >> "${sample_id}.txt"
        strandedness=\$(tail -n 1 ${sample_id}.txt | awk 'NF>1{print \$NF}')
        strand_info=\$(printf "%s\t%s\n" "${sample_id}" "\$strandedness")
        rm -r stranded_test_*
        """
    } else {
        println "checkStrand does not currently support single end reads"
        //TODO: the github page says it does support this 
        """
        check_strandedness \
            -g ${reference_gtf} \
            -n 1000000 \
            -r1 ${reads[0]} \
            -k "${kallisto_index}" >> "${sample_id}.txt"
        strandedness=\$(tail -n 1 ${sample_id}.txt | awk 'NF>1{print \$NF}')
        strand_info=\$(printf "%s\t%s\n" "${sample_id}" "\$strandedness")
        rm -r stranded_test_*
        """
    }
}