// Define process for checking strandedness
process checkStrand {
    label "qc"
    publishDir "${outdir}/check_strandedness", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) // Tuple containing sample id and read paths
    val kallisto_index                // Path to kallisto index dir
    val reference_gtf                 // Path to input reference gtf file
    val outdir                        // Path to output directory

    output:
    tuple val("${sample_id}"), env(strand_info), emit: strand
    path "${sample_id}.txt", emit: strand_file

    script:
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
}