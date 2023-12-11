

// Define process for stringtie
process stringtie {
    label "assembly"
    cpus 4 //chsange to 4 after testing
    time '24h' //change to 24 after testing
    memory '24 GB' //change to 24 after testing

    input: 
        val(strand)
        tuple val(sample_id), path(bam)
        val(chr)
    
    output:
        publishDir "${params.outDir}/stringtie", mode: 'copy'
        path("${sample_id}/${sample_id}.gtf")


    script:
    chromosome_exclusion = chr ? "-x ${chr}" : "" 

    """
        stringtie ${bam[0]} \
            -G ${params.referenceGTF} \
            --${strand} \
            -M 0.45 \
            -a 9 \
            -m 50 \
            -f 0.05 \
            -g 40 \
            -j 2 \
            -s 99999 \
            -p 4 ${chromosome_exclusion} \
            -o "${sample_id}/${sample_id}.gtf"
    """
}
