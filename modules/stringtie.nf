

// Define process for stringtie
process stringtie {
    label "assembly"

    input:
        val(strand)
        tuple val(sample_id), path(bam)
        val(chr)

    output:
        publishDir "${params.outdir}/stringtie", mode: 'copy'
        path("${sample_id}/${sample_id}.gtf")

    script:
    def chromosome_exclusion = chr ? "-x ${chr}" : ""


    """
        stringtie ${bam[0]} \
            -G ${params.reference_gtf} \
            --${strand} \
            -M 0.45 \
            -a 9 \
            -m 50 \
            -f 0.05 \
            -g 40 \
            -j 2 \
            -s 99999 \
            -p ${task.cpus} \
            ${chromosome_exclusion} \
            -o "${sample_id}/${sample_id}.gtf"
    """
}
