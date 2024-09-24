

// Define process for stringtie
process stringtie {
    label "assembly"
    publishDir "${outdir}/stringtie", mode: 'copy'

    input:
    tuple val(stringtie), val(featurecounts) // Strandedness of input reads
    tuple val(sample_id), path(bam) // Tuple of sample id and star result bam
    val chr                         // Chromosome exclusion list 
    val reference_gtf               // Reference gtf file location
    val outdir

    output:
    path "${sample_id}/${sample_id}.gtf" // Path to the output gtf of stringtie

    script:
    // Check if reads are unstranded and exit if they are
    if (stringtie == "unstranded") {
        println "Data must be stranded"
        exit 1 
    }

    // Groovy function to add the chromosome exclusion list to the stringtie command
    def chromosome_exclusion = chr ? "-x ${chr}" : ""

    """
    stringtie ${bam[0]} \
        -G ${reference_gtf} \
        --${stringtie} \
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
