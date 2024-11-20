include { indexLength; STAR; samtools } from '../modules/starAlign'

/*
* Runs star mapper with given input files, sort resulting bam files with samtools
*/
workflow ALIGN {
    take:
    reads              // Tuple, sample_id and location of input read files
    paired_end         // Bool, True if data is paired end
    reference_gtf      // Location of reference gtf
    star_index_basedir // Location of star index
    outdir             // Location of output directory

    main:
    index_length = indexLength(reads)

    // Map reads to genome using star
    STAR(reads,
        paired_end,
        index_length,
        reference_gtf,
        star_index_basedir,
        outdir)
    
    // Sort star output using samtools
    samtools(STAR.out.bam, outdir)
    bam = samtools.out.sorted_bam
    bam_list = samtools.out.sorted_bam.toList()

    emit:
    bam // Tuple of sample id and sorted bam file
    bam_list // All bam tuples combined into one big tuple
}
