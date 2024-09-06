include { indexLength; STAR; samtools } from '../modules/starAlign'

workflow ALIGN {
    take:
    reads
    paired_end
    outdir

    main:
    STAR(reads,
        paired_end,
        indexLength(reads),
        outdir,
        params.reference_gtf,
        params.star_index_basedir)
    
    bam = samtools(STAR.out.bam, outdir)
        .map({key, file ->
        tuple( key,
        file.findAll({ it =~ ~/.*\.Aligned\.sortedByCoord\.out\.bam$/ })
        )})
    bam.view()
    emit:
    bam

}