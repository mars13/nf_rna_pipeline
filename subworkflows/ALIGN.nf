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
        samtools(starAlign.out.bam)
        samtools.out
        .map({key, file ->
            tuple( key,
                    file.findAll({ it =~ ~/.*\.Aligned\.sortedByCoord\.out\.bam$/ })
                )
        }).set{ bam }

    emit:
        bam = bam

}