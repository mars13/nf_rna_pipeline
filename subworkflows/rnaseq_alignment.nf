include { indexLength; starAlign; samtools } from '../modules/align'

workflow rnaseq_alignment {
    take:
        reads
        paired_end
        outdir

    main:
        starAlign(reads, "${params.paired_end}", indexLength(reads))
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