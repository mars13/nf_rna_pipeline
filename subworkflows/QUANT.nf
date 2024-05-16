include { salmon_index; salmon_quasi; salmon_bam } from '../modules/salmon'

workflow EXPRESSION {
    take:
        reads
        bam
        mode
        paired_end
        transcriptome
        outdir

    main:
        if (mode ~= "sq") {
            salmon_index(transcriptome)
            salmon_quasi(reads, paired_end, salmon_index.out, outdir)
        } else if  (mode ~= "sa" ) {
            salmon_bam(bam, transcriptome, outdir)
        } else if (mode ~= "fc") {
            //TODO Gene expression with featureCounts if bam available
            println "FeatureCounts not yet implemented"
        }

    //TODO create salmon counts and expression tables per cohort

}