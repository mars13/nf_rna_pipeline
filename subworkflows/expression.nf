include { salmon_index; salmon_quasi} from "../modules/salmon"
//TODO make process for salmon_quasi and salmon_bam
workflow EXPRESSION {
    take:
        reads
        bam
        mode
        paired_end
        transcriptome

    main:
        if (mode ~= "sq") {
            salmon_index(transcriptome)
            salmon_quant(reads, paired_end, salmon_index.out, outdir)
        } else if  (mode ~= "sa" ) {
            salmon_bam(bam, outdir)
        } else if (mode ~= "fc") {
            println "FeatureCounts not yet implemented"
        }


        //Gene expression with featureCounts if bam available



}