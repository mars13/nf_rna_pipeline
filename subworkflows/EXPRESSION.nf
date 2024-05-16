include { salmon_index; salmon_quant} from "../modules/salmon"
//TODO make process for salmon_quasi and salmon_bam
workflow EXPRESSION {
    take:
        reads
        bam
        mode
        paired_end
        transcriptome

    main:
        if (mode == "sq" || !bam) {
            salmon_index(transcriptome)
            salmon_quant(reads, paired_end, salmon_index.out)
        } else {

        }


        //Gene expression with featureCounts if bam available



}