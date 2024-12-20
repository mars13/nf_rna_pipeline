include { salmon_index; salmon_quasi; salmon_bam } from '../modules/salmon'

workflow IMMUNE {

    take:
        bams
        tpms

    main:
    //TODO implement CIBERSORTx, IPASS and arcasHLA
    println("Workflow skeleton")


}