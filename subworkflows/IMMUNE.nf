include { salmon_index; salmon_quasi } from '../modules/salmon'

workflow IMMUNE {

    take:
        _bams
        _tpms

    main:
    //TODO implement CIBERSORTx, IPASS and arcasHLA
    println("Workflow skeleton")


}