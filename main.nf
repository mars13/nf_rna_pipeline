include { paramsHelp } from 'plugin/nf-schema'
include { printHeader } from "./modules/helperFunctions.nf"
include { RNASEQ } from "./workflows/RNASEQ"

workflow {

    printHeader()

    if (params.help) {

        log.info paramsHelp("nextflow run mars13/nf_rna_pipeline -c params.config")

    } else {

        RNASEQ()
    }

}

