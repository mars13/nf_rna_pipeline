include { trimGalore; checkStrand } from '../modules/qc'

workflow rnaseq_qc {
    take:
        reads
        paired_end
        outdir


    main:
        //Run trimGalore
        trimGalore(reads, "${params.paired_end}")

        //Collect trimGalore stats (changed to within trimgalore folder)
        trimGalore.out.map{key, files -> files }
            .flatten()
            .filter(~/.*trim_stats\.txt/)
            .collectFile(
                name: 'trim_stats.txt',
                storeDir: "${params.outdir}/trimgalore/",
                newLine: false, sort: true)

        //Run strandedness
        strandedness = checkStrand(reads, "${params.paired_end}").map { keys, files -> keys }
        strandedness
            .collectFile(
                name: 'strandedness_all.txt',
                storeDir: "${params.outdir}/check_strandedness/",
                newLine: true, sort: true)

    emit:
        strand = strandedness
        trimmed_reads = trimGalore.out
                                    .map({key, file ->
                                        tuple( key,
                                            file.findAll({ it =~ /.*(?:R1|R2)_trimmed.*\.fastq\.gz$/ })
                                            )
                                        })
}