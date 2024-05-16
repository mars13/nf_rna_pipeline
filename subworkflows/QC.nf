include { fastp } from '../modules/fastp'
include { checkStrand } from '../modules/strandedness'

workflow QC {
    take:
        reads
        paired_end
        outdir

    main:
        //Run fastp
        fastp(reads, paired_end)

        //Collect fastp stats (changed to within fastp folder)
        fastp.out.map{key, files -> files }
            .flatten()
            .filter(~/.*trim_stats\.txt/)
            .collectFile(
                name: 'trim_stats.txt',
                storeDir: "${params.outdir}/fastp/",
                newLine: false, sort: true)

        //Run strandedness
        strandedness = checkStrand(reads, paired_end).map { keys, files -> keys }
        strandedness
            .collectFile(
                name: 'strandedness_all.txt',
                storeDir: "${params.outdir}/check_strandedness/",
                newLine: true, sort: true)

    emit:
        strand = null
        strand = strandedness
        trimmed_reads = fastp.out
                            .map({key, file ->
                                tuple( key,
                                    file.findAll({ it =~ /.*(?:R1|R2)_trimmed.*\.fastq\.gz$/ })
                                    )
                                })
}