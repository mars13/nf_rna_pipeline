include { fastp } from '../modules/fastp'
include { checkStrand } from '../modules/strandedness'

//Run fastp on the input reads and obatains their strandedness
workflow QC {
    take:
    reads      // input fastq files
    paired_end // bool, true if paired end data
    outdir     // path to output dir

    main:
    // Run fastp and sets the output fastq file to variable trimmed_reads
    trimmed_reads = fastp(reads, paired_end, outdir).fastq_files
    
     // Run strandedness
    strandedness = checkStrand(reads, paired_end).strand
    // Create a file containing the strand of all files
    checkStrand.out
            .strand
            .collectFile(
            name: 'strandedness_all.txt',
            storeDir: "${outdir}/check_strandedness/",
            newLine: true, sort: true)

    emit:
    strandedness  // File, strandedness of the data
    trimmed_reads // Tuple, sample and trimmed fastq files

}