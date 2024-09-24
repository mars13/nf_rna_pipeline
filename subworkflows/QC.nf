include { fastp } from '../modules/fastp'
include { checkStrand } from '../modules/strandedness'

// Groovy function to check if strand type exists and to set the strand type if found
def getStrandtype(strand) {
    if (strand ==~ /.*RF\/fr-firststrand/) {
        strand_type = tuple("rf", 1)
    } else if (strand ==~ /.*FR\/fr-secondstrand/) {
        strand_type = tuple("fr", 2)
    } else {
        strand_type = tuple("unstranded", 0)
        println "WARNING: Data is unstranded"
       // exit 1
    }
}

//Run fastp on the input reads and obatains their strandedness
workflow QC {
    take:
    reads      // Input fastq files
    paired_end // Bool, true if paired end data
    outdir     // Path to output dir

    main:
    // Run fastp and sets the output fastq file to variable trimmed_reads
    trimmed_reads = fastp(reads, paired_end, outdir).fastq_files
    
     // Run strandedness
    strand = checkStrand(reads, paired_end).strand
    // Create a file containing the strand of all files
    checkStrand.out
            .strand
            .collectFile(
            name: 'strandedness_all.txt',
            storeDir: "${outdir}/check_strandedness/",
            newLine: true, sort: true)

    strandedness = strand.map{ it -> getStrandtype(it) }

    emit:
    strandedness  // File, strandedness of the data
    trimmed_reads // Tuple, sample and trimmed fastq files

}