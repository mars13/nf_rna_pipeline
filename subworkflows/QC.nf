include { fastp } from '../modules/fastp'
include { checkStrand } from '../modules/strandedness'

// Groovy function to check if strand type exists and to set the strand type if found
def getStrandtype(strand_tuple) {
    def (sample_id, strand) = strand_tuple
    def strand_type = null

    if (strand ==~ /.*RF\/fr-firststrand/) {
        strand_type = tuple(sample_id, "rf", 1)
    } else if (strand ==~ /.*FR\/fr-secondstrand/) {
        strand_type = tuple(sample_id, "fr", 2)
    } else {
        strand_type = tuple(sample_id, "unstranded", 0)
        // Optionally warn or exit here
        // println "WARNING: Data for ${sample_id} is unstranded"
        // exit 1
    }

    return strand_type
}

//Run fastp on the input reads and obatains their strandedness
workflow QC {
    take:
    reads          // Input fastq files
    paired_end     // Bool, true if paired end data
    kallisto_index // Path to the kallisto index file
    reference_gtf  // Path to the reference gtf 
    outdir         // Path to output dir

    main:
    // Run fastp and sets the output fastq file to variable trimmed_reads
    trimmed_reads = fastp(reads, paired_end, outdir).fastq_files
    fastp_json = fastp.out.fastp_json
    
    // Run strandedness
    strand = checkStrand(reads, paired_end, kallisto_index, reference_gtf, outdir).strand
    
    // Create a file containing the strand of all files
    checkStrand.out
            .strand
            .collectFile(
            name: 'strandedness_all.txt',
            storeDir: "${outdir}/check_strandedness/",
            newLine: true, sort: true)

    strandedness = strand.map{ it -> getStrandtype(it) }
    strand_file = checkStrand.out.strand_file

    emit:
    strandedness  // File, strandedness of the data
    trimmed_reads // Tuple, sample and trimmed fastq files
    fastp_json    // Fastp stats
    strand_file   // Text file containing the strandeness of all input samples
}