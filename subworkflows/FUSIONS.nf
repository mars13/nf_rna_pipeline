include { starAlignChimeric; runArriba} from '../modules/arriba'

workflow FUSIONS {
    take:
    reads             // Tuple: (sample_id, [fastq1, fastq2])
    vcf               // Tuple: (sample_id, vcf) or (sample_id, null) if no vcf available
    paired_end        // Channel bool, paired end or not
    arriba_reference  // Path arriba reference location
    outdir            // Path to output dir

    main:
    // Location to the Arriba reference
    reference = channel.fromPath("${arriba_reference}**")

    // Check if both fasta and index exists for Arriba reference 
    fa = channel.fromPath("${arriba_reference}**.fa", checkIfExists: true).first()
    channel.fromPath("${arriba_reference}**.fa.fai", checkIfExists: true).first()

    // Get gtf from the reference folder
    gtf = reference.filter(~/.*\.gtf/).first()

    // Take optional files and passing empty file if empty
    blacklist = reference.filter(~/.*blacklist.*/).ifEmpty{ file("EMPTY_BLACKLIST")}.first()
    whitelist = reference.filter(~/.*known_fusions.*/).ifEmpty{ file("EMTPY_WHITELIST")}.first()
    protein_domains = reference.filter(~/.*protein_domains.*/).ifEmpty{ file("EMPTY_DOMAINS")}.first()

    // Run star mapper
    starAlignChimeric(reads, paired_end, arriba_reference, outdir)
    fusion_bam = starAlignChimeric.out.bam

    arriba_input = fusion_bam
                    .join(vcf)

    // Run Arriba
    runArriba(arriba_input, 
        fa, 
        gtf,
        blacklist, 
        whitelist, 
        protein_domains, 
        outdir
    )
}