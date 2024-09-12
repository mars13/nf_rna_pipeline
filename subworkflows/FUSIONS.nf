include { starAlignChimeric; runArriba } from '../modules/arriba'

workflow FUSIONS {
    take:
    reads
    paired_end
    arriba_reference
    outdir

    main:
    // Location to the Arriba reference
    reference = Channel.fromPath("${arriba_reference}**")

    // Using fromFilePairs to make sure a fa_fai index exists for reference fasta
    fa = Channel.fromFilePairs("${arriba_reference}**{.fa,.fa.fai}", checkIfExists: true).flatten().filter(~/.*\.fa/).first()
    gtf = reference.filter(~/.*\.gtf/).first()

    // Take optional files and passing empty file if empty
    blacklist = reference.filter(~/.*blacklist.*/).ifEmpty{ file("EMPTY_BLACKLIST")}.first()
    whitelist = reference.filter(~/.*known_fusions.*/).ifEmpty{ file("EMTPY_WHITELIST")}.first()
    protein_domains = reference.filter(~/.*protein_domains.*/).ifEmpty{ file("EMPTY_DOMAINS")}.first()

    // TODO: test this function with actual data
    wgs = params.wgs_sv ? Channel.fromFilePairs(params.wgs_sv).ifEmpty{ file("EMPTY_WGS")} : Channel.fromPath("EMPTY_WGS").first()

    // Run star mapper
    starAlignChimeric(reads, paired_end, arriba_reference, outdir)
    fusion_bam = starAlignChimeric.out.bam

    // Run Arriba
    runArriba(fusion_bam, fa, gtf, blacklist, whitelist, protein_domains, wgs, outdir)

}