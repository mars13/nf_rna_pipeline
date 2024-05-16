include { starAlignChimeric; runArriba } from '../modules/arriba'

workflow FUSIONS {
    take:
        reads
        paired_end
        arriba_reference


    main:
        reference = Channel.fromPath("${arriba_reference}**")

        //Using fromFilePairs to make sure a fa_fai index exists for reference fasta
        fa = Channel.fromFilePairs("${arriba_reference}**{.fa,.fa.fai}", checkIfExists: true).flatten().filter(~/.*\.fa/)
        gtf = reference.filter(~/.*\.gtf/)

        //Take optional files and assing empty file if empty
        blacklist = reference.filter(~/.*blacklist.*/).ifEmpty{ file("EMPTY_BLACKLIST")}
        whitelist = reference.filter(~/.*known_fusions.*/).ifEmpty{ file("EMTPY_WHITELIST")}
        protein_domains = reference.filter(~/.*protein_domains.*/).ifEmpty{ file("EMPTY_DOMAINS")}

        //DEV change to channel.fromFilePairs to get sampleId key
        wgs = params.wgs_sv ? Channel.fromFilePairs(params.wgs_sv).ifEmpty{ file("EMPTY_WGS")} : Channel.fromPath("EMPTY_WGS")

        starAlignChimeric(reads, paired_end, arriba_reference)
        //DEV fusion_bams = Channel.fromFilePairs("/hpc/pmc_vanheesch/projects/Marina/dev/nf_rna_pipeline/test/analysis/arriba/star/*/*.Aligned.out.bam", size: 1)
        fusion_bams = starAlignChimeric.out.bam

        fusion_bams
            .combine(fa)
            .combine(gtf)
            .combine(blacklist)
            .combine(whitelist)
            .combine(protein_domains)
            .combine(wgs, by: 0) //hint for when wgs is avail: .combine(wgs, by: 0)
            .set{ arriba_input }

        runArriba(arriba_input)

}