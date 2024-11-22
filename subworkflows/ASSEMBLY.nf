include { aletsch } from '../modules/aletsch'
include { stringtie } from '../modules/stringtie'
include { mergeGTF; filterAnnotate; customAnotation; transcriptome_fasta } from '../modules/mergeTranscriptome'

workflow ASSEMBLY {
    take:
    strand
    bam
    bam_list
    sample_gtf_list
    reference_gtf
    refseq_gtf
    chr_exclusion_list
    masked_fasta
    output_basename
    min_occurrence
    min_tpm
    outdir

    main:
    // Run stringtie unless paths to precomputed individual sample GTF are provided
    if (!sample_gtf_list & params.assembly) {

        // Locate chromosome exclusion list
        if(chr_exclusion_list) {
        // Groovy command to join the chr exclusion list on ","
            chromosome_exclusion_list = file("${projectDir}/${chr_exclusion_list}").readLines().join(",")
        } else {
            chromosome_exclusion_list = null
        }

        // Run aletsch
        aletsch(bam_list,outdir)

        // Run stringtie
        stringtie(strand, bam, chromosome_exclusion_list, reference_gtf, outdir)

        // Wait untill all stringtie runs are completed
        gtf_paths = stringtie.out.collect().flatten()
                    .map { it -> it.toString() } // Change paths to strings

        // Store gtflist in outputdir
        // replaceFirst is a groovy statement to replace part of a string based on a pattern
        gtf_paths.map { it -> it.replaceFirst("${workDir}/[^/]*/[^/]*/", "${outdir}/stringtie/") } //Replace the work dir path with the output dir
        .collectFile(
            name: 'gtflist.txt',
            storeDir: "${outdir}/stringtie/",
            newLine: true, sort: true )

        // Store gtflist to workdir
        gtf_list = gtf_paths.collectFile(
            name: 'gtflist.txt',
            newLine: true, sort: true )

    }

    // Merges the gtf files created by stringtie
    if (params.merge) {

        if (sample_gtf_list) {
            gtf_list = Channel.fromPath("${sample_gtf_list}")
            gtf_list
            .splitText(){ it.trim() }
            .take(1)
            .ifEmpty { error "Could not find sample GTF files in: ${sample_gtf_list}" }
        }

        // Run merge process
        mergeGTF(gtf_list, masked_fasta, reference_gtf, output_basename, outdir)

        gtf_novel = mergeGTF.out.flatten().filter(~/.*\.combined\.gtf/)
        gtf_tracking = mergeGTF.out.flatten().filter(~/.*\.tracking/)

        scripts_dir = Channel.fromPath("${workflow.projectDir}/bin/")

        // Run filter annotate r script
        filterAnnotate(reference_gtf,
                        refseq_gtf,
                        gtf_novel,
                        gtf_tracking,
                        min_occurrence,
                        min_tpm,
                        output_basename,
                        "${projectDir}/bin/",
                        outdir)

        merged_gtf = filterAnnotate.out.gtf

        //if (params.custom_annotation) {
        //   if ( params.custom_annotation ==~ /orfquant/ ) {
        //        customAnotation(merged_gtf)
        //   }
        //}


        transcriptome_fasta(merged_gtf, masked_fasta, outdir)
        stringtie_transcriptome = transcriptome_fasta.out
    }else{
        merged_gtf = null
        stringtie_transcriptome = null
    }
    emit:
    merged_gtf
    stringtie_transcriptome
}
