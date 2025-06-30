include { stringtie } from '../modules/stringtie'
include { mergeGTF; filterAnnotate; transcriptome_fasta } from '../modules/mergeTranscriptome'

workflow ASSEMBLY {
    take:
    stringtie_input
    sample_gtf_list    // List of previously created gtf files to be merged
    reference_gtf      // Path to input reference gtf file
    refseq_gtf         // Path to input refseq gtf file
    chr_exclusion_list // Path to chromosome exclustion list
    masked_fasta       // Path to input reference fasta file
    output_basename    // Val containing the id/name given to the output files
    min_occurrence     // Val contatining the minimum occurence of transcripts for filtering
    min_tpm            // Val containing the minium tpm of transcripts for filtering
    outdir             // Path to output directory

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


        // Run stringtie
        stringtie(stringtie_input, chromosome_exclusion_list, reference_gtf, outdir)

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

        // Load gtf list file if not null
        if (sample_gtf_list) {
            gtf_list = Channel.fromPath("${sample_gtf_list}")
            gtf_list
            .splitText(){ it.trim() }
            .take(1)
            .ifEmpty { error "Could not find sample GTF files in: ${sample_gtf_list}" }
        }

        // Run merge process
        mergeGTF(gtf_list, masked_fasta, reference_gtf, output_basename, outdir)

        gtf_novel = mergeGTF.out.merged_gtf
        gtf_tracking = mergeGTF.out.tracking

        // Set channel containing paths to r scripts
        //scripts_dir = Channel.fromPath("${workflow.projectDir}/bin/")
        
        // Run filter annotate r script
        filterAnnotate(reference_gtf,
                        refseq_gtf ?: "",
                        gtf_novel,
                        gtf_tracking,
                        min_occurrence,
                        min_tpm,
                        output_basename,
                        "${projectDir}/bin/",
                        outdir)

        merged_filtered_gtf = filterAnnotate.out.gtf

        transcriptome_fasta(merged_filtered_gtf, masked_fasta, outdir)
        assembled_transcriptome_fasta = transcriptome_fasta.out
    }else{
        merged_filtered_gtf = null
        assembled_transcriptome_fasta = null
    }
    emit:
    merged_filtered_gtf
    assembled_transcriptome_fasta
}
