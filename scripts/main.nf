include { rnaseq_alignment } from './subworkflows/rnaseq_alignment'
include { transcriptome_assembly } from './subworkflows/transcriptome_assembly'


workflow {
   // Create reands input channel 
    reads =  Channel
            .fromFilePairs( params.readsPath, size: params.pairedEnd ? 2 : 1 )
            .ifEmpty { "Could not find any reads matching pattern: ${params.readsPath}" }
    // assign R1 and resolve symlinks
    r1 = reads.map { keys, files -> new File(files[0].toString()).canonicalPath } 
    // assign R2 and resolve symlinks
    r2 = params.pairedEnd ? reads.map { keys, files -> new File(files[1].toString()).canonicalPath } : null
    
    //write R1 file
    r1.collectFile(
        name: 'r1_files.txt',
        storeDir: "${params.projectFolder}/documentation/",
        newLine: true, sort: true)

    //write R2 file
    if (params.pairedEnd) {   
        r2.collectFile(
            name: 'r2_files.txt',
            storeDir: "${params.projectFolder}/documentation/", 
            newLine: true, sort: true)
    }

    //write samples file
    reads
    .map { keys, files -> keys }
    .set { sample_ids }

    sample_ids
    .collectFile(
        name: 'sample_ids.txt',
        storeDir: "${params.projectFolder}/documentation/",
        newLine: true, sort: true)

    // Step 01: Align reads
    rnaseq_alignment(reads, params.pairedEnd, params.qc, params.align, params.outDir)
    
    strand = rnaseq_alignment.out.strand
    bam = rnaseq_alignment.out.bam

    // Step 02: Aseemble transcriptome

    //Check if strand info available
    if (!params.qc && !params.strandInfo) {
        // Look for strandedness summary file
        strand_summary = "${params.outDir}/check_strandedness/strandedness_all.txt"
        if (file(strand_summary).isEmpty()) {
            println "No strand information found in ${strand_summary}"
            println "When running with `qc = false`, you can define your strandedness info file as `strandInfo=[path_to_file]` in `params.config`"
            println "The file is expected to be plain text file (tab separated, no headers) with sample_id and strand information (RF/fr-firststrand, FR/fr-secondstrand, unstranded)."

            exit 1
        } else {
            Channel.fromPath(strand_summary)
            .splitText(){ it.trim() }
            .set { strand }
        }
    } 

    //Check if strand info available
        if (!params.align && !params.bamFiles) {
            // Look for strandedness summary file
            bamFiles = "${params.outDir}/star/**/*.Aligned.sortedByCoord.out.bam"
            if (file(bamFiles).isEmpty()) {
                println "No .bam found in ${bamFiles}"
                println "When running with `align = false`, you can define your bam location as `bamFiles=[path_to_file]` in `params.config`"
                println "Please, remember to index your bams. "
                exit 1
            } else {
                Channel.fromFilePairs("${params.outDir}/star/**/*.Aligned.sortedByCoord.out.bam", size: 1, checkIfExists: true)
                .set { bam }

            }
        } 
    //Dev prints
    strand.view()
    bam.view()

    //Run transcriptome assembly
    transcriptome_assembly(strand, bam)
}

