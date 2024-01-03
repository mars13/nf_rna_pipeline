include { printStartLogs; check_params} from "./modules/helperFunctions.nf"
include { rnaseq_alignment } from './subworkflows/rnaseq_alignment'
include { transcriptome_assembly } from './subworkflows/transcriptome_assembly'

workflow {
    // Initialise workflow
    printStartLogs()
    check_params()

   // Create reads input channel 
    reads =  Channel
            .fromFilePairs(params.reads_path, size: params.paired_end ? 2 : 1, checkIfExists: true)
    // assign R1 and resolve symlinks
    r1 = reads.map { keys, files -> new File(files[0].toString()).canonicalPath } 
    // assign R2 and resolve symlinks
    r2 = params.paired_end ? reads.map { keys, files -> new File(files[1].toString()).canonicalPath } : null
    
    //write R1 file
    r1.collectFile(
        name: 'r1_files.txt',
        storeDir: "${params.project_folder}/documentation/",
        newLine: true, sort: true)

    //write R2 file
    if (params.paired_end) {   
        r2.collectFile(
            name: 'r2_files.txt',
            storeDir: "${params.project_folder}/documentation/", 
            newLine: true, sort: true)
    }

    //write samples file
    reads
    .map { keys, files -> keys }
    .set { sample_ids }

    sample_ids
    .collectFile(
        name: 'sample_ids.txt',
        storeDir: "${params.project_folder}/documentation/",
        newLine: true, sort: true)

    // Step 01: Align reads
    rnaseq_alignment(reads, params.paired_end, params.qc, params.align, params.outdir)
    
    strand = rnaseq_alignment.out.strand
    bam = rnaseq_alignment.out.bam

    // Step 02: Assemble transcriptome
    //Check if pre-computed sample GTFs are available 
    if ( params.sample_gtf_list) {
        strand = null
        bam = null
    } else {
        //Check if qc was executed
        if (!params.qc) {
            // Look for strandedness summary file
            if (params.strand_info) {
                Channel.fromPath(strand_info, checkIfExists: true)
                .splitText(){ it.trim() }
                .set { strand }
            } else {
                strand_summary = "${params.outdir}/check_strandedness/strandedness_all.txt"
                Channel.fromPath(strand_summary, checkIfExists: true)
                .splitText(){ it.trim() }
                .set { strand }
            }
        } 

        //Check if alignment available 
            if (!params.align) {
                // Look for alignment files 
                if (params.bam_files) {
                    Channel.fromFilePairs(params.bam_files, size: 1, checkIfExists: true)
                    .set { bam }
                    } else {
                    Channel.fromFilePairs("${params.outdir}/star/**/*.Aligned.sortedByCoord.out.bam", size: 1, checkIfExists: true)
                    .set { bam }
                }
            } 
    }

    //Run transcriptome assembly
    transcriptome_assembly(strand, bam)
}

