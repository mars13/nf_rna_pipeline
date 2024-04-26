include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'
include { printHeader; checkInputFiles} from "./modules/helperFunctions.nf"
include { rnaseq_qc } from './subworkflows/rnaseq_qc'
include { rnaseq_alignment } from './subworkflows/rnaseq_alignment'
include { transcriptome_assembly } from './subworkflows/transcriptome_assembly'
include { fusion_calling } from './subworkflows/fusion_calling'
include { expression } from './subworkflows/expression'


workflow {
    // Initialise workflow
    if (params.help) {
        log.info paramsHelp("nextflow run mars13/nf_rna_pipeline -c params.config")
        exit 0
    }

    printHeader()

    // Validate input parameters
    validateParameters()

    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    // Validate files
    checkInputFiles()

    //Create input channel from samplesheet
    ch_input = Channel.fromSamplesheet("samplesheet")
    //ch_input.view()
    ch_input
        .filter { meta, filename_1, filename_2 ->
            meta.filetype == "fastq" && meta.sequence == "rna"
        }
        .map { meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ paired_end:false ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ paired_end:true ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .set{ ch_rna_reads }

        ch_rna_reads.view()
    
    /*
    DEPRECATED IN FAVOUR OF SAMPLESHEET
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
    */

    // Step 01: QC
    if (params.qc) {
        rnaseq_qc(reads, params.paired_end, params.outdir)
        star_input = rnaseq_qc.out.trimmed_reads
        strand = rnaseq_qc.out.strand
    } else {
        // Look for trimmed reads at the usual location
        default_trimmed_reads = "${params.outdir}/trimgalore/**/*{R1,R2}_trimmed*.{fastq.gz,fq.gz}"
        if (file(default_trimmed_reads).isEmpty()) {
            star_input = reads
            println  "No trimmed reads found in path: ${default_trimmed_reads}, using ${params.reads_path}"
        } else {
            star_input = Channel
            .fromFilePairs("${default_trimmed_reads}", size: params.paired_end ? 2 : 1)
        }

        // Look for strandedness summary file
        default_strand_files = "${params.outdir}/check_strandedness/strandedness_all.txt"
        if (params.strand_info) {
            Channel.fromPath(strand_info, checkIfExists: true)
            .splitText(){ it.trim() }
            .set { strand }
        } else if (!file(default_strand_files).isEmpty()) {
            Channel.fromPath(default_strand_files, checkIfExists: true)
            .splitText(){ it.trim() }
            .set { strand }
        } else {
            strand = null
        }
    }

    // Step 02: Align
    if (params.align){
        rnaseq_alignment(star_input, params.paired_end, params.outdir)
        bam = rnaseq_alignment.out.bam
    } else {
         // Look for alignment files
        default_bams = "${params.outdir}/star/**/*.Aligned.sortedByCoord.out.bam"
        if (params.bam_files) {
            Channel.fromFilePairs(params.bam_files, size: 1, checkIfExists: true)
            .set { bam }
        } else if (!file(default_bams).isEmpty()) {
            Channel.fromFilePairs(default_bams, size: 1, checkIfExists: true)
            .set { bam }
        } else {
            bam = null
        }
    }

    // Step 03: Assemble transcriptome
    if (params.assembly || params.merge) {
        transcriptome_assembly(strand, bam)
    }

    // Step 04: Fusion calling
    if (params.fusions) {
        fusion_calling(star_input,
                        params.paired_end,
                        params.arriba_reference)
    }

    // Step 05: Expression
    if (params.expression) {
            expression(star_input,
                        bam,
                        params.expression_mode,
                        params.paired_end,
                        params.reference_transcriptome)
    }
}

