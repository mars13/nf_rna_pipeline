include { printStartLogs; check_params} from "./modules/helperFunctions.nf"
include { rnaseq_qc } from './subworkflows/rnaseq_qc'
include { rnaseq_alignment } from './subworkflows/rnaseq_alignment'
include { transcriptome_assembly } from './subworkflows/transcriptome_assembly'
include { fusion_calling } from './subworkflows/fusion_calling'
include { expression } from './subworkflows/expression'


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
        fusion_calling(star_input, params.paired_end, params.arriba_reference)
    }

    // Step 05: Expression
    if (params.expression) {
        expression(star_input, params.paired_end, params.reference_transcriptome)
    }
}

