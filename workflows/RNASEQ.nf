include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { checkInputFiles } from "../modules/helperFunctions.nf"
include { QC } from '../subworkflows/QC'
include { ALIGN } from '../subworkflows/ALIGN'
include { ASSEMBLY } from '../subworkflows/ASSEMBLY'
include { FUSIONS } from '../subworkflows/FUSIONS'
include { EXPRESSION } from '../subworkflows/EXPRESSION'

workflow RNASEQ {
  
    // Give output directory location
    //params.outdir_rnaseq = ${params.outdir}/rnaseq/
  
    // Initialise workflow

    // Validate input parameters
    validateParameters()

    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    // Validate files
    checkInputFiles()

    //Create input channel from samplesheet
    ch_input = Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))

    //Get channel for RNA reads (paired or single)
    ch_input
        .filter { meta, filename_1, filename_2 ->
            meta.filetype == "fastq" && meta.sequence == "rna"
        }
        .map { meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta + [ paired_end:false ], [ fastq_1 ] ]
                } else {
                    return [ meta + [ paired_end:true ], [ fastq_1, fastq_2 ] ]
                }
        }
        .set{ ch_fastq }

    ch_fastq
        .multiMap { meta, reads ->
            reads: tuple(meta.id, reads)
            paired_end:   meta.paired_end
        }
        .set { ch_reads }

    reads = ch_reads.reads
    paired_end = ch_reads.paired_end

    //TODO handle bamfiles from samplesheet
    //TODO handle WGS vcfs from samplesheet

    // Step 01: QC
    if (params.qc) {
        QC(reads, paired_end, params.outdir)
        star_input = QC.out.trimmed_reads
        strand = QC.out.strandedness
    } else {
        default_trimmed_reads = "${params.outdir}/**/*{R1,R2}*_trimmed.{fastq.gz,fq.gz}"
        if (file(default_trimmed_reads).isEmpty()) {
            star_input = reads
            log.info "No trimmed reads found from QC step, using input reads."
        } else {
            //log.info "Using trimmed reads in: ${default_trimmed_reads}" MISLEADING WHEN BAM IS PRESENT, TODO change or remove
            star_input = Channel
                .fromFilePairs(default_trimmed_reads, size: paired_end ? 2 : 1, checkIfExists: true)
        }

        // Look for strandedness summary file
        default_strand_files = "${params.outdir}/check_strandedness/strandedness_all.txt"
        if (params.strand_info) {
            Channel.fromPath(strand_info, checkIfExists: true)
            .splitText(){ it.trim() }
            .set { strand }

            log.info
        } else if (!file(default_strand_files).isEmpty()) {
            Channel.fromPath(default_strand_files, checkIfExists: true)
            .splitText(){ it.trim() }
            .set { strand }
        } else {
            strand = null
            log.info "No strandedness provided, setting strand to `null`. WARNING: It will conflict with assembly step."
        }
    }

    // Step 02: Align
    if (params.align){
        //TODO provide all params from main workflow
        ALIGN(star_input, paired_end, params.outdir)
        bam = ALIGN.out.bam
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
        if (paired_end) {
            //TODO provide al params from main workflow
            ASSEMBLY(strand, bam)
        } else {
            println "Transcriptome assembly not suported for single stranded data."
            exit 1
        }
    }

    // Step 04: Fusion calling
    // TODO could we combine mapping for fusions and for txome assembly in one?
    if (params.fusions) {
        FUSIONS(star_input,
            params.paired_end,
            params.arriba_reference)
    }

    // Step 05: Expression
    if (params.expression) {
        EXPRESSION(star_input,
            bam,
            params.expression_mode,
            paired_end,
            params.reference_transcriptome,
            params.reference_gtf,
            params.outdir,
            params.output_basename)
    }

    //Step 06: Immune landscape

    //Step 07: MultiQC

}

