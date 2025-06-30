include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { checkInputFiles } from "../modules/helperFunctions.nf"
include { QC } from '../subworkflows/QC'
include { ALIGN } from '../subworkflows/ALIGN'
include { ASSEMBLY } from '../subworkflows/ASSEMBLY'
include { FUSIONS } from '../subworkflows/FUSIONS'
include { EXPRESSION } from '../subworkflows/EXPRESSION'
include { MULTIQC } from '../modules/multiqc.nf'

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

    multiqc_files = Channel.empty()

    //TODO handle bamfiles from samplesheet
    //TODO handle WGS vcfs from samplesheet
    
    /*
    * Step 01: QC
    */
    if (params.qc) {
        QC(reads, 
        paired_end,
        params.kallisto_index,
        params.reference_gtf,
        params.outdir)
        
        star_input = QC.out.trimmed_reads
        strand = QC.out.strandedness
        
        multiqc_files = multiqc_files.mix(QC.out.fastp_json)
        multiqc_files = multiqc_files.mix(QC.out.strand_file)

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
            Channel.fromPath(params.strand_info, checkIfExists: true)
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

    /*
    *Step 02: Align
    */
    if (params.align){
        ALIGN(star_input, paired_end, params.reference_gtf, params.star_index_basedir, params.outdir)
        bam = ALIGN.out.bam
        multiqc_files = multiqc_files.mix(ALIGN.out.star_log)
        multiqc_files = multiqc_files.mix(ALIGN.out.samtools_stats)

    } else {
        // Look for alignment files in case star has been run previously
        default_bams = "${params.outdir}/star/**/*.Aligned.sortedByCoord.out.bam"
        if (params.bam_files) {
            bam = Channel.fromFilePairs(params.bam_files, size: 1, checkIfExists: true)
        } else if (!file(default_bams).isEmpty()) { // TODO: do we need this function?
            bam = Channel.fromFilePairs(default_bams, size: 1, checkIfExists: true)
        } else {
            bam = null
         }
    }

    /*
    * Step 03: Assemble transcriptome
    */
    if (params.assembly || params.merge) {
        if (paired_end) {
            // Join the strand info with the bam file to prevent sample mixing
            stringtie_input = strand.join(bam)
            .map { row ->
                    def sample_id = row[0]    // Sample Id
                    def strand_label = row[1] // Strand value
                    def file_path = row[3]    // BAM file path
                    [sample_id, strand_label, file_path]
                }

            ASSEMBLY(stringtie_input,
            params.sample_gtf_list,
            params.reference_gtf,
            params.refseq_gtf,
            params.chr_exclusion_list,
            params.masked_fasta,
            params.output_basename,
            params.min_occurrence,
            params.min_tpm,
            params.outdir)

            assembled_gtf = ASSEMBLY.out.merged_filtered_gtf
            assembled_fasta = ASSEMBLY.out.assembled_transcriptome_fasta

        } else {
            assembled_gtf = null
            println "Transcriptome assembly not supported for single stranded data."
            exit 1
        }
    } else{
        assembled_gtf = null
        assembled_fasta = null
    }
    
    /*
    * Step 04: Fusion calling
    * TODO: could we combine mapping for fusions and for txome assembly in one?
    */
    if (params.fusions) {
        FUSIONS(star_input,
            paired_end,
            params.arriba_reference,
            params.outdir)
    }

    /*
    * Step 05: Expression
    */
    if (params.expression) {
        // Join the strand info with the bam file to prevent sample mixing
        featurecounts_input = strand.join(bam)
        .map { row ->
                def sample_id = row[0]    // Sample Id
                def strand_code = row[2]  // Strand code value
                def file_path = row[3]    // BAM file path
                [sample_id, strand_code, file_path]
            }

        EXPRESSION(star_input,
            featurecounts_input,
            assembled_gtf,
            assembled_fasta,
            params.expression_mode,
            paired_end,
            params.reference_transcriptome,
            params.reference_gtf,
            params.output_basename,
            params.outdir)
    }

    /*
    * Step 06: Immune landscape
    */
    
    /*
    * Step 07: MultiQC
    */
    multiqc_config_file = file("${projectDir}/${params.multiqc_config}")
    MULTIQC (multiqc_files.collect(), multiqc_config_file, params.outdir) 
}

