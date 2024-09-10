// Define process for trimming and quality control (fastp)
process fastp {
    label "qc"
    publishDir "${outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) // Tuple of sample id and fastq read files
    val paired_end                    // Bool, true if data is paired end
    val outdir                        // Path to outputdir

    output:
    tuple val("${sample_id}"), path("fastp/${sample_id}/*_trimmed.fastq.gz"), emit: fastq_files // Tuple of sample id, and all fastq files
    path "fastp/${sample_id}/${sample_id}_fastp_report.json" // Output the json files for publishDir

    script:
    def r1_out = "fastp/${sample_id}/${reads[0].getBaseName(2)}"
    def json_out = "fastp/${sample_id}/${sample_id}_fastp_report.json"
    
    if (paired_end == true){
        // Set paired end second output for fastp command if paired end
        def r2_out = "fastp/${sample_id}/${reads[1].getBaseName(2)}"
        fastp_input_2 = """ -I "${reads[1]}" """
        fastp_output_2 = """ -O "${r2_out}_trimmed.fastq.gz" --detect_adapter_for_pe """
    } else {
        // Set paired end second output for fastp command to empty if single end
        fastp_input_2 = ""
        fastp_output_2 = ""
    }

    // Run fastp with given parameters
    """
    mkdir -p fastp/${sample_id}
    fastp \
    -i "${reads[0]}" \
    ${fastp_input_2} \
    -o "${r1_out}_trimmed.fastq.gz" \
    ${fastp_output_2} \
    --verbose \
    --json ${json_out} \
    --report_title '${sample_id} fastp report' \
    --thread $task.cpus
    """

}