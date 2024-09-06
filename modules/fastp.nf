// Define process for trimming and quality control (fastp)
process fastp {
    label "qc"

    publishDir "${outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val(paired_end)
    val(outdir)

    output:
    tuple val("${sample_id}"), path("fastp/${sample_id}/*_trimmed.fastq.gz"), emit: fastq_files // tuple of sample id, and all fastq files
    path "fastp/${sample_id}/${sample_id}_fastp_report.json" //output the json files for publishDir

    script:
    def r1_out = "fastp/${sample_id}/${reads[0].getBaseName(2)}"
    def json_out = "fastp/${sample_id}/${sample_id}_fastp_report.json"

    if (paired_end == true){

        def r2_out = "fastp/${sample_id}/${reads[1].getBaseName(2)}"

        """
        mkdir -p fastp/${sample_id}
        fastp \
        -i "${reads[0]}" \
        -I "${reads[1]}" \
        -o "${r1_out}_trimmed.fastq.gz" \
        -O "${r2_out}_trimmed.fastq.gz" \
        --detect_adapter_for_pe \
        --verbose \
        --json ${json_out} \
        --report_title '${sample_id} fastp report' \
        --thread $task.cpus
        """
    } else {
        """
        mkdir -p fastp/${sample_id}
        fastp \
        -i "${reads[0]}" \
        -o "${r1_out}_trimmed.fastq.gz" \
        --verbose \
        --json ${json_out} \
        --report_title '${sample_id} fastp report' \
        --thread $task.cpus
        """
    }
}