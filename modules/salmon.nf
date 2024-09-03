process salmon_index {
    label "salmon"

    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index -k 31
    """
}

//Salmon genomic alignment
process salmon_quasi {
    label "salmon"
    publishDir "${params.outdir}/salmon", mode: 'copy', pattern: "${sample_id}/*"

    input:
    tuple(val(sample_id), path(reads))
    val paired_end
    path salmon_index
    path outdir

    output:
    path("${sample_id}/*", emit: salmon_output_dir)
    path("${sample_id}/quant.sf", emit: quant)

    script:
        if (paired_end == true){
            """
            salmon quant \
            --libType "A" \
            --validateMappings \
            --gcBias \
            --quiet \
            --numGibbsSamples 30 \
            --threads $task.cpus \
            -i "${salmon_index}" \
            -1 "${reads[0]}" \
            -2 "${reads[1]}" \
            --output "${sample_id}"
            """
        } else {
            """
            salmon quant \
            --libType "A" \
            --validateMappings \
            --gcBias \
            --quiet \
            --numGibbsSamples 30 \
            --threads $task.cpus \
            -i "${salmon_index}" \
            -r "${reads}" \
            --output "${sample_id}"
            """
        }
}


process salmon_bam {
    label "salmon"

    input:
    tuple(val(sample_id), path(bam))
    path transcriptome
    val outdir

    output:
    publishDir "${outdir}/salmon", mode: 'copy', pattern: "${sample_id}/*"
    path("${sample_id}/*")

    script:
    """
    salmon quant \
    -t "${transcriptome}" \
    -a "${bam}" \
    --libType "A" \
    --gcBias \
    --quiet \
    --numGibbsSamples 30 \
    --threads $task.cpus \
    --output "${sample_id}"
    """
}


process salmon_tables {
    label "salmon_tables"

    input:
    val quant_paths
    path gtf
    val prefix

    output:
    publishDir "${params.outdir}/salmon", mode: 'copy', pattern: "${prefix}*"
    path "${prefix}*"

    script:
    """
    salmon_cohort_tables.R \
    ${quant_paths} \
    ${params.reference_gtf} \
    ${prefix}
    """
}