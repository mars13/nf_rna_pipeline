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

    input:
    tuple(val(sample_id), path(reads))
    val paired_end
    path salmon_index
    path outdir

    output:
    publishDir "${params.outdir}/salmon", mode: 'copy', pattern: "${sample_id}/*"
    path("${sample_id}/*")

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
    path outdir

    output:
    publishDir "${params.outdir}/salmon", mode: 'copy', pattern: "${sample_id}/*"
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
    path base_dir
    path gtf
    val prefix
    val testing

    output:
    publishDir "${params.outdir}/salmon", mode: 'copy', pattern: "${prefix}*"
    path "${prefix}*"

    script:
    """
    salmon_cohort_tables.R \
    ${base_dir} \
    ${params.reference_gtf} \
    ${prefix}
    """
}