process salmon_index {
    label "expression"

    cpus 12
    time '48h'
    memory '64 GB'

    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index -k 31
    """
}

process salmon_quasi {
    label "expression"

    cpus 12
    time '48h'
    memory '64 GB'

    input:
    tuple(val(sample_id), path(reads))
    val paired_end
    path salmon_index
    val expression_mode
    path outdir

    output:
    publishDir "${outdir}/salmon", mode: 'copy', pattern: "${sample_id}/*"
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
    label "expression"

    cpus 12
    time '48h'
    memory '64 GB'

    input:
    path bam
    path transcriptome
    path outdir

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