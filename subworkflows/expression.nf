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

process salmon_quant {
    label "expression"

    cpus 12
    time '48h'
    memory '64 GB'

    input:
    tuple(val(sample_id), path(reads))
    val(paired_end)
    path(salmon_index)

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

workflow expression {
    take:
        reads
        paired_end
        transcriptome

    main:
        salmon_index(transcriptome)
        salmon_quant(reads, paired_end, salmon_index.out)
}