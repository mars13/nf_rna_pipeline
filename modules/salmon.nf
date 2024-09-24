process salmon_index {
    label "salmon"

    input:
    path transcriptome // Existing transcriptome file

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index -k 31
    """
}

process salmon_quasi {
    label "salmon"
    publishDir "${outdir}/salmon", mode: 'copy', pattern: "${sample_id}/*"

    input:
    tuple(val(sample_id), path(reads))
    val paired_end
    path salmon_index
    val outdir

    output:
    path "${sample_id}/quant.sf", emit: quant
    path "${sample_id}/*"

    script:
        if (paired_end == true){
            // Set salmon quant paired end input arguments
            quant_input = """-1 "${reads[0]}" -2 "${reads[1]}" """
        } else {
            // Set salmon quant single end input arguments
            quant_input = """-i "${salmon_index}" """
        }

        """
        salmon quant \
        --libType "A" \
        --validateMappings \
        --gcBias \
        --quiet \
        --numGibbsSamples 30 \
        --threads $task.cpus \
        -i "${salmon_index}" \
        ${quant_input} \
        --output "${sample_id}"
        """
}

process salmon_bam {
    label "salmon"
    publishDir "${outdir}/salmon", mode: 'copy', pattern: "${sample_id}/*"

    input:
    tuple(val(sample_id), path(bam))
    path transcriptome
    val outdir

    output:
    path "${sample_id}/*"

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
    publishDir "${outdir}/salmon", mode: 'copy', pattern: "${prefix}*"

    input:
    val quant_paths
    path gtf
    val prefix
    val outdir

    output:
    path "${prefix}*"

    script:
    """
    salmon_cohort_tables.R \
    ${quant_paths} \
    ${gtf} \
    ${prefix}
    """
}

process featurecounts {
    label "featurecounts"
    publishDir "${outdir}/featurecounts", mode: 'copy'

    input:
    tuple(val(sample_id), path(bam))
    val reference_gtf
    val outdir

    output:
    path "${sample_id}/*"

    /*
    optional parameters to be considered:
    -s : perform strand-specific read counting
    -t : specify feature type in GTF annotation (exon by default)
    -g : specify attribute type in GTF annotation (gene_id by default)
    */
    //TODO: Add strand-specific read counting

    script:
    """
    mkdir -p ${sample_id}
    featureCounts \
    -p \
    -t gene \
    -a ${reference_gtf} \
    -o ${sample_id}/featurecounts_result.out \
    -T $task.cpus \
    ${bam}
    """

}