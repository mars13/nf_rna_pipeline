process starAlignChimeric {
    label "fusions"

    cpus 16
    time '48h'
    memory '64 GB'

    input:
        tuple val(sample_id), path(reads)
        val(paired_end)
        path reference

    output:
        publishDir "${params.outdir}/arriba/${sample_id}/star", mode: 'copy', pattern: "${sample_id}*.out*"
        path("${sample_id}.*"), emit: files
        tuple val("${sample_id}"), path("${sample_id}.*.bam"), emit: bam

    script:
    if (paired_end == true){
        """
        STAR --runThreadN 16 \
        --genomeDir ${reference}/STAR_index* \
        --genomeLoad NoSharedMemory \
        --readFilesIn "${reads[0]}" "${reads[1]}" --readFilesCommand zcat \
        --outSAMtype BAM Unsorted --outSAMunmapped Within \
        --outFileNamePrefix "${sample_id}." \
        --outBAMcompression 0 --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 \
        --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 \
        --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3 --chimMultimapNmax 50
        """
    } else {
        error  "`--fusions` not supported when `--paired_end false`"
    }
}

process runArriba{
    label "fusions"
    clusterOptions '--mem=104G --cpus-per-task=16 --gres=tmpspace:10G --time=12:00:00'

    input:
        tuple val(sample_id),
              path(bam),
              path(fa),
              path(gtf),
              path(blacklist),
              path(whitelist),
              path(protein_domains),
              path(wgs)

    output:
        publishDir "${params.outdir}/arriba/${sample_id}", mode: 'copy', pattern: "${sample_id}_fusions*"
        path("${sample_id}_fusions*"), emit: fusions

    script:

    def blacklist_options = blacklist.toString().contains("EMPTY") ?  "-f blacklist" : "-b ${blacklist}"
    def whitelist_options = whitelist.toString().contains("EMPTY") ? "" : "-k ${whitelist} -t ${whitelist}"
    def domains_options = protein_domains.toString().contains("EMPTY") ? "" : "-p ${protein_domains}"
    def wgs_options = wgs.toString().contains("EMPTY") ?  "" : "-d ${wgs}"

    """
    mkdir -p ${sample_id}/
    /arriba_v2.4.0/arriba -x ${bam} \
    -o ${sample_id}_fusions.tsv \
    -O ${sample_id}_fusions.discarded.tsv \
    -a ${fa} \
    -g ${gtf} \
    ${blacklist_options} ${whitelist_options} ${domains_options} ${wgs_options}
    """
}


