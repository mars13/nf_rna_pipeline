// Define process to assign best index length
process indexLength{
    label "alignment"

    input:
        tuple val(sample_id), path(reads)
    output:
        env used_index
    script:

    """
    # Check first 10k reads for read length for star index
    read_length=\$(gunzip -c "${reads[0]}" | \
    head -n 10000 | \
    awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count - 1}')

    # Make sure read length picks are correct
    if [[ "\${read_length}" =~ "." ]]; then
    read_length=\${read_length%.*}
    fi

    if (( \${read_length} >= 70 && \${read_length} <= 85 )); then
    used_index=74
    elif (( \${read_length} >= 86 && \${read_length} <= 111 )); then
    used_index=99
    elif (( \${read_length} >= 112 && \${read_length} <= 137 )); then
    used_index=124
    elif (( \${read_length} >= 138 && \${read_length} <= 163 )); then
    used_index=149
    else
    used_index=99
    fi
    """


}

// Define process for alignment with STAR
process STAR {
    label "alignment"

    input: 
        tuple val(sample_id), path(reads)
        val(paired_end)
        val(usedIndex)
        path(outdir)
        val(reference_gtf)
        val(star_index_basedir)

    output:
        publishDir "${params.outdir}/star/", mode: 'copy'
        path("${sample_id}/${sample_id}.*"), emit: files
        tuple val("${sample_id}"), path("${sample_id}/${sample_id}.*.bam"), emit: bam

    script:
    // STAR params defined as concatenated strings singe triple quote multiline declaration generates newline
    def star_params = "--readFilesCommand zcat " +
                      "--twopassMode Basic --runDirPerm All_RWX " +
                      "--outFilterType BySJout --outSAMunmapped Within " +
                      "--outSAMattributes NH HI AS nM NM MD jM jI MC ch " +
                      "--outSAMstrandField intronMotif --outSAMtype BAM Unsorted " +
                      "--outFilterMismatchNmax 6 --outTmpKeep None " +
                      "--alignSJoverhangMin 10 --outFilterMultimapNmax 10 " +
                      "--outFilterScoreMinOverLread 0.75"

    if (paired_end == true){
        """
        # Use STAR for mapping the reads
        STAR --genomeDir "${star_index_basedir}/${usedIndex}nt" \
        --sjdbGTFfile ${reference_gtf} \
        --readFilesIn "${reads[0]}" "${reads[1]}" \
        --outSAMattrRGline ID:${sample_id} LB:${sample_id} PL:IllUMINA SM:${sample_id} \
        --outFileNamePrefix "${sample_id}/${sample_id}." \
        --runThreadN $task.cpus ${star_params}
        """
    } else {
        println "STAR does not currently support single end reads"
        //TODO add single end commands

    }
}

// Define process for getting mapping stats, sorted bam and .bai with Samtools
process samtools {
    label "alignment"

    input: 
        tuple val(sample_id), path(bam)

    output:
        publishDir "${params.outdir}/star/", mode: 'copy'
        tuple val(sample_id), path("${sample_id}/${sample_id}*")

    script:
    def new_bam = "${bam.name.replaceFirst('.Aligned.out.bam', '.Aligned.sortedByCoord.out.bam')}"


    """
    mkdir -p ${sample_id}
    mkdir -p tmp/
    # Sort BAM
    samtools sort \
    -@ $task.cpus \
    -l 9 \
    -o "${sample_id}/${new_bam}" \
    -T "tmp/" \
    "${bam}"

    rm -r tmp/

    # Create mapping statistics with samtools
    samtools stats -@ $task.cpus "${sample_id}/${new_bam}" > "${sample_id}/${sample_id}_stats.txt"

    # Index the bam with samtools
    samtools index -@ $task.cpus "${sample_id}/${new_bam}"
    """
}
