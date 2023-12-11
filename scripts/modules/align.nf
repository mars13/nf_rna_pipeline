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
process starAlign {
    label "alignment"

    cpus 16 
    time '24h' 
    memory '100 GB' 

    // Define input/output as needed
    input: 
        tuple val(sample_id), path(reads)
        val(pairedEnd)
        val(usedIndex)

    output:
        publishDir "${params.outDir}/star/", mode: 'copy', pattern: "${sample_id}/${sample_id}*.out"    //fix to add .out.tab
        path("${sample_id}/${sample_id}.*"), emit: files
        tuple val("${sample_id}"), path("${sample_id}/${sample_id}.*.bam"), emit: bam

    script:
    // STAR params defined as concatenated strings singe triple quote multiline declaration generates newline
    def star_params = "--readFilesCommand zcat " +
                      "--twopassMode Basic --runThreadN 16 --runDirPerm All_RWX " +
                      "--outFilterType BySJout --outSAMunmapped Within " +
                      "--outSAMattributes NH HI AS nM NM MD jM jI MC ch " +
                      "--outSAMstrandField intronMotif --outSAMtype BAM Unsorted " +
                      "--outFilterMismatchNmax 6 --outTmpKeep None " +
                      "--alignSJoverhangMin 10 --outFilterMultimapNmax 10 " +
                      "--outFilterScoreMinOverLread 0.75"
    
    if (pairedEnd == true){
        """
        # Use STAR for mapping the reads
        STAR --genomeDir "${params.starIndexBasedir}/${usedIndex}nt" \
        --sjdbGTFfile ${params.referenceGTF} \
        --readFilesIn "${reads[0]}" "${reads[1]}" \
        --outSAMattrRGline ID:${sample_id} LB:${sample_id} PL:IllUMINA SM:${sample_id} \
        --outFileNamePrefix "${sample_id}/${sample_id}." \
        ${star_params}
        """
    } else {
        println "STAR does not currently support single end reads"
    }
}

// Define process for getting mapping stats, sorted bam and .bai with Samtools
process samtools {
    label "alignment"

    cpus 16 
    time '24h' 
    memory '48 GB'
    clusterOptions '--gres=tmpspace:100G'

    input: 
        tuple val(sample_id), path(bam)

    output:
        publishDir "${params.outDir}/star/", mode: 'copy'
        tuple val(sample_id), path("${sample_id}/${sample_id}*")

    script:
    def new_bam= "${bam.name.replaceFirst('.Aligned.out.bam', '.Aligned.sortedByCoord.out.bam')}"


    """
    mkdir -p ${sample_id}
    mkdir -p tmp/
    # Sort BAM
    samtools sort \
    -@ 8 \
    -l 9 \
    -o "${sample_id}/${new_bam}" \
    -T "tmp/" \
    "${bam}"

    rm -r tmp/

    # Create mapping statistics with samtools
    samtools stats -@ 8 "${sample_id}/${new_bam}" > "${sample_id}/${sample_id}_stats.txt"

    # Index the bam with samtools
    samtools index -@ 8 "${sample_id}/${new_bam}" 
    """
}
