
// Define process for quality control (TrimGalore)
process trimGalore {
    label "alignment"
    cpus 4 
    time '24h' 
    memory '24 GB' 
    

    input: 
        tuple val(sample_id), path(reads)
        val(pairedEnd)
    
    output:
        publishDir "${params.outDir}", mode: 'copy'
        tuple val("${sample_id}"), path("trimgalore/${sample_id}/*")


    script:
    def r1_fname = "trimgalore/${sample_id}/${reads[0].getBaseName(2)}_val_1.fq.gz"
    def r1_trimmed = "trimgalore/${sample_id}/${reads[0].name.replaceFirst('R1', 'R1_trimmed')}"

    if (pairedEnd == true){
        """   
        mkdir -p trimgalore 
        trim_galore \
            "${reads[0]}" "${reads[1]}" \
            --cores 2 \
            --paired \
            --gzip \
            --fastqc \
            --fastqc_args "--outdir trimgalore/${sample_id}/" \
            --output_dir "trimgalore/${sample_id}"

        mv "${r1_fname}" "${r1_trimmed}"
        mv "${r1_fname.replaceFirst("_R1_", "_R2_").replaceFirst("_val_1.fq.gz", "_val_2.fq.gz")}" "${r1_trimmed.replaceFirst("_R1_", "_R2_")}"
        tot_reads=\$(zcat "${reads[0]}" | echo \$((`wc -l`/4)))
        trimmed_reads=\$(zcat "${r1_trimmed}" | echo \$((`wc -l`/4)))
        trimmed_percentage=`awk -vn=248 "BEGIN{print(\${trimmed_reads}/\${tot_reads}*100)}"`
            
        printf '%s\t%s\t%s\t%s\n' "${sample_id}" "Trimmed" \$trimmed_reads \$trimmed_percentage > "trimgalore/${sample_id}/${sample_id}_trim_stats.txt"
            
        """
    } else {
        println "trimgalore does not currently support single end reads"
    }

}

// Define process for checking strandedness
process checkStrand {
    label "alignment"

    cpus 2 
    time '24h' 
    memory '10 GB' 

    input: 
        tuple val(sample_id), path(reads)
        val(pairedEnd)
    
    output:
        publishDir "${params.outDir}/check_strandedness", mode: 'copy'
        tuple env(strand_info), path("**")

    script:
    if (pairedEnd == true){
        """
        check_strandedness \
            -g ${params.referenceGTF} \
            -n 1000000 \
            -r1 ${reads[0]} \
            -r2 ${reads[1]} \
            -k "${params.kallistoIndex}" >> "${sample_id}.txt" 
        strandedness=\$(tail -n 1 ${sample_id}.txt | awk 'NF>1{print \$NF}')
        strand_info=\$(printf "%s\t%s\n" "${sample_id}" "\$strandedness")
        """
    } else {
        println "checkStrand does not currently support single end reads"
    }

}


process indexLength{
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
        publishDir "${params.outDir}/STAR/", mode: 'copy'
        path("${sample_id}/${sample_id}.*")

    script:
    // STAR params defined as concatenated strings singe triple quote multiline declaration generates newline
    def star_params = "--readFilesCommand zcat" +
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

// Define process for generating mapping statistics with Samtools
process samtools {
    label "alignment"

    cpus 16 
    time '24h' 
    memory '48 GB'

    // Define input/output as needed
    // ...

    script:
    """
    # Your code for Samtools
    """
}
