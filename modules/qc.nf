// Define process for quality control (TrimGalore)
//TODO update to fastp
process trimGalore {
    label "alignment"
    cpus 4 
    time '24h' 
    memory '24 GB' 
    

    input: 
        tuple val(sample_id), path(reads)
        val(paired_end)
    
    output:
        publishDir "${params.outdir}", mode: 'copy'
        tuple val("${sample_id}"), path("trimgalore/${sample_id}/*")


    script:
    def r1_fname = "trimgalore/${sample_id}/${reads[0].getBaseName(2)}_val_1.fq.gz"
    def r1_trimmed = "trimgalore/${sample_id}/${reads[0].name.replaceFirst('R1', 'R1_trimmed')}"
    
    if (paired_end == true){
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
        //TODO add single end commands

    }

}


process fastp {

//fastp commands


}

// Define process for checking strandedness
process checkStrand {
    label "alignment"

    cpus 2 
    time '24h' 
    memory '10 GB' 

    input: 
        tuple val(sample_id), path(reads)
        val(paired_end)
    
    output:
        publishDir "${params.outdir}/check_strandedness", mode: 'copy'
        tuple env(strand_info), path("**")

    script:
    if (paired_end == true){
        """
        check_strandedness \
            -g ${params.reference_gtf} \
            -n 1000000 \
            -r1 ${reads[0]} \
            -r2 ${reads[1]} \
            -k "${params.kallisto_index}" >> "${sample_id}.txt" 
        strandedness=\$(tail -n 1 ${sample_id}.txt | awk 'NF>1{print \$NF}')
        strand_info=\$(printf "%s\t%s\n" "${sample_id}" "\$strandedness")
        """
    } else {
        println "checkStrand does not currently support single end reads"
        //TODO add single end commands
    }

}
