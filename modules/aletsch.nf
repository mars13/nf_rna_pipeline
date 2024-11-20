process aletsch{
    publishDir "${outdir}/aletsch/", mode: 'copy'

    input:
    val bam_list_tuples
    val outdir

    output: 
    path "aletsch_input.txt"
    path "aletsch_transcriptome.gtf"

    script:
   // Groovy function to write each tuple's second value to output parameter
   // Output file format is bam file path, bai file path, paired_end
    def write_aletsch_input = bam_list_tuples.collect { itemTuple -> 
        def bam_file = itemTuple[1]  // Access the second element of each tuple
        return "${bam_file} ${bam_file}.bai paired_end" // Write output to aletsch input format
    }.join('\n')
    
    //Groovy function to set max number of splice graphs to recommended twice the number of samples
    def max_splice_graphs = bam_list_tuples.size() * 2

    """
    echo "$write_aletsch_input" > aletsch_input.txt
    mkdir aletsch_profile
    mkdir aletsch_gtf

    aletsch --profile -i aletsch_input.txt -p aletsch_profile
    aletsch -i aletsch_input.txt -o aletsch_transcriptome.gtf -p aletsch_profile -d aletsch_gtf -t ${task.cpus} -c ${max_splice_graphs}
    """
}
