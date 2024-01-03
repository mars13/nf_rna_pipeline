include { trimGalore; checkStrand } from '../modules/qc'
include { indexLength; starAlign; samtools } from '../modules/align'

workflow rnaseq_alignment {
    take: 
        reads
        paired_end
        qc
        align
        outdir


    main:
        if (params.qc) {
            //Run trimGalore
            trimGalore(reads, "${params.paired_end}")

            //Collect trimGalore stats (changed to within trimgalore folder)
            trimGalore.out.map{key, files -> files }
            .flatten()
            .filter(~/.*trim_stats\.txt/)
            .collectFile(
                name: 'trim_stats.txt',
                storeDir: "${params.outdir}/trimgalore/",
                newLine: false, sort: true)
            
            //Run strandedness
            strandedness = checkStrand(reads, "${params.paired_end}").map { keys, files -> keys }
            strandedness
            .collectFile(
                name: 'strandedness_all.txt',
                storeDir: "${params.outdir}/check_strandedness/",
                newLine: true, sort: true)
        }   
        
        if (params.align) {
            if (params.qc){
                //Collect trimGalore outputs 
                trimGalore.out
                .map({key, file ->
                    tuple( key, 
                        file.findAll({ it =~ /.*(?:R1|R2)_trimmed.*\.fastq\.gz$/ })
                        )
                    }).set{ star_input }
            } else {
                // Look for trimmed reads at the usual location
                trimmed_reads = "${params.outdir}/trimgalore/**/*{R1,R2}_trimmed*.{fastq.gz,fq.gz}"
                if (file(trimmed_reads).isEmpty()) {
                    star_input = reads
                    println  "No trimmed reads found in path: ${trimmed_reads}, using ${params.reads_path}"

                } else {
                    star_input = Channel
                    .fromFilePairs("${trimmed_reads}", size: params.paired_end ? 2 : 1)   
                }

            }
            starAlign(star_input, "${params.paired_end}", indexLength(star_input))
            samtools(starAlign.out.bam)
            samtools.out
            .map({key, file ->
                tuple( key, 
                        file.findAll({ it =~ ~/.*\.Aligned\.sortedByCoord\.out\.bam$/ })
                    )
            }).set{ bam }
        }


    emit:
        strand = qc ? strandedness : null
        bam = align ? bam : null
        
    
}