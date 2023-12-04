nextflow.enable.dsl = 2

include { trimGalore; checkStrand } from './modules/qc'
include { indexLength; starAlign; samtools } from './modules/align'

process getVersions {
    label "rnaseq"
    cpus 1
    input:
        path container_dir
    output:
        publishDir "${params.projectFolder}/documentation/", mode: 'copy'
        path "versions.txt", emit: versions
    script:
    """
    """
}


workflow {
    
    //check params

    reads =  Channel
            .fromFilePairs( params.readsPath, size: params.pairedEnd ? 2 : 1 )
            .ifEmpty { "Could not find any reads matching pattern: ${params.readsPath}" }
    // assign R1 and resolve symlinks
    r1 = reads.map { keys, files -> new File(files[0].toString()).canonicalPath } 
    // assign R2 and resolve symlinks
    r2 = params.pairedEnd ? reads.map { keys, files -> new File(files[1].toString()).canonicalPath } : null
    
    //write R1 file
    r1.collectFile(
        name: 'r1_files.txt',
        storeDir: "${params.projectFolder}/documentation/",
        newLine: true, sort: true)

    //write R2 file
    if (params.pairedEnd) {   
        r2.collectFile(
            name: 'r2_files.txt',
            storeDir: "${params.projectFolder}/documentation/", 
            newLine: true, sort: true)
    }

    //write samples file
    reads
    .map { keys, files -> keys }
    .set { sample_ids }

    sample_ids
    .collectFile(
        name: 'sample_ids.txt',
        storeDir: "${params.projectFolder}/documentation/",
        newLine: true, sort: true)

    //start pipeline

    if (params.qc) {
        //Run trimGalore
        trimGalore(reads, "${params.pairedEnd}")

        //Collect trimGalore stats (changed to within trimgalore folder)
        trimGalore.out.map{key, files -> files }
        .flatten()
        .filter(~/.*trim_stats\.txt/)
        .collectFile(
            name: 'trim_stats.txt',
            storeDir: "${params.outDir}/trimgalore/",
            newLine: false, sort: true)
        
        //Run strandedness
        strandedness = checkStrand(reads, "${params.pairedEnd}").map { keys, files -> keys }
        strandedness
        .collectFile(
            name: 'strandedness_all.txt',
            storeDir: "${params.outDir}/check_strandedness/",
            newLine: true, sort: true)
        strandedness.view()
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
            trimmed_reads = "${params.outDir}/trimgalore/**/*{R1,R2}_trimmed*.{fastq.gz,fq.gz}"
            if (file(trimmed_reads).isEmpty()) {
                star_input = reads
                println  "No trimmed reads found in path: ${trimmed_reads}, using ${params.readsPath}"

            } else {
                star_input = Channel
                .fromFilePairs("${trimmed_reads}", size: params.pairedEnd ? 2 : 1)   
            }

        }
        starAlign(star_input, "${params.pairedEnd}", indexLength(star_input))
        samtools(starAlign.out.bam)
    }

}

