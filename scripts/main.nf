nextflow.enable.dsl = 2

include { trimGalore; checkStrand; indexLength; starAlign; samtools} from './modules/rnaseq_alignment'



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
    apptainer exec -B "/hpc:/hpc" --env "LC_ALL=C.UTF-8" ${container_dir}/trimgalore-0.6.6.sif cutadapt --version >> versions.txt
    apptainer exec -B "/hpc:/hpc" --env "LC_ALL=C.UTF-8" ${container_dir}/trimgalore-0.6.6.sif fastqc --version >> versions.txt
    apptainer exec -B "/hpc:/hpc" --env "LC_ALL=C.UTF-8" ${container_dir}/trimgalore-0.6.6.sif trim_galore --version >> versions.txt

    """
    
}


workflow {
    
    //check params

    reads =  Channel
            .fromFilePairs( params.readsPath, size: params.pairedEnd ? 2 : 1 )
            .ifEmpty { error "Could not find any reads matching pattern: ${params.readsPath}" }
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
        
        //Collect trimGalore outputs 
        trimGalore.out
        .map({key, file ->
            tuple( key, 
                file.findAll({ it =~ /.*(?:R1|R2)_trimmed.*\.fastq\.gz$/ })
                )
            }).set{ star_input }

        //Run strandedness
        strandedness = checkStrand(reads, "${params.pairedEnd}")
        strandedness.map { keys, files -> keys }
        .collectFile(
            name: 'strandedness_all.txt',
            storeDir: "${params.outDir}/check_strandedness/",
            newLine: true, sort: true)
    }   
    
    if (params.align) {
        if (params.qc){

        } else {


        
        }
        starAlign(star_input, "${params.pairedEnd}", indexLength(star_input))
        samtools(starAlign.out.bam)
    }

}

