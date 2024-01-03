def printStartLogs () {
    log.info """ 
        .-------------------------------------------------------.
        |                _____                 _      _     _   |
        | _ _ ___ ___   |  |  |___ ___ ___ ___| |_   | |___| |_ |
        || | | .'|   |  |     | -_| -_|_ -|  _|   |  | | .'| . ||
        | \\_/|__,|_|_|  |__|__|___|___|___|___|_|_|  |_|__,|___||
        '-------------------------------------------------------'

        ${workflow.manifest.name} ${workflow.manifest.version}
        ==========================
        run as       : ${workflow.commandLine}
        started at   : ${workflow.start}
        config files : ${workflow.configFiles}
        --
        project folder   : ${params.project_folder}
        input from       : ${params.reads_path}
        output to        : ${params.outdir}
        reference files  : ${params.resource_folder}/GENOMES/${params.species}.${params.genome_version}/${params.annot_version}/
        sample gtf list  : ${params.sample_gtf_list}
        --
        qc         : ${params.qc}
        align      : ${params.align}
        assembly   : ${params.assembly}
        merge      : ${params.merge}
        ==========================
        """
        .stripIndent()
}


def check_params() {
    //Check inputs
    if (file(params.reads_path).isEmpty()) {
        error  "--reads_path: File doesn't exist, check path ${params.reads_path}"
    }

    default_strand =  "${params.outdir}/check_strandedness/strandedness_all.txt"    
    if (!params.qc && ( params.align || params.assembly )) {            
        if (!params.strand_info && file(default_strand, type : "file").exists()) {
            log.info "strand info   : ${default_strand}"
        } else if (params.strand_info && file(params.strand_info, type : "file").exists()) {
            log.info "strand info   : ${params.strand_info}"
        } else {
            error  """
                No strand information found in ${default_strand} and no alternative path provided. 
                When running with `qc = false`, you can define your strandedness info file as --strand_info option from command line  `strand_info=[path_to_file]` in `params.config`
                The file is expected to be plain text file (tab separated, no headers) with two columns: sample_id and strand information (RF/fr-firststrand, FR/fr-secondstrand, unstranded).
                """
        }
    }

    default_bams = "${params.outdir}/star/**/*.Aligned.sortedByCoord.out.bam"
    if (!params.align && params.assembly) {            
        if (!params.bam_files && !file(default_bams).isEmpty()) {
            log.info "bam files     : ${default_bams}"
        } else if (params.bam_files && !file(params.bam_files).isEmpty()) {
            log.info "bam files     : ${params.bam_files}"
        } else {
            error  """
                No bam found in ${default_bams}
                When running with `align = false`, you can define your bam location as --bam_files option from command line or `bam_files=[path_to_file]` in `params.config`.
                Please, remember to index your bams.
                """
        }
    }

    if (!params.assembly && params.merge) {            
        if(!params.sample_gtf_list) {
            error "When running merge without assembly you must provide `sample_gtf_list`."
        } else {
            sample_gtf_list = file(params.sample_gtf_list, type: "file")
            if (!sample_gtf_list.exists()) {
                error  "--sample_gtf_list: file doesn't exist, check path ${params.sample_gtf_list}"
            }
        }
    }

    //Check reference_gtf
    if (!params.reference_gtf) {
            error  "You must provide `reference_gtf`."
    } else {
        reference_gtf = file(params.reference_gtf, type: "file")
        if (!reference_gtf.exists()) {
            error  "--reference_gtf: file doesn't exist, check path ${params.reference_gtf}"
        }
    }
    
    //Check kallisto index
    if (params.qc){
        if (!params.kallisto_index) {
            error  "You must provide `kallisto_index`."
        } else {
            kallisto_index = file(params.kallisto_index, type: "file")
            if (!kallisto_index.exists()) {
                error  "--kallisto_index: File doesn't exist, check path ${params.kallisto_index}"
            }
        }   
    }
    
    //Check star_index_basedir
    if (params.align){
        if (!params.star_index_basedir) {
            error  "You must provide `star_index_basedir`."
        } else {
            star_index_basedir = file(params.star_index_basedir, type: "dir")
            if (!star_index_basedir.exists()) {
                error  "--star_index_basedir: Directory doesn't exist, check path ${params.star_index_basedir}"
            }
        }   
    }

    //Check refseq_gtf and masked_fasta for assembly steps
    if (params.assembly){
        if (!params.refseq_gtf) {
            error  "You must provide `refseq_gtf`."
        } else {
            refseq_gtf = file("${params.refseq_gtf}*")
            if (refseq_gtf.isEmpty()) {
                error  "--refseq_gtf: File doesn't exist, check path ${params.refseq_gtf}*"
            }
        } 

        if (!params.masked_fasta) {
            error  "You must provide `masked_fasta`."
        } else {
            masked_fasta = file(params.masked_fasta, type: "file")
            if (!masked_fasta.exists()) {
                error  "--masked_fasta: File doesn't exist, check path ${params.masked_fasta}"
            }
        }     
    }

    //TODO add check for twobit and reference_genome for build_custom_annotation

}

