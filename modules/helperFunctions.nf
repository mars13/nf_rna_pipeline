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
        reference dir    : ${params.resource_folder}/GENOMES/${params.species}.${params.genome_version}/${params.annot_version}/
        reference gtf    : ${params.reference_gtf}
        sample gtf list  : ${params.sample_gtf_list}
        --
        qc          : ${params.qc}
        align       : ${params.align}
        assembly    : ${params.assembly}
        merge       : ${params.merge}
        build annot : ${params.build_annotation}
        expression  : ${params.expression}
        fusions     : ${params.fusions}
        ==========================
        """
        .stripIndent()
}


def check_files(name, path, type) {
    if(!path) {
        error "When running merge without assembly you must provide `${name}`."
    } else {
            file_to_check = file(path, type: "${type}")
            if (!file_to_check.exists()) {
                error  "--${name}: ${type} doesn't exist, check path ${path}"
            }
        }
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
    check_files("reference_gtf", params.reference_gtf, "file")

    //Check kallisto index
    if (params.qc){
        check_files("kallisto_index", params.kallisto_index, "file")
    }
    
    //Check star_index_basedir
    if (params.align){
        check_files("star_index_basedir", params.star_index_basedir, "dir")
    }

    //Check refseq_gtf and masked_fasta for assembly steps
    if (params.assembly){
        check_files("refseq_gtf", "${params.refseq_gtf}*", "file")
        check_files("masked_fasta", ${params.masked_fasta}, "file")
    }

    if (params.build_annotation) {
        check_files("twobit", "${params.twobit}*", "file")
    }

}

