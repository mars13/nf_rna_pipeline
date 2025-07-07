def printHeader () {
    def logMessage =  """
        .-------------------------------------------------------.
        |                _____                 _      _     _   |
        | _ _ ___ ___   |  |  |___ ___ ___ ___| |_   | |___| |_ |
        || | | .'|   |  |     | -_| -_|_ -|  _|   |  | | .'| . ||
        | \\_/|__,|_|_|  |__|__|___|___|___|___|_|_|  |_|__,|___||
        '-------------------------------------------------------'

        ${workflow.manifest.name} ${workflow.manifest.version}
        ==========================
        """
        log.info logMessage.stripIndent()
}


def check_duplicate_samples(String filePath) {
    // Empty set to store unique sample ids
    def unique_values = new HashSet()
    def has_duplicates = false

    // Read the CSV file 
    new File(filePath).eachLine { line ->
        // Skip empty lines
        if (line.trim()) {
            def columns = line.split(',')
            // Obtain values from first columns
            def first_column_value = columns[1].trim()

            // Check if the value is unique
            if (unique_values.contains(first_column_value)) {
                // If duplicate found end loop and return true
                has_duplicates = true
                return
            } else {
                // Add value to set
                unique_values.add(first_column_value)
            }
        }
    }

    // Return the result: true if duplicates found, false otherwise
    return has_duplicates
}


def check_files(name, path, type) {
    if(!path) {
        error "When running merge without assembly you must provide `${name}`."
    } else {
        if (type == "dir") {
            file_to_check = file(path, type: "dir")
        } else {
            file_to_check = file(path, type: "file")
        }
        if (file_to_check != null) {
            // Check if it's a list of files
            if (file_to_check instanceof List) {
                // Check each file in the list
                file_to_check.flatten().each { file ->
                    if (!file.exists()) {
                        error "--${name}: ${type} doesn't exist, check path ${file}"
                    }
                }
            } else {
                // Check the single file
                if (!file_to_check.exists()) {
                    error "--${name}: ${type} doesn't exist, check path ${path}"
                }
            }
        } else {
            error "--${name}: No files found at ${path}"
        }
    }
}

def checkInputFiles() {
    // Check inputs
    default_strand =  "${params.outdir}/check_strandedness/strandedness_all.txt"
    if (!params.qc && ( params.align || params.assembly )) {
        if (!params.strand_info && file(default_strand, type : "file").exists()) {
            log.info "strand info   : ${default_strand}".stripIndent()
        } else if (params.strand_info && file(params.strand_info, type : "file").exists()) {
            log.info "strand info   : ${params.strand_info}".stripIndent()
        } else {
            error  """
                No strand information found in ${default_strand}
                If `qc = false`, define your strandedness info file with `strand_info` in `params.config`
                """.stripIndent()
        }
    }

    // Locate bams
    default_bams = "${params.outdir}/star/**/*.Aligned.sortedByCoord.out.bam"
    bam_avail = true
    if (!params.align) {
        if (!params.bam_files && !file(default_bams).isEmpty()) {
            log.info "bam files     : ${default_bams}".stripIndent()
        } else if (params.bam_files && !file(params.bam_files).isEmpty()) {
            log.info "bam files     : ${params.bam_files}".stripIndent()
        } else {
            bam_avail = null
        }
    }

    // Check reference_gtf
    check_files("reference_gtf", params.reference_gtf, "file")

    // Check kallisto index
    if (params.qc){
        check_files("kallisto_index", params.kallisto_index, "file")
    }

    // Check star_index_basedir
    if (params.align){
        check_files("star_index_basedir", params.star_index_basedir, "dir")
    }

    // Check bams and references for assembly steps
    if (params.assembly){
        if (!bam_avail) {
            error  """
                No bam found in ${default_bams}
                If `align = false`, define your bam location with `bam_files=[path_to_file]` in `params.config`
                """.stripIndent()
        }

        check_files("refseq_gtf", "${params.refseq_gtf}*", "file")
        check_files("masked_fasta", "${params.masked_fasta}", "file")
    }

    // Check sample gtf list for merge without assembly
    if (!params.assembly && params.merge) {
        if(!params.sample_gtf_list) {
            error "When running merge without assembly you must provide `sample_gtf_list`.".stripIndent()
        } else {
            sample_gtf_list = file(params.sample_gtf_list, type: "file")
            if (!sample_gtf_list.exists()) {
                error  "--sample_gtf_list: file doesn't exist, check path ${params.sample_gtf_list}".stripIndent()
            }
        }
    }

    // Check references for build_annotaiton
    if (params.build_annotation) {
        check_files("twobit", "${params.twobit}*", "file")
    }

    // Check fusion parameters
    if (params.fusions){
        check_files("arriba_reference", params.arriba_reference, "dir")
        check_files("arriba_reference STAR_index_*", "${params.arriba_reference}STAR_index*", "dir")
        check_files("arriba_reference *.gtf", "${params.arriba_reference}*.gtf", "file")
        check_files("arriba_reference *.fa", "${params.arriba_reference}*.fa", "file")

        if (params.wgs_sv) {
            check_files("wgs_sv", "${params.wgs_sv}", "dir")
        }
    }

    // Check expression parameters
    if (params.expression) {
        assert params.expression_mode in ["sq", "sa", "sqfc", "safc"], "`expression_mode` must be one of the following: sq, sa, sqfc, safc"

        if (!bam_avail & (params.expression_mode != "sq")) {
            log.info "expression mode  : ${params.expression_mode} --> sq [forced to quasi-mapping, no bam avail]"
        }
    }

    // Check for duplicate sample ids
    if (check_duplicate_samples(params.input)) {
        error "\nERROR: Duplicate sample_ids found in the input samplesheet\n"
    } 

    log.info "\n==========================\n"
}

/*
 * Checks whether all samples in a samplesheet are paired-end or single-end.
 * Throws an error if there's a mix.
 *
 * @param samplesheet_path path to CSV samplesheet
 * @return true if all samples are paired-end, false if all are single-end
 */
def is_paired_end(String filePath) {
    def uniqueFlags = new HashSet()
    def headerParsed = false
    int idx_filename_1 = -1
    int idx_filename_2 = -1

    new File(filePath).eachLine { line ->
        if (!line.trim()) return  // skip empty lines

        def cols = line.split(",")*.trim()

        if (!headerParsed) {
            idx_filename_1 = cols.indexOf("filename_1")
            idx_filename_2 = cols.indexOf("filename_2")

            if (idx_filename_1 == -1 || idx_filename_2 == -1) {
                throw new IllegalArgumentException("Samplesheet missing required columns 'filename_1' or 'filename_2'")
            }
            headerParsed = true
        } else {
            def f1 = (idx_filename_1 < cols.size()) ? cols[idx_filename_1] : null
            def f2 = (idx_filename_2 < cols.size()) ? cols[idx_filename_2] : null

            if (!f1) {
                throw new IllegalArgumentException("Missing filename_1 in row: ${line}")
            }

            uniqueFlags.add(f2 ? true : false)

            if (uniqueFlags.size() > 1) {
                throw new IllegalArgumentException("Mixed paired-end and single-end data detected in samplesheet! Please only use one type of data")
            }
        }
    }

    if (!headerParsed) {
        throw new IllegalArgumentException("Samplesheet is empty or missing header.")
    }
    if (uniqueFlags.isEmpty()) {
        throw new IllegalArgumentException("Samplesheet contains no data rows.")
    }

    return uniqueFlags.iterator().next()
}
