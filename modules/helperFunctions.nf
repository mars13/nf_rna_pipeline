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
        project folder   : ${params.projectFolder}
        input from       : ${params.readsPath}
        output to        : ${params.outDir}
        reference files  : ${params.resourceDir}/GENOMES/${params.species}.${params.genome_version}/${params.annot_version}/
        sample gtf list  : ${params.sampleGTFList}
        --
        qc         : ${params.qc}
        align      : ${params.align}
        assembly   : ${params.assembly}
        merge      : ${params.merge}
        ==========================
        """
        .stripIndent()
}

