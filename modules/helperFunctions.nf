def printStartLogs () {
    log.info """                                                      
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

