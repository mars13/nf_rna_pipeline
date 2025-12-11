process MULTIQC {
    label "multiqc"
    publishDir "${outdir}/multiqc/", mode: 'copy'

    input:
        path multiqc_files, stageAs: "?/*"
        val multiqc_config
        val outdir

    output:
        path "*multiqc_report.html", emit: multiqc_report

    script:
        """
        multiqc . -c ${multiqc_config}
        """
}