workflow RNASEQ_ALIGNMENT {
    take: 
        readsPath
        pairedEnd
        qc
        align
        outDir
    emit:
        strand = strandedness.out
        bam = samtools.out

    main:
        
    
}