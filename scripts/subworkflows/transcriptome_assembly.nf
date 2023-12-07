include { stringtie } from '../modules/sampleGTF'
include { mergeGTF } from '../modules/mergeTranscriptome'

def getStrandtype(strand) {
        if (strand ==~ /.*RF\/fr-firststrand/) {
        strand_type = "rf"
        } else if (strand ==~ /.*FR\/fr-secondstrand/) {
            strand_type = "fr"
        } else {
            println "Data must be stranded"
            exit 1
        }
}

workflow transcriptome_assembly {
    take: 
        strand
        bam

    main:
        //set strand
        strand.map{ it -> getStrandtype(it) }.view().set{ strand }
        
        //locate chromosome exclusion list
        if(params.chrExclusionList) {
            chromosome_exclusion_list = file(params.chrExclusionList)
                                        .readLines().join(",")
        } else {
            chromosome_exclusion_list = null
        }
        
        //Run stringtie unless paths to precomputed individual sample GTF are provided
        if (!params.sampleGTFList) {
            stringtie(strand, bam, chromosome_exclusion_list)
            
            stringtie.out
            .map { it -> it.toString() }
            .collectFile(
                name: 'gtflist.txt',
                storeDir: "${params.outDir}/stringtie/",
                newLine: true, sort: true)
            .set { gtf_list }
        } else {
            gtf_list = Channel.fromPath("${params.sampleGTFList}")
        }
        gtf_list.view()
        if (merge) {
            mergeGTF(gtf_list)
            mergeGTF.out.view()
        }



}