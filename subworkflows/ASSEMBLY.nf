include { stringtie } from '../modules/stringtie'
include { mergeGTF; filterAnnotate; customAnotation } from '../modules/mergeTranscriptome'

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

workflow ASSEMBLY {
    take:
        strand
        bam

    main:
        //Run stringtie unless paths to precomputed individual sample GTF are provided
        if (!params.sample_gtf_list & params.assembly) {
             //set strand
            strand.map{ it -> getStrandtype(it) }.set{ strand }

            //locate chromosome exclusion list
            if(params.chr_exclusion_list) {
            //TODO now only able to access default assets, not other paths in the file system.
                chromosome_exclusion_list = file("${workflow.projectDir}/${params.chr_exclusion_list}")
                                            .readLines().join(",")
            } else {
                chromosome_exclusion_list = null
            }

            stringtie(strand, bam, chromosome_exclusion_list)

            //Wait untill all stringtie runs are completed
            stringtie.out.collect()
            .set { gtf_list }

            //Combine all paths of stringties results into one file
            gtf_list = gtf_list.flatten()
            .map { it -> it.toString() }
            .map { it -> it.replaceFirst("${workDir}/[^/]*/[^/]*/", "${params.outdir}/stringtie/") } //Replace the work dir path with the output dir
            .collectFile(
                name: 'gtflist.txt',
                storeDir: "${params.outdir}/stringtie/",
                newLine: true, sort: true )
        }

        if (params.merge) {

            if (params.sample_gtf_list) {
                gtf_list = Channel.fromPath("${params.sample_gtf_list}")
                gtf_list
                .splitText(){ it.trim() }
                .take(1)
                .ifEmpty { error "Could not find sample GTF files in: ${params.sample_gtf_list}" }
            }

            mergeGTF(gtf_list)

            gtf_novel = mergeGTF.out.flatten().filter(~/.*\.combined\.gtf/)
            gtf_tracking = mergeGTF.out.flatten().filter(~/.*\.tracking/)

            scripts_dir = Channel.fromPath("${workflow.projectDir}/bin/")

            filterAnnotate(params.reference_gtf,
                           params.refseq_gtf,
                           gtf_novel,
                           gtf_tracking,
                           1,
                           params.output_basename,
                           scripts_dir)

            merged_gtf = filterAnnotate.out.filter(~/.*_novel_filtered\.gtf/)

            //if (params.custom_annotation) {
            //   if ( params.custom_annotation ==~ /orfquant/ ) {
            //        customAnotation(merged_gtf)
            //   }
            //}
        }
}
