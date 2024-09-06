include { salmon_index; salmon_quasi; salmon_bam; salmon_tables} from '../modules/salmon'

workflow EXPRESSION {

    take:
        reads
        bam  // bam file created by star align
        mode  // salmon mode to run
        paired_end
        transcriptome
        reference_gtf
        outdir
        output_basename // file prefix given to the salmon_tables results

    main:
        if (mode =~ /sq/) {
            salmon_index(transcriptome)
            //Run salmon and write paths of quant.sf output files to a text file
            quant_paths = salmon_quasi(reads, paired_end, salmon_index.out, outdir)
                .quant
                .map { it -> it.toString() }
                .collectFile(
                name: 'quant_paths.txt',
                newLine: true, sort: true )
            quant_paths.view()

        } else if  (mode =~ /sa/ ) {

            //salmon_bam(bam, transcriptome, outdir)
            //Does work but needs to filter the bam file with the chr_exclusion_list

            println "Container for salmon alignment mode not available."
            exit 1
        }
        salmon_tables(quant_paths, reference_gtf, output_basename)

}

//featurecounts