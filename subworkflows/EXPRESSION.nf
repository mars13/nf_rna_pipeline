include { salmon_index; salmon_quasi; salmon_bam; salmon_tables} from '../modules/salmon'

workflow EXPRESSION {

    take:
        reads
        bam
        mode
        paired_end
        transcriptome
        outdir

    main:
        if (mode =~ /sq/) {
            salmon_index(transcriptome)
            salmon_quasi(reads, paired_end, salmon_index.out, outdir)
        } else if  (mode =~ /sa/ ) {

            //salmon_bam(bam, transcriptome, outdir)
            //Does work but needs a new transcriptome with corrects IDs
            //or needs to filter the bam file with the chr_exclusion_list


            println "Container for salmon alignment mode not available."
            exit 1
        }
        //TODO implement salmon_tables
        salmon_tables("$outdir/salmon", transcriptome, params.output_basename, salmon_quasi.out)

}

//featurecounts