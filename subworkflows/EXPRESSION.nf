include { salmon_index; salmon_quasi; salmon_bam; salmon_tables; featurecounts} from '../modules/salmon'

/*
* Run salmon steps to obtain expression statistics
*/
workflow EXPRESSION {

    take:
    reads           // Input read file(s)
    bam             // Bam file created by star align
    assembled_gtf
    stringtie_transcriptome
    mode            // Salmon mode to run
    paired_end      // Bool, is data paired end or not
    transcriptome   // Path to the input transcriptome file
    reference_gtf   // Path to the reference gtf file
    output_basename // File prefix given to the salmon_tables results
    outdir          // Path to output directory

    main:
    if (params.stringtie_expression && assembled_gtf != null){
        input_gtf = assembled_gtf.first()
        input_transcriptome = stringtie_transcriptome.first()
    }
    else{
        input_gtf = reference_gtf
        input_transcriptome = transcriptome
    }
    
    if (mode =~ /sq/) {
        // Create the salmon index of the given transcriptome
        salmon_index(input_transcriptome)

        // Run salmon and write paths of quant.sf output files to a text file
        salmon_quasi(reads, paired_end, salmon_index.out, outdir)

        // Write the paths of the salmon_quasi output files to a text file
        quant_paths = salmon_quasi.out.quant
            .map { it -> it.toString() }
            .collectFile(
            name: 'quant_paths.txt',
            newLine: true, sort: true )

    } else if (mode =~ /sa/ ) {
        //salmon_bam(filtered_bam, transcriptome, outdir)
        println "This mode is not supported yet."
        exit 1
    }

    // Run the salmon_tables R script to obtain salmon statistics
    salmon_tables(quant_paths, input_gtf, output_basename, outdir)

    // Run featurecounts, TODO: figure out what subworkflow to put this 
    featurecounts(bam, input_gtf, outdir)

}
