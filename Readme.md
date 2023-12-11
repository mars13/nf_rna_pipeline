# Nextflow RNA-seq Pipeline

This Nextflow-based pipeline performs comprehensive RNA-seq analysis, including mapping and assembly (expression quantification and fusion calling coming up soon). 

## Overview

1. **Quality Control and Adapter Trimming**: Utilizes FastQC for assessing the quality of raw sequencing data and trims adapters and filters low-quality reads using Trimmomatic.
3. **Alignment**: Aligns processed reads to a reference genome using STAR.
4. **Transcriptome Assembly**: Assembles transcripts in individual samples with stringtie, merges all samples and checks against a reference with GFFcompare, and finally filters the merged transcriptome.

## Usage

1. **Clone the Repository**:

    ```bash
    git clone https://github.com/mars13/nf_rna_pipeline.git
    cd nf_rna_pipeline
    ```

2. **Configure Parameters**:
    
    Customize the parameters using the `params.config` file to suit your analysis.

    - **Inputs**: It accepts input reads (`readsPath`) and parameters specified in `params.config`.Channels are created for the input reads and processed to handle paired-end (single-end to be added as a future feature).
    - **Steps**: Adjust settings in `params.config` to toggle pipeline steps (`qc`, `align`, `assembly`, `merge`).
    
    ```
        params {
            //Toggle pipeline steps
            qc = false
            align = false
            assembly = false
            merge = true
        }
    ```

    - **Reference and other files**: Define paths to reference files, input data, and other settings in `params.config`.
    - **Containers**: Modify process containers in `nextflow.config` for tools like Trimmomatic, STAR, StringTie, etc., as per your environment.

    **Note: to use the pipeline for merginging indivitual samples GTF without prior steps you are require to provide a file containing the paths to indivitual GTFs `sampleGTFList="./stringtie/gtflist.txt"`**

    ```
    ./analysis/stringtie/sample_01/sample_01.gtf
    ./analysis/stringtie/sample_02/sample_02.gtf
    ```

3. **Run the Pipeline**:

    ```bash
    nextflow run main.nf -c params.config
    ```

    Additional nextflow run options can be provided. See [nextflow docs](https://www.nextflow.io/docs/latest/cli.html#run) for more information.
## Output

The pipeline generates output files including quality reports, trimmed reads, alignment results and assembled transcripts. The output directory is specified by the parameter `outdir` and has the following structure:

```
    check_strandedness/
    customannotation/
    gffcompare/
    star/
    stringtie/
    trimgalore/
```

## Support and Contributions

- For issues or questions, [create an issue](https://github.com/mars13/nf_rna_pipeline/issues).
- Contributions to enhance this pipeline are welcomed through pull requests.

## License

This pipeline is licensed under [MIT License](LICENSE).
