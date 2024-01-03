# Nextflow RNA-seq Pipeline

This Nextflow-based pipeline performs comprehensive RNA-seq analysis, including mapping and assembly (expression quantification and fusion calling coming up soon). 

## Overview

1. **Quality Control and Adapter Trimming**: Utilizes FastQC for assessing the quality of raw sequencing data and trims adapters and filters low-quality reads using Trimmomatic.
3. **Alignment**: Aligns processed reads to a reference genome using STAR.
4. **Transcriptome Assembly**: Assembles transcripts in individual samples with stringtie, merges all samples and checks against a reference with GFFcompare, and finally filters the merged transcriptome.

## Usage

1. **Configure Parameters**:
    
    Customize the parameters using the `params.config` file to suit your analysis.

    - **Inputs**: It accepts input reads (`reads_path`) and parameters specified in `params.config`. See an example [here](https://github.com/mars13/nf_rna_pipeline/blob/main/test/documentation/params.config). Channels are created for the input reads and processed to handle paired-end (single-end to be added as a future feature).

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

    **Note: to use the pipeline for merginging individual samples GTF without prior steps you are require to provide a file containing the paths to indivitual GTFs `sample_gtf_list="./stringtie/gtflist.txt"`**

    ```
    ./analysis/stringtie/sample_01/sample_01.gtf
    ./analysis/stringtie/sample_02/sample_02.gtf
    ```

2. **Run the Pipeline**:

    ```bash
    nextflow run mars13/nf_rna_pipeline -c params.config -profile [local/slurm]
    ```

    Additional nextflow run options can be provided. See [nextflow docs](https://www.nextflow.io/docs/latest/cli.html#run) for more information.


3. **Pull repository**:
    To keep up with the latest version pull the repository if it has already been downloaded.

    ```bash
    nextflow pull mars13/nf_rna_pipeline
    ```

## Output

The pipeline generates output files including quality reports, trimmed reads, alignment results and assembled transcripts. The output directory is specified by the parameter `outdir` and has the following structure:

```
{outdir}
├── check_strandedness/
├── customannotation/
├── gffcompare/
├── star/
├── stringtie/
└── trimgalore/
```

Additionally it produces:

- `r1_files.txt`, `r2_files.txt` and `sample_ids.txt` in `{project_folder}/documentation`
- Nextflow execution reports in `{project_folder}/log`

## Support and Contributions

- For issues or questions, [create an issue](https://github.com/mars13/nf_rna_pipeline/issues).
- Contributions to enhance this pipeline are welcomed through pull requests.
