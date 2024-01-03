# Nextflow RNA-seq Pipeline

This Nextflow-based pipeline performs comprehensive RNA-seq analysis, including mapping and assembly (expression quantification and fusion calling coming up soon). 

## Overview

1. **Quality Control and Adapter Trimming**: Utilizes FastQC for assessing the quality of raw sequencing data and trims adapters and filters low-quality reads using Trimmomatic.
3. **Alignment**: Aligns processed reads to a reference genome using STAR.
4. **Transcriptome Assembly**: Assembles transcripts in individual samples with stringtie, merges all samples and checks against a reference with GFFcompare, and finally filters the merged transcriptome.

## Usage

1. **Configure Parameters**:
    
    Customize the parameters using the `params.config` file to suit your analysis.

    - **Inputs**: Parameters are specified in `params.config`. See an example [here](https://github.com/mars13/nf_rna_pipeline/blob/main/test/documentation/params.config). Alternatively, they can be provided as command-line arguments if double-dashed (`--reads_path`). A detailed description of the required inputs can be found in the [Inputs](##inputs) section.

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

2. **Run the Pipeline**:

    ```bash
    nextflow run mars13/nf_rna_pipeline -c params.config -profile [local/slurm]
    ```

    Additional nextflow run options can be provided. See [nextflow docs](https://www.nextflow.io/docs/latest/cli.html#run) for more information.

    Add or overwrite params in the config files with double dashed arguments:
    
    ```bash
    nextflow run mars13/nf_rna_pipeline --align false --bam_files precomputed_aligment/**.Aligned.sortedByCoord.out.bam -c params.config -profile [local/slurm]
    ```



3. **Pull repository**:
    To keep up with the latest version pull the repository if it has already been downloaded.

    ```bash
    nextflow pull mars13/nf_rna_pipeline
    ```

## Inputs

### Required inputs

**General inputs**:
- `reads_path`: regular expression pointing to the reads location in `.fastq.gz` or `.fq.gz` extension. Example: "project_dir/data/*_{R1,R2}*.{fastq.gz,fq.gz}"
- `reference_gtf`: path to reference gtf.

**QC inputs**:
- `kallisto_index`: path to precomputed kallisto index. 

**Align inputs**:
- `star_index_basedir`: path to precomputed STAR index base directory. 

**Merge inputs**:
- `masked_fasta`: Masked fasta file for GFFcompare.
- `refseq_gtf`: Refseq GTF used in the filtering step.

### Optional inputs

- `strand_info`: When not running the qc step, the path to the file containing strand information for each sample can be provided. The file is expected to be plain text file (tab separated, no headers) with two columns: sample_id and strand information (RF/fr-firststrand, FR/fr-secondstrand, unstranded). If `strand_info` is not provided the software will check the default location: `"${outdir}/check_strandedness/strandedness_all.txt"`.
- `bam_files`: When not running the alignment step before assembly, the path to bam files can be provided as a regular expression. If `bam_files` is not provided the software will check the default location: `"${outdir}/star/**/*.Aligned.sortedByCoord.out.bam"`.
- `sample_gtf_list`: For merginging GTF files from individual samples without prior assembly steps, you are require to provide a file containing the paths to individual GTFs. Example:

```
./analysis/stringtie/sample_01/sample_01.gtf
./analysis/stringtie/sample_02/sample_02.gtf
```



## Outputs

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
