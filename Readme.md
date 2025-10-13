# Nextflow RNA-seq Pipeline

This pipeline is designed for the analysis of RNA-seq data, covering various modules such as alignment, assembly, QC, fusion detection, and expression analysis. It is flexible and can currently handle paired-end sequencing data (with intention to expand to single-end for the modules that support it).

## Overview

1. **Quality Control and Adapter Trimming**: Utilizes FastP for assessing the quality of raw sequencing data and trims adapters and filters low-quality reads. Additionaly confirms strandedness for downstream transcriptome assembly.
2. **Alignment**: Aligns processed reads to a reference genome using STAR.
3. **Transcriptome Assembly**: Assembles transcripts of individual samples with stringtie, merges all samples and checks against a reference with GFFcompare, and finally filters the merged transcriptome.
4. **Fusions**: Detects gene fusions using Arriba, optionally takes vcf files from WGS to refine confidence thresholds.
5. **Expression Analysis**: Quantifies gene expression using Salmon.

## Usage

1. **Configure Parameters**:

    Customize the parameters using the `params.config` file to suit your analysis.

    - **Inputs**: Parameters are specified in `params.config`. See an example [here](https://github.com/mars13/nf_rna_pipeline/blob/dev/test/documentation/test_params.config). Alternatively, they can be provided as command-line arguments if double-dashed (`--samplesheet`). A detailed description of the required inputs can be found in the [Inputs](#inputs) section.

    - **Modules**: Adjust settings in `params.config` to toggle pipeline modules on and off (`qc`, `align`, `assembly`, `merge`, `fusion`, `expression`). All steps are set to `true` by default.
    
    ```
        params {
            //Toggle pipeline steps
            qc = false
            align = false
            assembly = false
            merge = true
        }
    ```

    - **Reference and other files**: Define paths to reference files, optional inputs, and other settings in `params.config`.

    - **Containers**: Modify process containers in `params.config` for tools like Trimmomatic, STAR, StringTie 2/3, etc., as per your environment. Note: This is expected to be provided with the pipeline in the future.

2. **Pull repository**:
    To keep up with the latest version pull the repository if it has already been downloaded.

    ```bash
    nextflow pull mars13/nf_rna_pipeline
    ```

3. **Run the Pipeline**:

    ```bash
    nextflow run mars13/nf_rna_pipeline -c params.config -profile [local/slurm]
    ```

    Additional nextflow run options can be provided. See [nextflow docs](https://www.nextflow.io/docs/latest/cli.html#run) for more information.


## Inputs

### Samplesheet
- The samplesheet is a critical input for the RNA pipeline. It is a structured CSV file that lists all samples and files to be processed. Each row represents a sample with detailed information required by the pipeline. Below is a breakdown of the expected columns and their constraints:

| **Column Name**  | **Description**                                                                     | **Required** | **Constraints**                                                                 |
|------------------|---------------------------------------------------------------------------------|--------------|---------------------------------------------------------------------------------|
| `subject_id`     | Unique identifier for the subject (e.g., patient).                              | Yes          | Must be a string without spaces.                                                |
| `sample_id`      | Unique identifier for the sample (e.g., sample barcode).                        | Yes          | Must be a string without spaces.                                                |
| `group_id`       | (Optional) Identifier for the sample group (e.g., cohort).                      | No           | Must be a string without spaces. If left empty, it must be omitted entirely.    |
| `sample_type`    | Type of sample: either `tumor` or `normal`.                                     | Yes          | Must be either `tumor` or `normal`.                                             |
| `sequence_type`  | Specifies the type of sequencing data: `rna` or `dna`.                          | Yes          | Must be either `rna` or `dna`.                                                  |
| `file_type`      | Format of the input files, e.g., `fastq`.   | Yes          | Must be one of the supported formats: `fastq`, `bam`, `cram`, `vcf`, `csv`, etc. *|
| `filename_1`     | Path to the first file (e.g., R1 FASTQ file for paired-end or single-end data).  | Yes          | Must be a valid file path with no spaces. File extension must match `file_type`. |
| `filename_2`     | (Optional) Path to the second file (e.g., R2 FASTQ file for paired-end data).    | No           | Must be a valid file path with no spaces. If not applicable, leave empty.       |

(*) Note: Not all suported formats can be used within the pipeline. I.e. bams/crams are still not implemented as input.

**Example samplesheet**:
```
subject_id,sample_id,group_id,sample_type,sequence_type,file_type,filename_1,filename_2
subject1,sampleA,cohort1,tumor,rna,fastq,/path/to/sampleA_R1.fastq.gz,/path/to/sampleA_R2.fastq.gz
subject2,sampleB,,normal,rna,bam,/path/to/sampleB.bam,
subject3,sampleC,cohort2,tumor,dna,vcf,/path/to/sampleC.vcf,
```

Notes:
- All file paths (filename_1 and filename_2) must not contain spaces and should have extensions that match the declared file_type.
- Subject and sample ids must contain charaters (only numeric values will not be read properly).
- For single-end data or non-FASTQ files (e.g., BAM, VCF), filename_2 can be omitted or left empty.

### Parameter specification
| **Parameter**               | **Description**                                                                                     | **Type**   | **Required** | **Used by Module**                | **Default**                          |
|-----------------------------|-----------------------------------------------------------------------------------------------------|------------|--------------|------------------------------------|--------------------------------------|
| **Path specifications**      |                                                                                                     |            |              |                                    |                                      |
| `input`                     | Path to the samplesheet (CSV) file containing sample information.                                    | string     | Yes          |                                    |                                      |
| `project_folder`            | Path to the folder where outputs and logs will be stored.                                            | string     | No           | All modules                        | `.`                                  |
| `outdir`                    | Name of the output folder.                                                                          | string     | No           | All modules                        | `${project_folder}/output`           |
| **Reference files**          |                                                                                                     |            |              |                                    |                                      |
| `reference_gtf`             | Path to the reference annotation GTF file.                                                          | string     | Yes          |                                    |                                      |
| `reference_transcriptome`    | Path to the reference transcriptome assembly.                                                       | string     | No           | `align`, `expression`, `assembly`  |                                      |
| `star_index_basedir`        | Path to the STAR index folder.                                                                       | string     | No           | `align`                            |                                      |
| `refseq_gtf`                | Path to the RefSeq GTF file.                                                                        | string     | No           | `align`, `expression`              |                                      |
| `masked_fasta`              | Path to the masked FASTA file.                                                                      | string     | No           | `align`                            |                                      |
| `species`                   | Species name.                                                                                       | string     | No           | All modules                        | `Homo_sapiens`                       |
| `genome_version`            | Genome version.                                                                                     | string     | No           | All modules                        | `GRCh38`                             |
| `annot_version`             | Annotation version.                                                                                 | string     | No           | `build_annotation`                 | `102`                                |
| `twobit`                    | Path to the 2bit file.                                                                              | string     | No           | `assembly`                         |                                      |
| `kallisto_index`            | Path to the Kallisto index file.                                                                    | string     | No           | `qc`                               |                                      |
| `arriba_reference`          | Path to the pre-built Arriba reference folder.                                                      | string     | No           | `fusions`                          |                                      |
| **Modules**                  |                                                                                                     |            |              |                                    |                                      |
| `qc`                        | Toggle for the quality control module.                                                              | boolean    | No           | `qc`                               | `true`                               |
| `align`                     | Toggle for the alignment module.                                                                    | boolean    | No           | `align`                            | `true`                               |
| `assembly`                  | Toggle for the assembly module.                                                                     | boolean    | No           | `assembly`                         | `true`                               |
| `expression`                | Toggle for the expression quantification module.                                                    | boolean    | No           | `expression`                       | `true`                               |
| `fusions`                   | Toggle for the fusion detection module.                                                             | boolean    | No           | `fusions`                          | `true`                               |
| `merge`                     | Toggle for merging individual GTF files into a single file.                                         | boolean    | No           | `merge`                             | `true`                               |
| `build_annotation`          | Toggle for building a custom annotation from merged GTFs.                                           | boolean    | No           | `build_annotation`                 | `true`                               |
| **Extended pipeline options**|                                                                                                     |            |              |                                    |                                      |
| `custom_annotation`         | Formats for custom annotation outputs (e.g., "orfquant").                                           | string     | No           | `build_annotation`                 | `"orfquant"`                         |
| `expression_mode`           | Expression quantification mode for tools like Salmon and FeatureCounts (e.g., "sqfc").              | string     | No           | `expression`                       | `"sqfc"`                             |
| `chr_exclusion_list`        | Path to the list of chromosomes to exclude during transcriptome assembly.                           | string     | No           | `assembly`                         | `assets/chr_exclusion_list.txt`      |
| `output_basename`           | Basename for merged GTF, custom annotation, and expression files.                                   | string     | No           | `merge`, `build_annotation`        | `transcriptome_full`                 |
| **Optional path definitions**|                                                                                                     |            |              |                                    |                                      |
| `resource_folder`           | Path to the annotation resource folder.                                                             | string     | No           | All modules                        |                                      |
| `container_folder`          | Path to the container location.                                                                     | string     | No           | All modules                        |                                      |
| `sample_gtf_list`           | Path to a file containing a list of individual sample GTFs for merging.                             | string     | No           | `merge`                            |                                      |
| `strand_info`               | Strandedness information, if available.                                                             | string     | No           | `assembly`                         |                                      |
| **Maximum specifications**   |                                                                                                     |            |              |                                    |                                      |
| `max_memory`                | Maximum memory allocation per process.                                                              | string     | No           | All modules                        |                                      |
| `max_cpus`                  | Maximum CPU allocation per process.                                                                 | integer    | No           | All modules                        |                                      |
| `max_time`                  | Maximum time allowed per process.                                                                   | string     | No           | All modules                        |                                      |                      |

### Optional inputs for modular execution

- `strand_info`: When not running the QC step, the path to the file containing strand information for each sample can be provided. The file is expected to be plain text file (tab separated, no headers) with two columns: sample_id and strand information (RF/fr-firststrand, FR/fr-secondstrand, unstranded). If `strand_info` is not provided the software will check the default location: `"${outdir}/check_strandedness/strandedness_all.txt"`.
- `bam_files`: When not running the ALIGN step before ASSEMBLY, the path to bam files can be provided as a regular expression. If `bam_files` is not provided the software will check the default location: `"${outdir}/star/**/*.Aligned.sortedByCoord.out.bam"`.
- `sample_gtf_list`: For merging GTF files from individual samples without prior assembly steps, you are require to provide a file containing the paths to individual GTFs. Example:

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
├── fastp/
├── gffcompare/
├── salmon/
├── star/
└── stringtie/
```

Additionally it produces Nextflow execution reports in `{project_folder}/log`

## Support and Contributions

- For issues or questions, [create an issue](https://github.com/mars13/nf_rna_pipeline/issues).
- Contributions to enhance this pipeline are welcomed through pull requests.
