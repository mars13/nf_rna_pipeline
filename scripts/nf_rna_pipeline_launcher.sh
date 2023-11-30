#!/bin/bash
#SBATCH --job-name=nfRNApipe
#SBATCH --output=/hpc/pmc_vanheesch/projects/Marina/test/nf_rna_pipeline/log/nfRNApipe.out
#SBATCH --error=/hpc/pmc_vanheesch/projects/Marina/test/nf_rna_pipeline/log/nfRNApipe.err
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --nodes=1
#SBATCH --open-mode=append



module load java
module load nextflow/23.04.4

nextflow run scripts/main.nf -c documentation/params.config -profile slurm