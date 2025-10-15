#!/bin/bash
#SBATCH --account=ctb-hussinju
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=laser_pca
#SBATCH --output=parallel_launch_nextflow_PCA.out
#SBATCH --error=parallel_launch_nextflow_PCA.err


module load nextflow


nextflow run parallel_PCA_on_HGDP1KG.nf -c parallel_nextflow.config -resume



echo "DONE"
