#!/bin/bash
#SBATCH --account=<user_account>
#SBATCH --time=<12:00:00>
#SBATCH --cpus-per-task=<1>
#SBATCH --mem=<16G>
#SBATCH --job-name=<laser_pca>
#SBATCH --output=launch_da.out
#SBATCH --error=launch_da.err


module load nextflow


nextflow run main.nf -c nextflow.config -resume


echo "DONE"
