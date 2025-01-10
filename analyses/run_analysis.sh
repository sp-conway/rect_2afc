#!/bin/bash 
#SBATCH -c 20  # Number of Cores per Task
#SBATCH --mem=30g  # Requested Memory
#SBATCH --partition=cpu
#SBATCH --account=pi_alc_umass_edu
#SBATCH -t 18:00:00  # Job time limit
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
module load r-rocker-ml-verse/4.4.2+apptainer
Rscript analysis_for_paper.R