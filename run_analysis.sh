#!/bin/bash 
#SBATCH -c 30  # Number of Cores per Task
#SBATCH --mem=200g  # Requested Memory
#SBATCH --partition=cpu
#SBATCH --account=pi_alc_umass_edu
#SBATCH -t 12:00:00  # Job time limit
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
module load r-rocker-ml-verse/4.4.0+apptainer
Rscript analyses/analysis_for_paper.R