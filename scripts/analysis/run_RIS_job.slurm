#!/bin/bash
########## Slurm preamble
#SBATCH --job-name=RIS_timbre
#SBATCH --mail-user=jofritzi@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --account=microb
#SBATCH --partition=standard
#SBATCH --time=8:00:00
#SBATCH --output=RIS_timbre_%j.out
#SBATCH --error=RIS_timbre_%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=24
########## End of preamble!

module load matlab/R2025b
matlab -nodisplay -r "save_RIS_thresholds.m; exit"