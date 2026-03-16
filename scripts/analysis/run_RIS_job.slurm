#!/bin/bash
########## Slurm preamble
#SBATCH --job-name=RIS_timbre
#SBATCH --mail-user=jofritzi@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --account=microb0
#SBATCH --time=08:00:00
#SBATCH --output=RIS_timbre_%j.out
#SBATCH --error=RIS_timbre_%j.err
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=standard
########## End of preamble!

module load matlab/R2025b
matlab -batch "save_RIS_thresholds"
