#!/bin/bash
#SBATCH --job-name=cavity_analysis
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err

# Purge loaded modules
module purge 

# Activate conda environment, see environment.yml
module load Miniconda3
source activate /path/to/env/biobb_vs   # e.g. /shared/work/BiobbWorkflows/envs/biobb_vs

# Launch workflow
cavity_analysis --structures_path ../../data/receptor/receptor.pdb \
                --filtering_selection "resid 31 or resid 21"