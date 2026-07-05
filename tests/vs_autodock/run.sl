#!/bin/bash
#SBATCH --job-name=vs_autodock
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err

# Purge loaded modules
module purge

# Activate conda environment, see environment.yml
module load Miniconda3
source activate /path/to/env/biobb_vs   # e.g. /shared/work/BiobbWorkflows/envs/biobb_vs

# Input files
INPUT_PATH=../../data
STRUCTURE_PATH=$INPUT_PATH/receptor/receptor.pdb
LIGAND_LIB=$INPUT_PATH/ligands/zinc_200_425_001_reduced.sdf

# Launch workflow 
vs_autodock --ligand_lib $LIGAND_LIB \
            --structure_path $STRUCTURE_PATH \
            --keep_poses \
            --pocket_selection "resid 37 or resid 49 or resid 112" \
            --box_offset 5 \
            --cpus 4 \
            --exhaustiveness 8