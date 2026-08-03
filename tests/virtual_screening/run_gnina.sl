#!/bin/bash
#SBATCH --job-name=virtual_screening_gnina
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err
# Add a GPU request here (e.g. #SBATCH --gres=gpu:1) - gnina's CNN rescoring is much
# faster on a GPU. Without one, pass --gnina_no_gpu and expect minutes per ligand.

# Purge loaded modules
module purge

# Activate conda environment, see environment.yml
module load Miniconda3
source activate /path/to/env/biobb_vs   # e.g. /shared/work/BiobbWorkflows/envs/biobb_vs

# gnina is not conda installable and the conda-forge biobb_vs ships no gnina module:
# this workflow needs a gnina release binary and a biobb_vs built from source.
# https://github.com/gnina/gnina/releases
GNINA_BIN=/path/to/gnina

# Input files
INPUT_PATH=../../data
STRUCTURE_PATH=$INPUT_PATH/receptor/receptor.pdb
LIGAND_LIB=$INPUT_PATH/ligands/imatinib_analogs_prepared.sdf

# Launch workflow
# --gnina_cnn fast is a single CNN model, ~4x faster than the default 3 model ensemble.
# Drop it for the default ensemble, or use --gnina_cnn_scoring none to skip the CNN entirely.
virtual_screening --ligand_lib $LIGAND_LIB \
                  --structure_path $STRUCTURE_PATH \
                  --keep_poses \
                  --pocket_selection "resid 37 or resid 49 or resid 112" \
                  --box_offset 5 \
                  --docking_engine gnina \
                  --gnina_bin $GNINA_BIN \
                  --gnina_cnn fast \
                  --gnina_seed 42 \
                  --cpus 4 \
                  --exhaustiveness 8 \
                  --output gnina_out
