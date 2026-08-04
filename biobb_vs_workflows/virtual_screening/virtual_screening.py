#!/shared/work/BiobbWorkflows/envs/biobb_sp_virtual_screening/bin/python

# Importing all the needed libraries
import os
import re
import sys
import json
import time
import shutil
import argparse
from typing import Dict, List, Pattern, Tuple, Union, Optional

# Load pdb parser from biopython
from Bio.PDB import PDBParser

import MDAnalysis as mda

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_vs_workflows.common import to_yaml
from biobb_structure_utils.utils.str_check_add_hydrogens import str_check_add_hydrogens
from biobb_structure_utils.utils.extract_residues import extract_residues
from biobb_structure_utils.utils.extract_molecule import extract_molecule
from openbabel import pybel
from biobb_vs.utils.box import box
from biobb_vs.gnina.gnina_run import gnina_run
from biobb_vs.fpocket.fpocket_select import fpocket_select
from biobb_vs.vina.autodock_vina_run import autodock_vina_run
from biobb_chemistry.babelm.babel_convert import babel_convert

# Docking step and biobb path keys of each engine
DOCKING_ENGINES = {
    'vina': {'step': 'step5_autodock_vina_run',
             'receptor': 'input_receptor_pdbqt_path',
             'outputs': ('output_pdbqt_path', 'output_log_path')},
    'gnina': {'step': 'step5b_gnina_run',
              'receptor': 'input_receptor_path',
              'outputs': ('output_sdf_path', 'output_summary_path')}}

# gnina score to rank by, mapped to whether the highest value is the best one
GNINA_SORT_ORDERS = {'CNNaffinity': True, 'CNNscore': True, 'minimizedAffinity': False}


def find_matching_str(pattern: Union[str, Pattern[str]], 
                      filepath: str
    ) -> Optional[str]:
    '''
    Finds the first match of a regular expression pattern in a file.
    Returns the matching string or None if there is no match.

    Parameters
    ----------
    pattern : Union[str, Pattern[str]]
        Regular expression pattern to search in file lines.
    filepath : str
        File path to search in.

    Returns
    -------
    Optional[str]
        String matching the pattern or None if there is no match.
    '''
    try:
        with open(filepath, 'r') as file:
            for line in file:
                match = re.search(pattern, line)
                if match:
                    return match.group(1)
    except FileNotFoundError:
        print(f"File not found: {filepath}")
    except IOError as e:
        print(f"An I/O error occurred: {e}")
    
    return None

def get_affinity(pdbqt_path: str) -> Optional[float]:
    '''
    Find best binding affinity from pdbqt file from AutoDock Vina. 

    The first line with the following structure contains the best affinity:

    REMARK VINA RESULT:    -4.259      0.000      0.000

    Where the first number is the binding affinity.
    
    Inputs
    ------

        pdbqt_path: path to pdbqt file

    Outputs
    -------

        affinity: best binding affinity
    '''

    affinity_pattern=r'REMARK VINA RESULT:\s+(-?\d+\.\d+)'

    affinity = find_matching_str(affinity_pattern, pdbqt_path)

    if affinity:
        return float(affinity)
    else:
        return None

def get_gnina_score(summary_path: str, rank_by: str) -> Optional[Dict]:
    '''
    Find the best scoring pose in the summary json file written by gnina.

    The file holds one entry per pose:

    "pose1": {"ligand_name": "CHEMBL1642", "ligand_index": 1, "pose": 1,
              "minimizedAffinity": -8.751, "CNNscore": 0.83, "CNNaffinity": 6.42}

    Best is the highest CNNaffinity/CNNscore or the lowest minimizedAffinity. Which scores
    are present depends on the cnn_scoring property. The poses are not sorted by rank_by,
    gnina orders them by its own pose_sort_order, so all of them are looked at.

    Inputs
    ------

        summary_path : path to gnina summary json file
        rank_by      : score to rank the poses by

    Outputs
    -------

        pose: entry of the best scoring pose, or None if there is no usable score
    '''

    try:
        with open(summary_path) as file:
            poses = json.load(file)
    except (OSError, ValueError):
        return None

    # a pose may be missing the score, gnina only writes some of them for some cnn_scoring values
    scored = [pose for pose in poses.values() if isinstance(pose.get(rank_by), (int, float))]

    if not scored:
        return None

    best = max if GNINA_SORT_ORDERS[rank_by] else min

    return best(scored, key = lambda pose: pose[rank_by])

def ligand_step_path(output_path: str, ligand_name: str, step_output_path: str) -> str:
    '''
    Path to a step output file inside a ligand subfolder, given the same path in the global paths dict

    Inputs
    ------

        output_path      : path to output directory
        ligand_name      : name of the ligand
        step_output_path : path to the step output file in the global paths dict

    Outputs
    -------

        path: path to the step output file inside the ligand subfolder
    '''

    step_name = os.path.basename(os.path.dirname(step_output_path))

    return os.path.join(output_path, ligand_name, step_name, os.path.basename(step_output_path))

def save_ranking(ranking: List[Tuple],
                 num_top_ligands: Union[int, None],
                 ranking_path: str,
                 score_columns: Tuple = ("Affinity",)
    ) -> List[str]:
    '''
    Create file with ranking of ligands according to affinity

    Inputs
    ------

        ranking          : list with tuples -> (score, index, id, scores), already ordered, best first
        num_top_ligands  : number of top ligands to save in the ranking, if None all ligands are saved
        ranking_path     : path to ranking file
        score_columns    : scores reported for each ligand, in the order they are written

    Output
    ------

        top_ligand_indices : list with indices of top ligands
    '''

    # Find number of top ligands to save
    if num_top_ligands is None:
        ranking_length = len(ranking)
    else:
        ranking_length = min(num_top_ligands, len(ranking)) 

    # Extract top ligands (ranking is already sorted)
    top_ligands = ranking[:ranking_length]
    top_ligand_indices = []
    
    # Create summary file with top ligands
    with open(ranking_path, 'w') as file:

        # Write header
        file.write("Rank," + ",".join(score_columns) + ",Index,Identifier \n")

        # For each ligand
        for rank, affinity_tuple in enumerate(top_ligands):

            _, ligand_index, ligand_id, scores = affinity_tuple

            # Write line
            file.write(f"{rank+1}," + ",".join(str(score) for score in scores) + f",{ligand_index},{ligand_id}\n")

            # Add ligand name to list
            top_ligand_indices.append(ligand_index)

    return top_ligand_indices

def validate_step(*output_paths: str) -> bool:
    '''
    Check all output files exist and are not empty
    
    Inputs
    ------

        *output_paths (str): variable number of paths to output file/s

    Output
    ------

        validation_result (bool): result of validation
    '''

    # Initialize value 
    validation_result = True

    # Check existence of files
    for path in output_paths:
        validation_result = validation_result and os.path.exists(path)

    # Check files are not empty if they exist
    if (validation_result):

        for path in output_paths:
            file_not_empty = os.stat(path).st_size > 0
            validation_result = validation_result and file_not_empty

    return validation_result

def read_ligand_lib(ligand_lib_path: str) -> Tuple[List[str], List[str]]:
    '''
    Read all ligand identifiers from ligand library file. The expected format is one of the following:

    - Format 1:

            ligand1_id \n
            ligand2_id

    - Format 2:

            ligand1_id  ligand1_name \n
            ligand2_id  ligand2_name

    Where ligand_id is a unique identifier (e.g. SMILES) and name_ligand is a string with the ligand name
    
    Inputs
    ------

        ligand_lib_path (str): path to ligand library file

    Output
    ------

        ligand_ids     : list of ligand identifiers
        ligand_names   : list with ligand names
    '''

    ligand_ids = []
    ligand_names = []

    # Open file
    with open(ligand_lib_path) as file:

        # Read all lines
        ligand_lib = file.readlines()

        # Process every line
        for index, line in enumerate(ligand_lib):

            line = line.split()

            # Append ligand ID to list
            ligand_ids.append(line[0])

            # If there is no name, use index as name
            if len(line)>1:
                ligand_names.append(line[1])
            else:
                ligand_names.append(str(index))

        # If there are no ligands, raise an error
        if len(ligand_ids) == 0:
            raise ValueError(f"No ligands found in ligand library file {ligand_lib_path}")
        
    return ligand_ids, ligand_names

def write_smiles(smiles: str, smiles_path: str):
    '''
    Writes a SMILES code into a file. If the file exists, it will be overwritten.

    Inputs
    ------
        smiles         :  SMILES code
        smiles_path    :  smiles file path
    '''

    # Save SMILES in tmp file inside step_path, overwrite if exists
    smiles_tmp_file = open(smiles_path, 'w')
    smiles_tmp_file.write(smiles)
    smiles_tmp_file.close()

def get_ranking(
    ligand_ids: List, 
    ligand_names: List, 
    autodock_vina_paths: Dict, 
    output_path: str
    ) -> List[Tuple]:
    """
    Takes the name of each ligand and finds the best affinity given by AutoDock Vina from the output pdbqt file.
    Returns a list of tuples ordered by affinity: (affinity, index, ligand_id, scores)

    Inputs
    ------

        ligand_ids          : list of ligand ids
        ligand_names        : list of ligand names
        autodock_vina_paths : paths dictionary for autodock vina step
        output_path         : path to output directory

    Output
    ------

        ranking         : list of tuples ordered by affinity: (affinity, index, ligand_id, scores)
    """

    # List where best affinity for each ligand will be stored
    ranking = []

    # Go through all ligands
    for ligand_index in range(len(ligand_ids)):

        ligand_id = ligand_ids[ligand_index]
        ligand_name = ligand_names[ligand_index]

        # Find path to output pdbqt file with poses
        pdbqt_path = ligand_step_path(output_path, ligand_name, autodock_vina_paths['output_pdbqt_path'])

        # Find best affinity among different poses
        affinity = get_affinity(pdbqt_path = pdbqt_path)

        if affinity:
            ranking.append((affinity, ligand_index, ligand_id, [affinity]))

    # Sort list according to affinity (first element of tuple), most negative first
    ranking = sorted(ranking)

    return ranking

def get_ranking_gnina(
    ligand_ids: List, 
    ligand_names: List, 
    gnina_paths: Dict, output_path: str,
    rank_by: str, score_columns: Tuple
    ) -> List[Tuple]:
    """
    Takes the name of each ligand and finds its best scoring pose in the summary json file from gnina.
    Returns a list of tuples ordered by rank_by: (score, index, ligand_id, scores)

    Inputs
    ------

        ligand_ids      : list of ligand ids
        ligand_names    : list of ligand names
        gnina_paths     : paths dictionary for gnina step
        output_path     : path to output directory
        rank_by         : score to rank the ligands by
        score_columns   : scores reported for each ligand

    Output
    ------

        ranking         : list of tuples ordered by rank_by: (score, index, ligand_id, scores)
    """

    ranking = []

    for ligand_index in range(len(ligand_ids)):

        summary_path = ligand_step_path(output_path, ligand_names[ligand_index],
                                        gnina_paths['output_summary_path'])

        pose = get_gnina_score(summary_path = summary_path, rank_by = rank_by)

        if pose:
            ranking.append((pose[rank_by], ligand_index, ligand_ids[ligand_index],
                            [pose.get(column, '') for column in score_columns]))

    # Sort by rank_by, highest first for the CNN scores and lowest first for minimizedAffinity
    return sorted(ranking, reverse = GNINA_SORT_ORDERS[rank_by])

def dock_ligand(docking_engine: str, ligand_name: str, ligand_paths: Dict, ligand_prop: Dict, global_log) -> bool:
    '''
    Dock one ligand with the selected engine. A failure is logged and swallowed so the
    screening carries on with the next ligand.

    Inputs
    ------

        docking_engine  : docking engine to use, vina or gnina
        ligand_name     : name of the ligand, used for logging
        ligand_paths    : paths dictionary for this ligand
        ligand_prop     : properties dictionary for this ligand
        global_log      : global log

    Output
    ------

        success (bool): whether the docking produced its output files
    '''

    engine = DOCKING_ENGINES[docking_engine]
    step = engine['step']
    run_docking = gnina_run if docking_engine == 'gnina' else autodock_vina_run

    global_log.info(f"{step}: Docking the ligand")

    try:
        run_docking(**ligand_paths[step], properties=ligand_prop[step])
        return validate_step(*(ligand_paths[step][key] for key in engine['outputs']))
    except (Exception, SystemExit) as err:
        global_log.info(f"{step}: failed to dock ligand {ligand_name}")
        global_log.exception(f"{step}: {type(err).__name__}: {err}")
        return False

def clean_output(ligand_names: List, output_path: str):
    """
    Removes all ligand sub folders in the output folder

    Inputs
    ------

        ligand_names    : list of ligand names
        output_path     : path to output directory
    """

    # Remove all ligand subdirectories
    for name in ligand_names:
        ligand_path = os.path.join(output_path, name)

        if os.path.exists(ligand_path):
            shutil.rmtree(ligand_path)
    
def check_arguments(global_log,
                    ligand_lib_path,
                    structure_path,
                    input_pockets_zip,
                    pocket_num,
                    pocket_selection: Optional[str],
                    box_offset,
                    docking_engine: str,
                    gnina_cnn_scoring: Optional[str],
                    gnina_rank_by: Optional[str]
    ) -> None:
    """
    Check the arguments provided by the user and values of configuration file
    """

    # Check the ligand library path exists and it's a file
    if not os.path.exists(ligand_lib_path):
        global_log.warning(f"Ligand library file {ligand_lib_path} does not exist")
    elif not os.path.isfile(ligand_lib_path):
        global_log.warning(f"Ligand library path {ligand_lib_path} is not a file")

    # Check we have a structure file
    if not os.path.exists(structure_path):
        global_log.error(f"Structure file {structure_path} does not exist")
        sys.exit(1)
    elif not os.path.isfile(structure_path):
        global_log.error(f"Structure path {structure_path} is not a file")
        sys.exit(1)

    if pocket_selection is None:
        # Check we have a pocket selection file
        if not os.path.exists(input_pockets_zip):
            global_log.error(f"Pocket selection file {input_pockets_zip} does not exist")
            sys.exit(1)
        elif not os.path.isfile(input_pockets_zip):
            global_log.error(f"Pocket selection path {input_pockets_zip} is not a file")
            sys.exit(1)
            
        # Check we have a pocket number
        if pocket_num is None:
            global_log.warning(f"Pocket number not provided. Using the first pocket in the pocket selection file")
    else:
        if input_pockets_zip is not None:
            global_log.error(f"Cannot provide both pocket selection file and residues to dock to")
            sys.exit(1)
    
    if box_offset < 0:
        global_log.error(f"Box offset must be a positive number. Provided value: {box_offset}")
        sys.exit(1)
    elif box_offset > 5:
        global_log.warning(f"Box offset is {box_offset} angstroms. This may be unnecessarily large when docking to a selection of residues surrounding the binding site. Consider using a smaller value to improve performance.")

    # Check the docking engine is usable before any step runs
    if docking_engine == 'gnina':
        if gnina_cnn_scoring == 'none' and gnina_rank_by in ('CNNaffinity', 'CNNscore'):
            global_log.error(f"--gnina_rank_by {gnina_rank_by} needs CNN scores, not available with --gnina_cnn_scoring none")
            sys.exit(1)

def check_pdb(residues_path: str, global_log):
    """
    Checks the pdb is not empty and contains residues using biopython
    """
    
    # Check the residue file exists
    if not os.path.exists(residues_path):
        global_log.error(f"Residues file {residues_path} does not exist")
        sys.exit(1)
        
    # Load structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', residues_path)
    
    residues = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                residues += 1
    
    if residues == 0:
        global_log.error(f"No residues found in residues file {residues_path}")
        sys.exit(1)
    
    return residues


# YML construction
def config_contents(
        structure_path: str,
        input_pockets_zip: Optional[str],
        pocket_num: Optional[int],
        pocket_residues: List[Dict],
        box_offset: float,
        docking_engine: str,
        vina_bin: str,
        gnina_bin: str,
        gnina_scoring: Optional[str],
        gnina_cnn_scoring: Optional[str],
        gnina_cnn: Optional[str],
        gnina_seed: Optional[int],
        gnina_no_gpu: bool,
        cpus: int,
        exhaustiveness: int,
        restart: bool = False,
        debug: bool = False
    ) -> str:
    """
    Returns the contents of the YAML configuration file as a string.

    The YAML file contains the configuration for the protein preparation workflow.

    Parameters
    ----------
    debug: bool
        Keep temporary files for debugging purposes

    Returns
    -------
    str
        The contents of the YAML configuration file.
    """
    
    if structure_path:
        structure_path = os.path.abspath(structure_path)

    # fpocket_select does not run in the residue selection branch, but ConfReader resolves
    # every path in the config and it cannot resolve a null
    input_pockets_zip = os.path.abspath(input_pockets_zip) if input_pockets_zip else 'pockets.zip'

    # gnina docks the sdf directly, vina needs it converted to pdbqt
    ligand_format = 'sdf' if docking_engine == 'gnina' else 'pdbqt'

    return f"""
# Global properties (common for all steps)
global_properties:
  working_dir_path: output
  can_write_console_log: False
  restart: {to_yaml(restart)}
  remove_tmp: {to_yaml(not debug)}

# Section 1: Pocket selection and receptor preparation

# Step 0: Extract only the protein from the receptor (drops water/ligand/ion) - skippable
step0_extract_protein:
  tool: extract_molecule
  paths:
    input_structure_path: {structure_path}
    output_molecule_path: protein.pdb
  properties:
    molecule_type: protein                 # keep only protein (drops water/ligand/ion/na/dna/rna)

step1_fpocket_select:
  tool: fpocket_select
  paths:
    input_pockets_zip: {to_yaml(input_pockets_zip)}
    output_pocket_pdb: fpocket_cavity.pdb
    output_pocket_pqr: fpocket_pocket.pqr
  properties:
    pocket: {to_yaml(pocket_num)}

step1b_extract_residues:
  tool: extract_residues
  paths:
    input_structure_path: {structure_path}
    output_residues_path: pocket_residues.pdb
  properties:
    residues: {to_yaml(pocket_residues)}

step2_box:
  tool: box
  paths:
    input_pdb_path: dependency/step1_fpocket_select/output_pocket_pqr
    output_pdb_path: box.pdb 
  properties:
    offset: {to_yaml(box_offset)}
    box_coordinates: True

step3_str_check_add_hydrogens:
  tool: str_check_add_hydrogens
  paths:
    input_structure_path: {structure_path}
    output_structure_path: prep_receptor.pdbqt
  properties:
    charges: True
    mode: 'auto'                

# Section 2: Source each ligand and dock it to receptor

step4_babel_protonate:
  tool: babel_convert
  paths:
    input_path: ligand.smi
    output_path: ligand.{ligand_format}
  properties:
    coordinates: 3
    ph: 7.4

step4b_babel_convert:
  tool: babel_convert
  paths:
    input_path: ligand.sdf
    output_path: ligand.pdbqt
  properties:

step5_autodock_vina_run:
  tool: autodock_vina_run
  paths:
    input_ligand_pdbqt_path: dependency/step4_babel_protonate/output_path
    input_receptor_pdbqt_path: dependency/step3_str_check_add_hydrogens/output_structure_path
    input_box_path: dependency/step2_box/output_pdb_path
    output_pdbqt_path: output_vina.pdbqt
    output_log_path: output_vina.log
  properties:
    exhaustiveness: {to_yaml(int(exhaustiveness))}
    cpu: {to_yaml(int(cpus))}
    binary_path: {vina_bin}

# Docked with the box from step2 and not with gnina's own autobox, so that both engines
# search an identical box. The receptor is the step3 pdbqt, which gnina takes as given
# instead of reprocessing it with Open Babel on every ligand.
step5b_gnina_run:
  tool: gnina_run
  paths:
    input_ligand_path: dependency/step4_babel_protonate/output_path
    input_receptor_path: dependency/step3_str_check_add_hydrogens/output_structure_path
    input_box_path: dependency/step2_box/output_pdb_path
    output_sdf_path: output_gnina.sdf
    output_summary_path: output_gnina.json
    output_log_path: output_gnina.log
  properties:
    exhaustiveness: {to_yaml(int(exhaustiveness))}
    cpu: {to_yaml(int(cpus))}
    scoring: {to_yaml(gnina_scoring)}
    cnn_scoring: {to_yaml(gnina_cnn_scoring)}
    cnn: {to_yaml(gnina_cnn)}
    seed: {to_yaml(gnina_seed)}
    no_gpu: {to_yaml(gnina_no_gpu)}
    binary_path: {gnina_bin}

step6_babel_prepare_pose:
  tool: babel_convert
  paths:
    input_path: dependency/step5_autodock_vina_run/output_pdbqt_path
    output_path: output_vina.pdb
  properties:
"""

def create_config_file(output_path: str, 
                       **config_args) -> str:
    """
    Create a YAML configuration file for the workflow in the output path.
    Return the path to the configuration file.
    
    Parameters
    ----------
    output_path : str
        Path to the output folder
    config_args : dict
        Arguments to be used in the configuration file.
    
    Returns
    -------
    
    str
        Path to configuration file
    """
    
    config_path = os.path.join(output_path, 'config.yml')
    
    # Write the contents to the file
    with open(config_path, 'w') as f:
        f.write(config_contents(**config_args))

    print(f"Configuration file created at {config_path}.")
    
    return config_path
        
def virtual_screening(ligand_lib_path: str,
                      structure_path: str,
                      input_pockets_zip: str,
                      pocket_num: int = 1,
                      num_top_ligands: Optional[int] = None,
                      keep_poses: bool = False,
                      pocket_selection: Optional[str] = None,
                      box_offset: float = 5.0,
                      docking_engine: str = "vina",
                      vina_bin: str = "vina",
                      gnina_bin: str = "gnina",
                      gnina_scoring: Optional[str] = None,
                      gnina_cnn_scoring: Optional[str] = None,
                      gnina_cnn: Optional[str] = None,
                      gnina_rank_by: Optional[str] = None,
                      gnina_seed: Optional[int] = None,
                      gnina_no_gpu: bool = False,
                      cpus: int = 1,
                      exhaustiveness: int = 8,
                      debug: bool = False,
                      skip_extraction: bool = False,
                      restart: bool = True,
                      output_path: str = "output"
                      ) -> Tuple[Dict, Dict]:
    '''
    Main VS workflow. This workflow takes a ligand library, a pocket (defined by the output of a cavity analysis or some residues)
    and a receptor to screen the cavity using the ligand library (with AutoDock Vina or gnina).

    Inputs
    ------

        ligand_lib_path: 
            path to ligand library. Either a SMILES file or a SDF file
        structure_path: 
            path to receptor structure
        input_pockets_zip: 
            path to input pockets zip file
        pocket_num: 
            pocket number to be used from the input_pockets_zip. Default: 1
        num_top_ligands: 
            number of top ligands to be saved
        keep_poses: 
            keep poses of top ligands
        pocket_selection:
            list of residues to dock to. If provided, the input_pockets_zip and pocket_num will be ignored
        box_offset:
            extra distance (Angstroms) between the last residue atom and the box boundary. Default: 5.0
        docking_engine:
            docking engine to use, vina or gnina. Default: vina
        vina_bin:
            path to AutoDock Vina binary
        gnina_bin:
            path to gnina binary
        gnina_scoring:
            empirical scoring function used by gnina. None leaves gnina's own default
        gnina_cnn_scoring:
            where gnina uses the CNN. None leaves gnina's own default (rescore)
        gnina_cnn:
            CNN model used by gnina. None leaves gnina's own default (a 3 model ensemble)
        gnina_rank_by:
            gnina score to rank the ligands by. None resolves to CNNaffinity, or to
            minimizedAffinity when gnina_cnn_scoring is none
        gnina_seed:
            random seed for gnina, docking is stochastic
        gnina_no_gpu:
            force gnina to run on CPU even when a GPU is available
        cpus:
            number of cpus to use for each docking
        exhaustiveness:
            exhaustiveness of the docking
        debug:
            keep intermediate files for debugging
        skip_extraction:
            skip protein extraction from the receptor structure (keep cofactors/ligands/ions)
        restart:
            whether to restart the workflow from the last completed 
            step or start from the beginning.
        output_path:  
            path to output folder

    Outputs
    -------

        global_paths    : dictionary with all workflow paths
        global_prop     : dictionary with all workflow properties
    '''

    start_time = time.time()
    
    # Determine final output path
    output_path = fu.get_working_dir_path(output_path, restart=restart)

    # Initializing a global log file
    global_log, _ = fu.get_logs(path=output_path, light_format=True)
    
    # Check input files
    check_arguments(global_log,
                    ligand_lib_path,
                    structure_path,
                    input_pockets_zip,
                    pocket_num,
                    pocket_selection,
                    box_offset,
                    docking_engine,
                    vina_bin,
                    gnina_bin,
                    gnina_cnn_scoring,
                    gnina_rank_by)

    # Resolve the residue selection into the list of residues to extract (baked into the config).
    # Done here (not in config_contents) to keep the config builder a pure template.
    pocket_residues = []
    if pocket_selection is not None:
        u = mda.Universe(structure_path)
        unique_residues = set(atom.resid for atom in u.select_atoms(pocket_selection))
        pocket_residues = [{'res_id': str(res_id), 'model': '1'} for res_id in unique_residues]

    # Resolve the gnina ranking score, only minimizedAffinity exists without CNN scoring
    if gnina_rank_by is None:
        gnina_rank_by = 'minimizedAffinity' if gnina_cnn_scoring == 'none' else 'CNNaffinity'

    # Scores reported in scores.csv. Affinity is vina's affinity or gnina's minimizedAffinity
    if docking_engine == 'gnina' and gnina_cnn_scoring != 'none':
        score_columns = ('minimizedAffinity', 'CNNaffinity', 'CNNscore')
    elif docking_engine == 'gnina':
        score_columns = ('minimizedAffinity',)
    else:
        score_columns = ('Affinity',)

    # Create and load the configuration
    config_args = {
        'structure_path': structure_path,
        'input_pockets_zip': input_pockets_zip,
        'pocket_num': pocket_num,
        'pocket_residues' : pocket_residues,
        'box_offset' : box_offset,
        'docking_engine': docking_engine,
        'vina_bin': vina_bin,
        'gnina_bin': gnina_bin,
        'gnina_scoring': gnina_scoring,
        'gnina_cnn_scoring': gnina_cnn_scoring,
        'gnina_cnn': gnina_cnn,
        'gnina_seed': gnina_seed,
        'gnina_no_gpu': gnina_no_gpu,
        'cpus' : cpus,
        'exhaustiveness' : exhaustiveness,
        'restart' : restart,
        'debug' : debug}
    configuration_path = create_config_file(output_path, **config_args)

    conf = settings.ConfReader(configuration_path)
    conf.working_dir_path = output_path

    # Parsing the input configuration file (YAML);
    # Dividing it in global properties and global paths
    global_prop  = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # Docking step and biobb path keys for the selected engine
    engine = DOCKING_ENGINES[docking_engine]
    step = engine['step']

    # STEP 0: Extract protein from receptor (unless skipped) so downstream steps only see the protein
    receptor_path = structure_path
    if not skip_extraction:
        global_log.info("step0_extract_protein: Extracting protein from receptor structure")
        extract_molecule(**global_paths["step0_extract_protein"], properties=global_prop["step0_extract_protein"])
        receptor_path = global_paths["step0_extract_protein"]["output_molecule_path"]
        # Point downstream receptor consumers at the cleaned protein
        global_paths["step1b_extract_residues"]["input_structure_path"] = receptor_path
        global_paths["step3_str_check_add_hydrogens"]["input_structure_path"] = receptor_path

    # STEP 1: Select pocket or extract residues
    if pocket_selection is not None:
        
        # Extract residues from structure
        global_log.info("step1b_extract_residues: Extracting residues from structure")
        extract_residues(**global_paths["step1b_extract_residues"], properties=global_prop["step1b_extract_residues"])
        
        output_residues_path = global_paths["step1b_extract_residues"]['output_residues_path']
        
        # Check the output residues file exists and is not empty
        check_pdb(output_residues_path, global_log)

        # Modify step2_box paths to use residues
        global_paths['step2_box']['input_pdb_path'] = output_residues_path
    else:

        # Pocket selection from filtered list (input_pockets_zip / pocket already set in the config)
        global_log.info("step1_fpocket_select: Extract pocket cavity")
        fpocket_select(**global_paths["step1_fpocket_select"], properties=global_prop["step1_fpocket_select"])

    # STEP 2: Generate box around selected cavity or residues
    global_log.info("step2_box: Generating cavity box")
    box(**global_paths["step2_box"], properties=global_prop["step2_box"])

    # STEP 3: Prepare target protein for docking
    global_log.info("step3_str_check_add_hydrogens: Preparing target protein for docking")
    str_check_add_hydrogens(**global_paths["step3_str_check_add_hydrogens"], properties=global_prop["step3_str_check_add_hydrogens"]) 

    docking_start_time = time.time()

    # STEP 4-5: Prepare each ligand and dock it with the selected engine

    # Option 1: SDF library with protonated ligands
    if ligand_lib_path.endswith('.sdf'):

        ligand_names = []
        ligand_ids = []

        global_log.info(f"Reading ligand library in SDF format")
        ligand_supplier = pybel.readfile('sdf', ligand_lib_path)

        for index, ligand in enumerate(ligand_supplier, start=0):
            
            ligand_title = "".join(x for x in ligand.title if x.isalnum())
            
            # Create unique name from title and index
            ligand_name = f"{ligand_title}_{index}" if ligand_title else f"ligand_{index}"
            ligand_id = ligand_title if ligand_title else ""

            # Add ligand name to properties and paths
            ligand_prop = conf.get_prop_dic(prefix=ligand_name)
            ligand_paths = conf.get_paths_dic(prefix=ligand_name)
            
            # Update common paths
            ligand_paths[step][engine['receptor']] = global_paths[step][engine['receptor']]
            ligand_paths[step]['input_box_path'] = global_paths[step]['input_box_path']

            ligand_names.append(ligand_name)
            ligand_ids.append(ligand_id)

            # Create ligand subfolder
            ligand_folder = os.path.join(output_path, ligand_name)
            if not os.path.exists(ligand_folder):
                os.makedirs(ligand_folder)

            # gnina docks the sdf as written, vina needs the pdbqt from step4b
            if docking_engine == 'gnina':
                ligand_paths[step]['input_ligand_path'] = os.path.join(ligand_prop[step]['path'], 'ligand.sdf')
                ligand_sdf_path = ligand_paths[step]['input_ligand_path']
                step_folder = ligand_prop[step]['path']
            else:
                ligand_paths[step]['input_ligand_pdbqt_path'] = ligand_paths['step4b_babel_convert']['output_path']
                ligand_sdf_path = ligand_paths['step4b_babel_convert']['input_path']
                step_folder = ligand_prop['step4b_babel_convert']['path']

            # Create step folder
            if not os.path.exists(step_folder):
                os.makedirs(step_folder)

            # Write ligand to sdf file - writing the pdbqt file directly with pybel discards hydrogens.
            # Overwritten so a --restart run does not trip over the file it wrote last time
            ligand.write(format='sdf', filename=ligand_sdf_path, overwrite=True)

            # STEP 4: Convert ligand from sdf to pdbqt without adding hydrogens - gnina reads the sdf as is
            lastStep_successful = True
            if docking_engine == 'vina':
                global_log.info("step4b_babel_convert: Convert ligand to pdbqt format")
                try:
                    babel_convert(**ligand_paths['step4b_babel_convert'], properties = ligand_prop["step4b_babel_convert"])
                    lastStep_successful = validate_step(ligand_paths['step4b_babel_convert']['output_path'])
                except Exception:
                    global_log.info(f"step4b_babel_convert: Open Babel failed to convert ligand {ligand_name} to pdbqt format")
                    lastStep_successful = False

            # STEP 5: Docking
            if lastStep_successful:
                lastStep_successful = dock_ligand(docking_engine, ligand_name, ligand_paths, ligand_prop, global_log)

    # Option 2: SMILES library with ligands to be prepared
    elif ligand_lib_path.endswith('.smi'):

        global_log.info(f"Reading ligand library in SMILES format")
        ligand_ids, ligand_names = read_ligand_lib(ligand_lib_path)

        for ligand_id, ligand_name in zip(ligand_ids, ligand_names):

            # Add ligand name to properties and paths
            ligand_prop = conf.get_prop_dic(prefix=ligand_name)
            ligand_paths = conf.get_paths_dic(prefix=ligand_name)
            
            # Update common paths
            ligand_paths[step][engine['receptor']] = global_paths[step][engine['receptor']]
            ligand_paths[step]['input_box_path'] = global_paths[step]['input_box_path']

            # Create ligand subfolder
            ligand_folder = os.path.join(output_path, ligand_name)
            if not os.path.exists(ligand_folder):
                os.makedirs(ligand_folder)

            # Create step folder
            if not os.path.exists(ligand_prop['step4_babel_protonate']['path']):
                os.makedirs(ligand_prop['step4_babel_protonate']['path'])

            # Write smiles to file
            write_smiles(smiles = ligand_id, smiles_path = ligand_paths['step4_babel_protonate']['input_path'])

            # STEP 4: Convert ligand from smiles adding hydrogens at a certain pH - sdf for gnina, pdbqt for vina
            global_log.info("step4_babel_protonate: Prepare ligand for docking")
            try:
                babel_convert(**ligand_paths['step4_babel_protonate'], properties = ligand_prop["step4_babel_protonate"])
                lastStep_successful = validate_step(ligand_paths['step4_babel_protonate']['output_path'])
            except Exception:
                global_log.info(f"step4_babel_protonate: Open Babel failed to prepare ligand {ligand_name}")
                lastStep_successful = False

            # STEP 5: Docking
            if lastStep_successful:
                lastStep_successful = dock_ligand(docking_engine, ligand_name, ligand_paths, ligand_prop, global_log)

    else:

        global_log.error(f"Ligand library file {ligand_lib_path} should be in SDF or SMILES format")
        return

    # Rank ligands: find the best score for each ligand
    if docking_engine == 'gnina':
        global_log.info(f"Ranking ligands by {gnina_rank_by}, "
                        f"{'highest' if GNINA_SORT_ORDERS[gnina_rank_by] else 'lowest'} first")
        ranking = get_ranking_gnina(ligand_ids, ligand_names, global_paths[step], output_path,
                                    gnina_rank_by, score_columns)
    else:
        ranking = get_ranking(ligand_ids, ligand_names, global_paths[step], output_path)

    # Find top ligands and create csv file with ranking
    global_log.info("Create ranking and save poses for top ligands")
    ranking_path = os.path.join(output_path, "scores.csv")
    top_ligand_indices = save_ranking(ranking, num_top_ligands, ranking_path, score_columns)

    # STEP 6: extract poses for top ligands if requested
    if keep_poses:

        # Create poses folder
        poses_folder = os.path.join(output_path, "poses")
        if not os.path.exists(poses_folder):
            os.makedirs(poses_folder)

        # Iterate over top ligands
        for index in top_ligand_indices:

            ligand_name = ligand_names[index]

            # Add ligand name to properties and paths
            top_ligand_prop = conf.get_prop_dic(prefix=ligand_name)
            top_ligand_paths = conf.get_paths_dic(prefix=ligand_name)

            try:
                # gnina already writes sdf poses with their scores as SD data, a conversion
                # to pdb would drop them. Copied and not moved, the sdf is the docking output
                if docking_engine == 'gnina':
                    shutil.copy(top_ligand_paths[step]['output_sdf_path'],
                                os.path.join(poses_folder, f"{ligand_name}_poses.sdf"))
                    continue

                # Convert pose from pdbqt to pdb
                global_log.info("step6_babel_prepare_pose: Converting ligand pose to PDB format")
                babel_convert(**top_ligand_paths['step6_babel_prepare_pose'], properties=top_ligand_prop["step6_babel_prepare_pose"])

                # Move pose to final location
                # Pose path inside ligand subfolder
                pose_path = top_ligand_paths['step6_babel_prepare_pose']['output_path']
                # New pose path in poses folder
                new_pose_path = os.path.join(poses_folder, f"{ligand_name}_poses.pdb")
                # Move pose to new location
                shutil.move(pose_path, new_pose_path)

            except Exception:
                global_log.info(f"step6_babel_prepare_pose: failed to save pose for ligand {ligand_name}")
    
    # Show success rate of screening
    success_rate = round(len(ranking)/len(ligand_names)*100, 2)
    global_log.info(f"Success rate: {success_rate}%")

    if not debug:
        # Clean up the output folder 
        clean_output(ligand_names, output_path)

    # Save receptor used for docking in output_path (cleaned protein when extraction ran)
    shutil.copy(receptor_path, os.path.join(output_path, 'receptor.pdb'))

    # Save absolute path to ligand library in a text file
    with open(os.path.join(output_path, 'ligand_library.txt'), 'w') as file:
        file.write(os.path.abspath(ligand_lib_path))

    # Timing information
    elapsed_time = time.time() - start_time
    docking_elapsed_time = time.time() - docking_start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow name: Virtual Screening')
    global_log.info('  Docking engine: %s' % docking_engine)
    global_log.info('  Output path: %s' % output_path)
    global_log.info('  Config File: %s' % configuration_path)
    global_log.info('  Ligand library: %s' % ligand_lib_path)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('Docking time: %.1f minutes' % (docking_elapsed_time/60))
    global_log.info('')

    return global_paths, global_prop

def main():
    
    parser = argparse.ArgumentParser(description="Simple High-throughput virtual screening (HTVS) pipeline using BioExcel Building Blocks")

    parser.add_argument('--ligand_lib', dest='ligand_lib', type=str,
                        help="""Path to file with the ligand library. The format should be SMILES (.smi) or SDF (.sdf). For .smi files, 
                        one ligand per line is expected: 'smiles name'. For sdf files, the file may contain one or more ligands.""",
                        required=True)
    
    parser.add_argument('--structure_path', dest='structure_path', type=str,
                        help="""Path to file with target structure (PDB format). By default only the protein is kept (waters/ligands/ions are stripped);
                                use --skip_extraction to keep cofactors/ligands/ions. Beware that hydrogens will be added to the target structure
                                using biobb_structure_checking tool and a pH of 7""",
                        required=True)
    
    parser.add_argument('--input_pockets_zip', dest='input_pockets_zip', type=str,
                        help="Path to file with pockets in a zip file. Provide this path or a list of residues in the configuration file.",
                        required=False)

    parser.add_argument('--pocket_num', dest='pocket_num', type=int,
                        help="Pocket number to be used from the input_pockets_zip. Default: 1",
                        required=False)

    parser.add_argument('--num_top_ligands', dest='num_top_ligands', type=int,
                        help="Number of top ligands to be saved. Default: all successfully docked ligands",
                        required=False)

    parser.add_argument('--keep_poses', dest='keep_poses', action='store_true',
                        help="Save docking poses for top ligands. Default: False",  
                        required=False)

    parser.add_argument('--pocket_selection', dest='pocket_selection', type=str, default=None,
                        help="""Residue selection to define the pocket to give as an alternative to the pocket selection file.
                                The residue selection should be provided as a string following the MDAnalysis syntax, e.g. 'resid 37 49 112'.
                                This option is mutually exclusive with the pocket selection file.""",
                        required=False)
    
    parser.add_argument('--box_offset', dest='box_offset', type=float,
                        help="Extra distance (Angstroms) between the last residue atom and the box boundary. Default: 12",
                        required=False, default=5.0)

    parser.add_argument('--docking_engine', dest='docking_engine', type=str, choices=['vina', 'gnina'],
                        help="""Docking engine. 'vina' scores with the AutoDock Vina empirical function. 'gnina' uses the same sampling
                                but rescores the poses with a convolutional neural network, which needs a gnina binary and a biobb_vs
                                built from source. Default: vina""",
                        required=False, default='vina')

    parser.add_argument('--vina_bin', dest='vina_bin', type=str,
                        help="Path to AutoDock Vina binary. Default: vina",
                        required=False, default='vina')

    parser.add_argument('--gnina_bin', dest='gnina_bin', type=str,
                        help="Path to the gnina binary. gnina is not conda installable, download a release binary. Default: gnina",
                        required=False, default='gnina')

    parser.add_argument('--gnina_cnn_scoring', dest='gnina_cnn_scoring', type=str,
                        choices=['none', 'rescore', 'refinement', 'metrorescore', 'metrorefine', 'all'],
                        help="""Where gnina uses the CNN. 'none' is by far the fastest and gives no CNN scores at all, 'rescore' only
                                re-ranks the final poses, 'refinement' is around 10 times slower. Default: gnina's own default (rescore)""",
                        required=False, default=None)

    parser.add_argument('--gnina_cnn', dest='gnina_cnn', type=str,
                        help="""Name of the CNN model used by gnina, or a NAME_ensemble to use an ensemble of such models. 'fast' is a single model, around 3 times
                                faster than the default ensemble. Default: gnina's default 3 model ensemble""",
                        required=False, default=None)

    parser.add_argument('--gnina_scoring', dest='gnina_scoring', type=str,
                        choices=['default', 'vina', 'vinardo', 'ad4_scoring', 'dkoes_fast', 'dkoes_scoring', 'dkoes_scoring_old'],
                        help="Empirical scoring function used by gnina. vinardo often does better in virtual screening. Default: gnina's own default",
                        required=False, default=None)

    parser.add_argument('--gnina_rank_by', dest='gnina_rank_by', type=str,
                        choices=['CNNaffinity', 'CNNscore', 'minimizedAffinity'],
                        help="""gnina score to rank the ligands by. CNNaffinity (pK units, higher is better) is what ranks compounds,
                                CNNscore answers whether a pose is right, minimizedAffinity is kcal/mol (lower is better).
                                Default: CNNaffinity, or minimizedAffinity with --gnina_cnn_scoring none""",
                        required=False, default=None)

    parser.add_argument('--gnina_seed', dest='gnina_seed', type=int,
                        help="Random seed for gnina, docking is stochastic. Default: none",
                        required=False, default=None)

    parser.add_argument('--gnina_no_gpu', dest='gnina_no_gpu', action='store_true',
                        help="Force gnina to run on CPU even when a GPU is available. Default: False",
                        required=False, default=False)

    parser.add_argument('--cpus', dest='cpus', type=int,
                        help="Number of CPUs to use for each docking. Default: 1",
                        required=False, default=1)

    parser.add_argument('--exhaustiveness', dest='exhaustiveness', type=int,
                        help="Exhaustiveness of the docking. Number of runs for the sampling algorithm. Choose 4 to optimize speed and 8 to optimize accuracy. Default: 8",
                        required=False, default=8)
    
    parser.add_argument('-d', '--debug', dest='debug', action='store_true',
                        help="Keep intermediate files for debugging. Default: False",
                        required=False)

    parser.add_argument('--skip_extraction', dest='skip_extraction', action='store_true', default=False,
                        help="Skip protein extraction from the receptor structure (keep cofactors/ligands/ions). Default: False",
                        required=False)
    
    parser.add_argument('--restart', action='store_true',
                        help="Restart the workflow from the last completed step. Default: False",
                        required=False, default=False)

    parser.add_argument('--output', dest='output_path',
                        help="Output path (default: working_dir_path in YAML config file)",
                        required=False)
    
    args = parser.parse_args()

    virtual_screening(ligand_lib_path = args.ligand_lib,
                      structure_path = args.structure_path,
                      input_pockets_zip = args.input_pockets_zip,
                      pocket_num = args.pocket_num,
                      num_top_ligands = args.num_top_ligands,
                      keep_poses = args.keep_poses,
                      pocket_selection = args.pocket_selection,
                      box_offset = args.box_offset,
                      docking_engine = args.docking_engine,
                      vina_bin = args.vina_bin,
                      gnina_bin = args.gnina_bin,
                      gnina_scoring = args.gnina_scoring,
                      gnina_cnn_scoring = args.gnina_cnn_scoring,
                      gnina_cnn = args.gnina_cnn,
                      gnina_rank_by = args.gnina_rank_by,
                      gnina_seed = args.gnina_seed,
                      gnina_no_gpu = args.gnina_no_gpu,
                      cpus = args.cpus,
                      exhaustiveness = args.exhaustiveness,
                      debug = args.debug,
                      skip_extraction = args.skip_extraction,
                      restart = args.restart,
                      output_path = args.output_path)


if __name__ == '__main__':
    main()
