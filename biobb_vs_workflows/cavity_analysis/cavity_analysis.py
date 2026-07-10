#!/shared/work/BiobbWorkflows/envs/biobb_sp_cavity_analysis/bin/python

# Conversion of the BioExcel building blocks Python tutorials
# to a command line workflow with two files: Python Script and YAML input configuration file

# Importing all the needed libraries
import os
import time
import glob
import argparse
import csv
import yaml
import json
import shutil
import logging
import numpy as np
from pathlib import Path
from typing import Dict, Tuple, Any, Optional, Literal

import MDAnalysis as mda

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu

from biobb_vs_workflows.common import to_yaml

from biobb_gromacs.gromacs.make_ndx import make_ndx
from biobb_analysis.gromacs.gmx_cluster import gmx_cluster
from biobb_analysis.ambertools.cpptraj_convert import cpptraj_convert
from biobb_structure_utils.utils.extract_model import extract_model
from biobb_vs.fpocket.fpocket_run import fpocket_run
from biobb_vs.fpocket.fpocket_filter import fpocket_filter

# Cluster count limit and retry policy for step4_gmx_cluster
# (too many clusters -> combined centroid PDB can exceed step5's atom limit, and extraction gets slow)
MAX_CLUSTER_ATOMS = 1000000
MAX_CLUSTERS = 100
MAX_CLUSTER_RETRIES = 10
CUTOFF_INCREASE_FACTOR = 1.5

def is_gromacs_format(traj_path: Optional[str]) -> bool:
    """
    Checks if the trajectory is in a GROMACS-compatible format 
    (xtc, trr, cpt, g96, gro, pdb, tng)

    Inputs
    ------

        traj_path: 
            path to the trajectory file
    
    Output
    -------

        bool: True if the trajectory is in a GROMACS-compatible format, False otherwise
    """

    # List of GROMACS-compatible formats
    gromacs_formats = ['.xtc', '.trr', '.cpt', '.g96', '.gro', '.pdb', '.tng']

    if traj_path:
        # Check if the trajectory file is in a GROMACS-compatible format
        if any([traj_path.endswith(format) for format in gromacs_formats]):
            return True
        else:
            return False
    else:
        raise ValueError("traj_path is None, please provide a valid trajectory path")

def get_clusters_population(log_path: str, output_path: str, global_log) -> list:
    '''
    Reads the centroids' ID and populations from the log, sorts the clusters by population 
    in descending order and writes the sorted list to a JSON file in the output folder.

    Inputs
    ------

        log_path          : path to log file of clustering step
        output_path       : path to the output folder where the JSON file will be saved
        global_log        : global log object

    Outputs
    -------

        clusters_population (list<tuple>): list with tuples containing (population, cluster_ID) sorted by population
    '''

    # Open log file if it exists
    if os.path.exists(log_path):
        file = open(log_path)
    else:
        global_log.error("Clustering log file not found")
        return

    # Read file
    csv_reader = csv.reader(file)

    cluster_ids = []
    populations = []

    # Parse the log file
    start = False
    for row in csv_reader:

        # When we have reached clustering information
        if start: 
            # Split the whole row into columns using '|' as separators
            col = row[0].split('|')

            # If the first column contains something more than whitespace
            if len(col[0].strip()): 
                # Save cluster id and population
                cluster_ids.append(col[0].strip())
                populations.append(col[1].strip().split()[0])

        # We have reached beginning of clustering information, start saving
        if len(row) and row[0].startswith('cl.'):
            start = True
        
    # Close log file
    file.close()
    
    # Convert to integers
    populations = [int(x) for x in populations]
    cluster_ids = [int(x) for x in cluster_ids]

    # Sort clusters by population
    clusters_population = sorted(zip(populations, cluster_ids), reverse=True)

    ordered_clusters = []

    # Save sorted cluster ids and populations in dictionaries and populations in list
    for population,index in clusters_population:
        ordered_clusters.append({
            'cluster': index,
            'population': population
            })
    
    # Save list in JSON file
    with open(os.path.join(output_path, "clusters.json"), 'w') as outfile:
        json.dump(ordered_clusters, outfile)

    outfile.close()

    # If the number of clusters is very large, issue a warning - the user might want to increase the RMSD cutoff
    if len(cluster_ids)>MAX_CLUSTERS:
        global_log.warning(f"   Warning: Large number of clusters found. Consider increasing the RMSD cutoff.")
        global_log.warning(f"   Warning: Number of clusters: {len(cluster_ids)}")
        global_log.warning(f"   Warning: Number of clusters with more than 1 member: {len([x for x in populations if x > 1])}")
    
    # Return list with cluster population and id sorted by population
    return clusters_population

def create_summary(cluster_names, cluster_populations, cluster_filtered_pockets, global_paths, output_path):
    '''
    Creates 3 sorted summary files with information for each pocket found and filtered for each model. 
    The summary file is a YAML file.

    
    Inputs
    ------

        cluster_names            (list): List with names of clusters
        cluster_populations      (list): List with tuples containing (population, cluster_ID) ordered by population
        cluster_filtered_pockets (dict): Dictionary with filtered pocket's IDs for each cluster
        global_paths             (dict): global paths dictionary
        output_path               (str): path to output folder
    
    Output
    ------

        Info dumped to yaml summary files
    '''

    # Find step names
    cavity_analysis_folder = 'step6_cavity_analysis'

    # Find unfiltered summary file name
    pockets_summary_filename = Path(global_paths[cavity_analysis_folder]['output_summary']).name

    # Dictionary where all available cluster_summary dictionaries will be saved
    global_summary = {}

    # For each cluster
    for cluster_index, cluster_name in enumerate(cluster_names):
        
        # Find list of filtered pocket IDs for this cluster
        filtered_pocket_names = cluster_filtered_pockets[cluster_name]

        # If any pockets are found 
        if len(filtered_pocket_names) > 0:

            # Dictionary with information for this cluster
            cluster_summary = {}

            # If clustering was done externally we might not have this information 
            if cluster_populations is not None:
                # Save population of cluster
                cluster_summary.update({'population' : cluster_populations[cluster_index][0]})

            # Path to this cluster's summary with information for all pocket found
            summary_path = os.path.join(output_path, cluster_name, cavity_analysis_folder, pockets_summary_filename)

            # Load all pockets summary as dictionary
            with open(summary_path) as json_file:
                all_pockets_summary = json.load(json_file)

            # For each filtered pocket
            for pocket_name in filtered_pocket_names:

                # Save entry with pocket information
                cluster_summary.update({pocket_name : all_pockets_summary[pocket_name]})

            # Save pocket IDs
            cluster_summary.update({'pockets' : filtered_pocket_names})

            # Update global_summary
            global_summary.update({cluster_name : cluster_summary})

    # Sort models by 3 criteria (volume, druggability score, score)
    sorted_pockets_by_volume, sorted_pockets_by_drug_score, sorted_pockets_by_score = sort_summary(global_summary)
    
    # Create file names for sorted summary files
    volume_summary_path = os.path.join(output_path, f"summary_by_volume.yml")
    drug_score_summary_path = os.path.join(output_path, f"summary_by_drug_score.yml")
    score_summary_path = os.path.join(output_path, f"summary_by_score.yml")

    # Write the sorted pockets by volume to a YAML file
    with open(volume_summary_path, 'w') as f:
        yaml.dump(sorted_pockets_by_volume, f, sort_keys = False)

    # Write the sorted pockets by druggability score to a YAML file
    with open(drug_score_summary_path, 'w') as f:
        yaml.dump(sorted_pockets_by_drug_score, f, sort_keys = False)
    
    # Write the sorted pockets by score to a YAML file
    with open(score_summary_path, 'w') as f:
        yaml.dump(sorted_pockets_by_score, f, sort_keys = False)

    return 

def filter_residue_com(input_pockets_zip: str, input_pdb_path: str, output_filter_pockets_zip: str, properties: dict, global_log):
    """
    Function that filters pockets by the distance of their center of mass to a group of residues.

    Inputs
    ------

        input_pockets_zip         (str): path to input pockets zip file
        input_pdb_path            (str): path to input pdb with the pocket model (pdb of the receptor)
        output_filter_pockets_zip (str): path to filtered pockets zip file
        properties               (dict): dictionary with properties for this step
        global_log                (log): global log object
    
    Output
    ------

        filtered_pocket_IDs (list(str)): list with pocket IDs that passed the filter
    """

    # To return and use to create the summary file
    filtered_pocket_IDs = []

    # Create step folder
    fu.create_dir(properties['path'])

    # Check if input pockets zip file exists
    if not os.path.exists(input_pockets_zip):
        global_log.warning("Input pockets zip file not found, previous step didn't find any pockets or failed")
        return filtered_pocket_IDs

    # Check if step should run
    if not properties['run_step']:

        # Copy input pockets zip file to output filtered pockets zip file
        shutil.copyfile(input_pockets_zip, output_filter_pockets_zip) 

        # Find list of filtered pocket IDs 
        filtered_pocket_IDs = get_pockets_IDs(output_filter_pockets_zip, properties, global_log)

        # Log warning
        global_log.warning("    Skipping step because run_step = False")

        return filtered_pocket_IDs
    
    # Extract all pockets in step folder
    pocket_paths = fu.unzip_list(zip_file=input_pockets_zip, dest_dir=properties['path'])

    # If no pockets are found, return
    if len(pocket_paths) == 0:
        global_log.warning("No pockets found after filtering in previous step")
        return filtered_pocket_IDs
    
    # Load input pdb
    model_universe = mda.Universe(input_pdb_path)

    # Select the residues of interest, e.g. residue number 42
    residue_selection = model_universe.select_atoms(properties['residue_selection'])

    # Compute the center of mass of the selected residues
    residue_com = np.array(residue_selection.center_of_mass())

    # Save paths to pqr files in another list
    pockets_pqr_paths = []
    
    # For each pocket
    for pocket_path in pocket_paths:
        # If path is a pqr file, append to list
        if pocket_path.endswith('.pqr'):
            pockets_pqr_paths.append(pocket_path)

    # Save pockets that pass the distance filter in another list
    filtered_pocket_paths = []

    # Iterate over all pqr files
    for pocket_pqr_path in pockets_pqr_paths:

        # Find pocket ID as file name without "_vert.pqr"
        pocket_ID = Path(pocket_pqr_path).stem.replace("_vert", "")

        # Load pocket
        pocket_universe = mda.Universe(pocket_pqr_path)

        # Select all atoms in pocket
        filtering_selection = pocket_universe.select_atoms('all')

        # Compute the center of mass of the pocket
        pocket_com = np.array(filtering_selection.center_of_mass())

        # Compute distance between pocket and residue center of mass using numpy
        distance = np.linalg.norm(pocket_com - residue_com)

        # Save pocket if distance is smaller than threshold
        if distance < properties['distance_threshold']:

            # Save pocket
            filtered_pocket_IDs.append(pocket_ID)
            filtered_pocket_paths.append(pocket_pqr_path)

    # If no pockets are found, return
    if len(filtered_pocket_paths) == 0:
        
        # Erase all the pockets remaining in the step folder
        fu.rm_file_list(file_list=pocket_paths)

        # Log warning
        global_log.warning("No pockets found after filtering")

        return filtered_pocket_IDs
    
    # Zip filtered pockets
    fu.zip_list(zip_file=output_filter_pockets_zip, file_list=filtered_pocket_paths)

    # Erase all the pockets remaining in the step folder
    fu.rm_file_list(file_list=pocket_paths)

    return filtered_pocket_IDs

def sort_summary(pockets_summary: dict):
    """
    Function that reads the dictionary with all models and pockets and sorts the models by:

        1. Volume of the largest pocket in that model
        2. Druggability score of the best pocket in that model
        3. Score of the best pocket in that model
    
    It returns the 3 dictionaries with sorted models and their pockets.

    Inputs
    ------

        pockets_summary (dict): dictionary with all models and pockets
    
    Outputs
    -------

        sorted_pockets_by_volume     (dict): dictionary with models sorted by volume
        sorted_pockets_by_drug_score (dict): dictionary with models sorted by druggability score
        sorted_pockets_by_score      (dict): dictionary with models sorted by score
    
    """

    # Sort the pockets by volume
    sorted_pockets_by_volume = dict(sorted(pockets_summary.items(), key = lambda x: largest_volume(x[1]), reverse = True))

    # Sort the pockets by druggability score
    sorted_pockets_by_drug_score = dict(sorted(pockets_summary.items(), key = lambda x: highest_drug_score(x[1]), reverse = True))

    # Sort the pockets by score
    sorted_pockets_by_score = dict(sorted(pockets_summary.items(), key = lambda x: highest_score(x[1]), reverse = True))
    
    return sorted_pockets_by_volume, sorted_pockets_by_drug_score, sorted_pockets_by_score

def largest_volume(model: dict):
    """
    Function to sort the pockets by volume
    """

    # Find the largest volume
    largest_volume = 0
    for pocket in model["pockets"]:
        if model[pocket]["volume"] > largest_volume:
            largest_volume = model[pocket]["volume"]

    return largest_volume

def highest_drug_score(model: dict):
    """
    Function to sort the pockets by druggability score
    """

    # Find the highest druggability score
    highest_drug_score = 0
    for pocket in model["pockets"]:
        if model[pocket]["druggability_score"] > highest_drug_score:
            highest_drug_score = model[pocket]["druggability_score"]

    return highest_drug_score

def highest_score(model: dict):
    """
    Function to sort the pockets by score
    """

    # Find the highest score
    highest_score = 0
    for pocket in model["pockets"]:
        if model[pocket]["score"] > highest_score:
            highest_score = model[pocket]["score"]

    return highest_score

def get_pockets_IDs(input_pockets_zip: str, properties: dict, global_log):
    """
    Function that retrieves all the pocket IDs from the pockets zip file.

    Inputs
    ------

        input_pockets_zip (str): path to input pockets zip file
        properties       (dict): dictionary with properties for this step
        global_log        (log): global log object
    
    Output
    ------

        filtered_pocket_IDs (list(str)): list with all pocket IDs found in the pockets zip file
    """

    # To return and use to create the summary file
    filtered_pocket_IDs = []

    # Check if input pockets zip file exists
    if not os.path.exists(input_pockets_zip):
        global_log.warning("Input pockets zip file not found, last step didn't find any pockets or failed")
        return filtered_pocket_IDs

    # Extract all pockets in step folder
    pocket_paths = fu.unzip_list(zip_file=input_pockets_zip, dest_dir=properties['path'])

    # If no pockets are found, return
    if len(pocket_paths) == 0:
        global_log.warning("No pockets found after filtering in last step")
        return filtered_pocket_IDs 
    
    # Save paths to pqr files in another list
    pockets_pqr_paths = []

    # For each pocket
    for pocket_path in pocket_paths:
        # If path is a pqr file, append to list
        if pocket_path.endswith('.pqr'):
            pockets_pqr_paths.append(pocket_path)
    
    # Iterate over all pqr files
    for pocket_pqr_path in pockets_pqr_paths:

        # Find pocket ID as file name without "_vert.pqr"
        pocket_ID = Path(pocket_pqr_path).stem.replace("_vert", "")
        filtered_pocket_IDs.append(pocket_ID)
    
    # Erase all the pockets remaining in the step folder
    fu.rm_file_list(file_list=pocket_paths)

    return filtered_pocket_IDs

def check_arguments(global_log: logging.Logger, 
                    traj_path: Optional[str], 
                    top_path: Optional[str], 
                    structures_path: Optional[str]
    ):
    """
    Check the arguments provided by the user
    
    Parameters
    ----------
    
    global_log : logging.Logger
        Global log object
    traj_path : Optional[str]
        Path to the trajectory file
    top_path : Optional[str]
        Path to the topology file
    structures_path : Optional[str]
        Path to a folder with the representative structures in pdb format
    """

    # If the user doesn't provide traj_path and top_path or structures_path 
    if (None in [traj_path, top_path]) and structures_path is None:

        global_log.error("ERROR: traj_path and top_path or structures_path must be provided")
        raise SystemExit

    # If the user provides traj_path and top_path and structures_path -> exit
    if (None not in [traj_path, top_path]) and structures_path is not None:
        global_log.error("ERROR: traj_path, top_path and structures_path are provided, provide either traj_path and top_path or structures_path")
        raise SystemExit

    # If the user provides traj_path and not top_path -> exit
    if traj_path is not None and top_path is None:
        global_log.error("ERROR: top_path must be provided if traj_path is provided")
        raise SystemExit
    
    # If the user provides top_path and not traj_path -> exit
    if top_path is not None and traj_path is None:
        global_log.error("ERROR: traj_path must be provided if top_path is provided")
        raise SystemExit
    
    # If the user provides traj_path and it doesn't exist -> exit
    if traj_path is not None and not os.path.exists(traj_path):
        global_log.error("ERROR: traj_path doesn't exist")
        raise SystemExit

    # If the user provides top_path and it doesn't exist -> exit
    if top_path is not None and not os.path.exists(top_path):
        global_log.error("ERROR: top_path doesn't exist")
        raise SystemExit
    
    # If the user provides structures_path and it doesn't exist -> exit
    if structures_path is not None and not os.path.exists(structures_path):
        global_log.error("ERROR: structures_path doesn't exist")
        raise SystemExit

# YML construction
def config_contents(
    gmx_bin: Optional[str],
    traj_path: Optional[str],
    top_path: Optional[str],
    clustering_method: Optional[Literal['linkage', 
                                        'jarvis-patrick', 
                                        'monte-carlo', 
                                        'diagonalization', 
                                        'gromos']],
    clustering_cutoff: Optional[float],
    filtering_selection: Optional[str],
    distance_threshold: Optional[float],
    restart: bool = False

    ) -> str:
    """
    Returns the contents of the YAML configuration file as a string.
    
    The YAML file contains the configuration for the protein preparation workflow.
    
    Paramters
    ---------
    
    gmx_bin: str
        Path to GROMACS binary
    traj_path: str
        Path to the trajectory file
    top_path: str
        Path to the topology file
    filtering_selection: str
        Residue selection to filter pockets by distance to center of mass
    distance_threshold: float
        Distance threshold to filter pockets by distance to center of mass
    clustering_method: str
        Clustering method to use (linkage, jarvis-patrick, monte-carlo, diagonalization, gromos)
    clustering_cutoff: float
        Clustering cutoff to use for the clustering method
    Returns
    -------
    str
        The contents of the YAML configuration file.
    """
    
    if traj_path:
        traj_path = os.path.abspath(traj_path)
        
    if top_path:
        top_path = os.path.abspath(top_path)
        
    if gmx_bin == None:
        gmx_bin = "gmx"

    # Quoted MDAnalysis selection string when given, null otherwise (avoids emitting "None")
    # Only run the residue center-of-mass filter when a selection is provided.
    if filtering_selection is not None:
        residue_selection_property = f'residue_selection: "{filtering_selection}"'
        filter_residue_com_run_step = True
    else:
        residue_selection_property = "residue_selection: null"
        filter_residue_com_run_step = False

    return f"""
# Global properties (common for all steps)
global_properties:                                # Wether to use GPU support or not
  working_dir_path: output                        # Workflow default output directory
  can_write_console_log: False                    # Verbose writing of log information
  restart: {to_yaml(restart)}                     # Skip steps already performed
  remove_tmp: True                                # Do not execute steps if output files are already created

# Step 0: Convert from Amber to Gromacs compatible format
# Optional step (will be executed if the trajectory is not in a Gromacs-compatible format)
step0_convert_amber_traj:
  tool: cpptraj_convert
  paths:
    input_traj_path: {traj_path}                                        # Amber compatible trajectory file
    input_top_path: {top_path}                                          # topology file
    output_cpptraj_path: trajectory.xtc
  properties:
    mask: "all-atoms"

# Step 3: Create index file to select the atoms for the RMSD calculation
step3A_rmsd_calculation_ndx:
  tool: make_ndx 
  paths:
    input_structure_path: {top_path}
    output_ndx_path: index.ndx
  properties:
    selection: "System"
    binary_path: {gmx_bin}      

# Steps 4-5: Cluster trajectory and extract centroids pdb
step4_gmx_cluster:
  tool: gmx_cluster
  paths:
    input_traj_path: {traj_path}
    input_structure_path: {top_path}
    input_index_path: dependency/step3A_rmsd_calculation_ndx/output_ndx_path
    output_pdb_path: output.cluster.pdb
    output_cluster_log_path: output.cluster.log
    output_rmsd_cluster_xpm_path: output.rmsd-clust.xpm
    output_rmsd_dist_xvg_path: output.rmsd-dist.xvg
  properties:
    fit_selection: Protein       
    output_selection: Protein        
    dista: False
    method: {to_yaml(clustering_method)}
    cutoff: {to_yaml(clustering_cutoff)}
    nofit: False
    binary_path: {gmx_bin} 

step5_extract_models:
  tool: extract_model
  paths:
    input_structure_path: dependency/step4_gmx_cluster/output_pdb_path
    output_structure_path: cluster.pdb     
  properties:

# Step 6-8: Cavity analysis with fpocket on centroids + filtering
step6_cavity_analysis:
  tool: fpocket_run
  paths:
    input_pdb_path: dependency/step5_extract_models/output_structure_path
    output_pockets_zip: all_pockets.zip
    output_summary: summary.json
  properties:
    min_radius: 3
    max_radius: 6
    num_spheres: 35
    sort_by: druggability_score

step7_filter_cavities:
  tool: fpocket_filter
  paths:
    input_pockets_zip: dependency/step6_cavity_analysis/output_pockets_zip
    input_summary: dependency/step6_cavity_analysis/output_summary
    output_filter_pockets_zip: filtered_pockets.zip
  properties:
    score: [0.4, 1]
    druggability_score: [0.4, 1]
    volume: [200, 5000]

step8_filter_residue_com:
  paths: 
    input_pockets_zip: dependency/step7_filter_cavities/output_filter_pockets_zip
    input_pdb_path: dependency/step5_extract_models/output_structure_path
    output_filter_pockets_zip: filtered_pockets.zip
  properties:
    {residue_selection_property}      # MDAnalysis selection string
    distance_threshold: {to_yaml(distance_threshold)}     # Distance threshold in Angstroms (6-8 are reasonable values if the residue/s are part of the pocket)
    run_step: {to_yaml(filter_residue_com_run_step)}      # Run only when a filtering_selection is given
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
        
# Main workflow   
def cavity_analysis(traj_path: Optional[str], 
                    top_path: Optional[str], 
                    structures_path: Optional[str], 
                    distance_threshold: Optional[float],  
                    filtering_selection: Optional[str], 
                    num_clusters: Optional[int],
                    clustering_method: Optional[Literal['linkage', 
                                                        'jarvis-patrick', 
                                                        'monte-carlo', 
                                                        'diagonalization', 
                                                        'gromos']],
                    clustering_cutoff: Optional[float],
                    gmx_bin: Optional[str],
                    restart: bool,
                    output_path: Optional[str]
                    ) -> Tuple[str, Dict[str, Any]]:
    '''
    Main clustering and cavity analysis workflow. This workflow clusters a given trajectory and 
    analyzes the cavities of the most representative structures. Then filters the cavities 
    according to a pre-defined criteria and outputs the pockets that passed the filter.

    Inputs
    ------

        traj_path:  
            path to trajectory file
        top_path:  
            path to topology file
        structures_path:  
            path to the folder with the most representative structures in pdb format from an external clustering
        distance_threshold:  
            distance threshold to filter pockets by distance to center of mass 
        filtering_selection:
            residue selection to filter pockets by distance to center of mass 
        num_clusters:
            number of clusters to extract from the trajectory and analyze with fpocket
        clustering_method:
            clustering method to use (linkage, jarvis-patrick, monte-carlo, diagonalization, gromos)
        clustering_cutoff:
            clustering cutoff to use for the clustering method
        gmx_bin:
            path to GROMACS binary
        restart:
            whether to restart the workflow from the last completed 
            step or start from the beginning.
        output_path:  
            path to output folder

    Outputs
    -------

        /output folder

        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary will all workflow properties
    '''

    start_time = time.time()

    # Determine final output path
    output_path = fu.get_working_dir_path(output_path, restart=restart)

    # Initialize a global log file
    global_log, _ = fu.get_logs(path=output_path, light_format=True)
    
    # Check input files
    check_arguments(global_log, traj_path, top_path, structures_path)

    # Create and load the configuration
    config_args = {
        'gmx_bin': gmx_bin,
        'traj_path': traj_path,
        'top_path': top_path,
        'clustering_method': clustering_method,
        'clustering_cutoff': clustering_cutoff,
        'filtering_selection': filtering_selection,
        'distance_threshold': distance_threshold,
        'restart': restart
    }
    configuration_path = create_config_file(output_path, **config_args)

    conf = settings.ConfReader(configuration_path)
    conf.working_dir_path = output_path

    # Parsing the input configuration file (YAML);
    # Dividing it in global properties and global paths
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # If clustering is not given externally -> cluster the input trajectory
    if structures_path is None:

        # If the trajectory is not in a GROMACS-compatible format, convert it
        if not is_gromacs_format(traj_path):

            # STEP 0: Convert the trajectory to xtc format
            global_log.info("step0_convert_amber_traj: Converting AMBER trajectory to xtc format")
            cpptraj_convert(**global_paths['step0_convert_amber_traj'], properties=global_prop['step0_convert_amber_traj'])

            # Change subsequent traj path 
            global_paths['step4_gmx_cluster']['input_traj_path'] = global_paths['step0_convert_amber_traj']['output_cpptraj_path']

        # STEP 3A: Create index file for rmsd calculation
        global_log.info(f"Paths: {global_paths['step3A_rmsd_calculation_ndx']}")
        global_log.info("step3A_rmsd_calculation_ndx: Creation of index file")
        make_ndx(**global_paths['step3A_rmsd_calculation_ndx'], properties=global_prop['step3A_rmsd_calculation_ndx'])

        # STEP 4: Cluster trajectory with gmx_cluster
        global_log.info("step4_gmx_cluster: Clustering structures from the trajectory")
        gmx_cluster(**global_paths["step4_gmx_cluster"], properties=global_prop["step4_gmx_cluster"])

        # Save centroid IDs and populations in JSON file
        global_log.info( "step4_gmx_cluster: Reading clustering outcome, generating clusters JSON file")
        cluster_populations = get_clusters_population(log_path = global_paths["step4_gmx_cluster"]['output_cluster_log_path'],
                                                      output_path = global_prop["step4_gmx_cluster"]['path'],
                                                      global_log = global_log)
        
        # Find the number of atoms in the combined centroid PDB file
        combined_centroid_universe = mda.Universe(global_paths["step4_gmx_cluster"]['output_pdb_path'])
        num_atoms = len(combined_centroid_universe.atoms)

        # If too many clusters were found, the combined centroid PDB from step4 can exceed step5's atom
        # limit and extraction becomes very slow. Increase the cutoff and re-cluster, up to a retry limit.
        retry = 0
        while num_atoms > MAX_CLUSTER_ATOMS and retry < MAX_CLUSTER_RETRIES:
            retry += 1
            old_cutoff = global_prop["step4_gmx_cluster"]["cutoff"]
            new_cutoff = old_cutoff * CUTOFF_INCREASE_FACTOR

            global_log.warning(
                f"step4_gmx_cluster: {num_atoms} atoms found between all cluster (> {MAX_CLUSTER_ATOMS} limit). "
                f"Too many clusters can make step5_extract_models fail (atom count limit on the combined "
                f"centroid PDB) or become very slow. Increasing clustering cutoff {old_cutoff} -> {new_cutoff} "
                f"and re-clustering (attempt {retry}/{MAX_CLUSTER_RETRIES})."
            )

            global_prop["step4_gmx_cluster"]["cutoff"] = new_cutoff
            # Force re-run: outputs already exist from the previous attempt, bypass the restart-skip
            global_prop["step4_gmx_cluster"]["restart"] = False

            gmx_cluster(**global_paths["step4_gmx_cluster"], properties=global_prop["step4_gmx_cluster"])
            cluster_populations = get_clusters_population(log_path = global_paths["step4_gmx_cluster"]['output_cluster_log_path'],
                                                          output_path = global_prop["step4_gmx_cluster"]['path'],
                                                          global_log = global_log)

            # Find the number of atoms in the combined centroid PDB file
            combined_centroid_universe = mda.Universe(global_paths["step4_gmx_cluster"]['output_pdb_path'])
            num_atoms = len(combined_centroid_universe.atoms)

        if num_atoms > MAX_CLUSTER_ATOMS:
            global_log.warning(
                f"step4_gmx_cluster: still {len(cluster_populations)} clusters after {retry} retries "
                f"(final cutoff={global_prop['step4_gmx_cluster']['cutoff']}). Proceeding anyway - "
                f"step5_extract_models may fail or be slow."
            )

        if num_clusters:
            # Number of clusters: minimum between number of clusters requested and number of clusters obtained
            num_clusters = min(num_clusters, len(cluster_populations))
        else:
            # Number of clusters: number of clusters obtained
            num_clusters = len(cluster_populations)

        # Cluster names are the cluster IDs
        cluster_names = [str(cluster_populations[i][1]) for i in range(num_clusters)]
    
    # If clustering is given externally
    else:

        # Obtain the full sorted list of pdb files from clustering path
        # If the clustering path is a file, we assume it is a single pdb file
        if os.path.isfile(structures_path):
            global_log.info("External clustering file provided")
            pdb_paths = [structures_path]
        else:
            pdb_paths = sorted(glob.glob(os.path.join(structures_path,"*.pdb")))

        # Population information will not be available in this case
        cluster_populations = None

        # Number of clusters: number of pdb files
        num_clusters = len(pdb_paths)

        # Cluster names are the pdb file names
        cluster_names = [Path(pdb_path).stem for pdb_path in pdb_paths]
    
    global_log.info(f"Number of models to analyze: {num_clusters}")

    # Dictionary to save the filtered pocket IDs for each model
    cluster_filtered_pockets = {}

    # For representative structure (model)
    for cluster_index, cluster_name in enumerate(cluster_names):

        # Create sub folder for the model
        cluster_prop = conf.get_prop_dic(prefix=cluster_name)
        cluster_paths = conf.get_paths_dic(prefix=cluster_name)
            
        # If clustering was done here, extract the model from the clustering results
        if structures_path is None:

            # Update input structures path and model index
            cluster_paths['step5_extract_models']['input_structure_path'] = global_paths['step4_gmx_cluster']['output_pdb_path']
            if cluster_populations is not None:
                cluster_prop['step5_extract_models']['models'] = [cluster_populations[cluster_index][1]]

            # STEP 5: Extract one model from the input structures path
            extract_model(**cluster_paths['step5_extract_models'], properties=cluster_prop['step5_extract_models'])

        # If clustering was done externally, just update the input pdb path
        else:

            cluster_paths['step6_cavity_analysis']['input_pdb_path'] = pdb_paths[cluster_index]
            cluster_paths['step8_filter_residue_com']['input_pdb_path'] = pdb_paths[cluster_index]
        
        # STEP 6: Cavity analysis
        global_log.info("step6_cavity_analysis: Compute protein cavities using fpocket")
        fpocket_run(**cluster_paths['step6_cavity_analysis'], properties=cluster_prop["step6_cavity_analysis"])

        # STEP 7: Filtering cavities
        global_log.info("step7_filter_cavities: Filter found cavities")
        fpocket_filter(**cluster_paths['step7_filter_cavities'], properties=cluster_prop["step7_filter_cavities"])
        
        # STEP 8: Filter by pocket center of mass
        global_log.info("step8_filter_residue_com: Filter cavities by center of mass distance to a group of residues") 
        filtered_pockets_IDs = filter_residue_com(**cluster_paths['step8_filter_residue_com'], properties=cluster_prop["step8_filter_residue_com"], global_log=global_log)

        # Update dictionary with filtered pockets
        cluster_filtered_pockets.update({cluster_name : filtered_pockets_IDs})

        # Save model pdb file in sub folder
        model_subfolder = os.path.join(output_path, cluster_name)
        shutil.copyfile(cluster_paths['step6_cavity_analysis']['input_pdb_path'], os.path.join(model_subfolder, 'model.pdb'))

    # Create summary with available pockets per cluster 
    global_log.info("    Creating YAML summary file...")
    create_summary(cluster_names, cluster_populations, cluster_filtered_pockets, global_paths, output_path)
        
    # Print timing information to log file
    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow name: Cavity analysis')
    global_log.info('  Output path: %s' % output_path)
    global_log.info('  Config File: %s' % configuration_path)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')

    return global_paths, global_prop

def main():

    parser = argparse.ArgumentParser(description="Simple clustering and cavity analysis pipeline using BioExcel Building Blocks")

    parser.add_argument('--traj_path', dest='traj_path',
                        help="Path to input trajectory (GROMACS or AMBER formats)", 
                        required=False)

    parser.add_argument('--top_path', dest='top_path',
                        help="Path to input structure (gro, pdb)",
                        required=False) 

    parser.add_argument('--structures_path', dest='structures_path',
                        help="Path to folder with structures in PDB format.", 
                        required=False)

    parser.add_argument('--filtering_selection', dest='filtering_selection',
                        help="Atom selection to filter pockets by distance to center of mass (MDAnalysis syntax)",
                        required=False)
    
    parser.add_argument('--distance_threshold', dest='distance_threshold', type=float,
                        help="Distance threshold to filter pockets by distance to center of mass",
                        required=False)
    
    parser.add_argument('--num_clusters', dest='num_clusters', type=int, default=20,
                        help="""Number of most populated clusters to extract from the trajectory and analyze with fpocket (default: 20).
                        if representative structures are given instead of a traj, this number is ignored""",
                        required=False)

    parser.add_argument('--clustering_method', dest='clustering_method', type=str, default='linkage',
                        help="Clustering method to use (linkage, jarvis-patrick, monte-carlo, diagonalization, gromos). Default: linkage",
                        required=False)
    
    parser.add_argument('--clustering_cutoff', dest='clustering_cutoff', type=float, default=0.1,
                        help="Clustering cutoff to use for the clustering method. Reduce to increase the number of clusters. Default: 0.1",
                        required=False)
    
    parser.add_argument('--gmx_bin', type=str,
                        help="Path to GROMACS binary (gmx for single node and gmx_mpi for multi-node). Default: gmx",
                        required=False, default='gmx')
    
    parser.add_argument('--restart', action='store_true',
                        help="Restart the workflow from the last completed step. Default: False",
                        required=False, default=False)

    parser.add_argument('--output', dest='output_path',
                        help="Output path (default: working_dir_path in YAML config file)",
                        required=False)
    
    args = parser.parse_args()

    cavity_analysis(traj_path = args.traj_path,
                    top_path = args.top_path,
                    structures_path = args.structures_path,
                    distance_threshold = args.distance_threshold,
                    filtering_selection = args.filtering_selection,
                    num_clusters = args.num_clusters,
                    clustering_method = args.clustering_method,
                    clustering_cutoff = args.clustering_cutoff,
                    gmx_bin = args.gmx_bin,
                    restart = args.restart,
                    output_path = args.output_path)


if __name__ == '__main__':
    main()