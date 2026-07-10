# Cavity analysis

Cluster an MD trajectory (or accept pre-computed representative structures) and run Fpocket cavity detection + filtering on each representative structure.

## Description

The input can be either (1) a trajectory and topology or (2) a folder of representative structures from an external clustering. With a trajectory, the workflow clusters it to find representative structures. With structures, it uses them directly.

The workflow then runs Fpocket on each representative structure. Fpocket detects pockets as clusters of alpha spheres. It reports a score, a druggability score, and a volume for every pocket.

- **Step 1**: Conversion of the trajectory from Amber-compatible formats to GROMACS `xtc` (automatic when the input trajectory is not GROMACS-compatible).
- **Step 2**: Creation of the atom-selection index file used for clustering.
- **Step 3**: Clustering of the trajectory (`gmx_cluster`).
- **Step 4**: Extraction of the most populated centroids.
- **Step 5**: Cavity analysis of each centroid with Fpocket (`fpocket_run`).
- **Step 6**: Filtering of cavities by score, druggability score and volume (`fpocket_filter`).
- **Step 7**: Filtering of cavities by distance of their center of mass to a residue selection (`filter_residue_com`).

## Usage

```bash
conda activate biobb_vs
cavity_analysis --structures_path data/receptor/receptor.pdb --filtering_selection "resid 31 or resid 21"
```

The `config.yml` is auto-generated from the CLI arguments into `--output`. `--restart` resumes from the last completed step when re-run against the same output folder. Run `cavity_analysis --help` for the full option list.

## Options

### Inputs

Provide either a trajectory (`--traj_path` + `--top_path`) or a folder of representative structures (`--structures_path`).

| Flag | Default | Description |
|------|---------|-------------|
| `--traj_path` | `None` | Input trajectory (GROMACS or AMBER formats); use with `--top_path`. |
| `--top_path` | `None` | Input topology/structure (`gro`, `pdb`); use with `--traj_path`. |
| `--structures_path` | `None` | Folder of representative structures (PDB); skips clustering. |

### Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--filtering_selection` | `None` | Atom selection (MDAnalysis syntax) used to filter pockets by the distance of their center of mass. |
| `--distance_threshold` | `None` | Distance threshold (Å) for the center-of-mass filter. |
| `--num_clusters` | `20` | Most-populated clusters to extract and analyze; ignored when `--structures_path` is given. |
| `--clustering_method` | `linkage` | `gmx cluster` method (`linkage`, `jarvis-patrick`, `monte-carlo`, `diagonalization`, `gromos`). |
| `--clustering_cutoff` | `0.1` | Clustering cutoff; reduce to increase the number of clusters. |
| `--gmx_bin` | `gmx` | GROMACS binary (`gmx` single-node, `gmx_mpi` multi-node). |
| `--restart` | `False` | Restart from the last completed step. |
| `--output` | `working_dir_path` | Output directory. |

## Recommendations

- **Clean the input first.** Remove ligands, ions, and waters from the trajectory. Clustering and Fpocket work on the protein. Extra molecules can cause errors or noisy pockets.
- **Clustering uses the whole protein.** To cluster on a specific region instead, do the clustering yourself and pass the representative structures with `--structures_path`.
- **Number of clusters.** Lower `--clustering_cutoff` to get more clusters, raise it to get fewer (default `0.1`). If too many clusters are obtained the workflow will re-launch the clustering increasing the cut-off value 50% up to 10 times. `--num_clusters` sets how many of the most populated clusters are analyzed (default `20`); it is ignored when you pass `--structures_path`.
- **Distance filter.** `--filtering_selection` and `--distance_threshold` keep only pockets whose center of mass is close to the chosen residues. A threshold of 10 Å is reasonable when the residues are part of the pocket.

## Output

The top-level output folder contains three YAML summaries of the results, one per ranking:

- `summary_by_drug_score.yml`
- `summary_by_score.yml`
- `summary_by_volume.yml`

Each summary lists the models (clusters). For every model it gives its pockets with their score, druggability score, and volume. Models are ranked by their best pocket. When the workflow did the clustering, each model also shows its population (how many trajectory frames it represents).

The druggability score goes from 0 to 1. A higher value means the pocket is more likely to bind a drug-like molecule. Volume is in Å³.

Each model also has its own subfolder with the structure (`model.pdb`) and the detected pockets.

## Limitations

- **Clustering uses the whole protein.** The fit and cluster selection is fixed to `Protein`. You cannot restrict it to a region from the command line.
- **Fpocket detection settings are fixed:** alpha-sphere radius 3–6 Å, and at least 35 spheres per pocket.
- **Filter thresholds are fixed:** score 0.4–1, druggability score 0.4–1, volume 200–5000 Å³. Pockets outside these ranges are dropped. These cannot be changed from the command line.
