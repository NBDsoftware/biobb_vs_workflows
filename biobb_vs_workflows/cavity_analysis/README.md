# Cavity analysis

Cluster an MD trajectory (or accept pre-computed representative structures) and run Fpocket cavity detection + filtering on each representative structure.

## Description

<!-- TODO: refine content -->

The input for the workflow can be either (1) a trajectory and topology or (2) a folder of representative structures from an external clustering. In the former case the workflow clusters the trajectory to find representative structures; in the latter it uses the given structures directly for the cavity analysis and filtering.

**Pipeline steps:**

- **Step 0**: Conversion of the trajectory from Amber-compatible formats to GROMACS `xtc` (automatic when the input trajectory is not GROMACS-compatible).
- **Steps 1–2**: Index-file creation and extraction of a selected atom subset from trajectory and topology (removes waters/ions/exotic atoms that break `gmx cluster`).
- **Step 3**: Index files defining the RMSD group (fit + RMSD calculation) and the output group (atoms written to the representative structures).
- **Step 4**: Clustering of the trajectory (`gmx_cluster`).
- **Step 5**: Extraction of the most populated centroids.
- **Step 6**: Cavity analysis of each centroid with Fpocket (`fpocket_run`).
- **Step 7**: Filtering of cavities by score, druggability score and volume (`fpocket_filter`).
- **Step 8**: Filtering of cavities by distance of their center of mass to a residue selection (`filter_residue_com`).

## Usage

```bash
conda activate biobb_vs
cavity_analysis --structures_path data/receptor/receptor.pdb --filtering_selection "resid 31 or resid 21"
```

The `config.yml` is auto-generated from the CLI arguments into `--output`. `--restart` resumes from the last completed step when re-run against the same output folder. Run `cavity_analysis --help` for the full option list.

## Options

Command-line arguments take priority over the auto-generated `config.yml`.

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

<!-- TODO: refine content -->

## Output

<!-- TODO: refine content -->

The output folder contains three YAML summaries of the results, ordered by drug score, score, and volume (`summary_by_drug_score.yml`, `summary_by_score.yml`, `summary_by_volume.yml`).

## Limitations

<!-- TODO: refine content -->
