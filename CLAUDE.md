# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# AI Assistant Context & Refactoring Guidelines

In all interactions and commit messages be extremely concise. Sacrifice grammar for the sake of concision.

## Project Overview
BioBB Workflows is a repository of bioinformatics pipelines built on top of the BioBB (BioExcel Building Blocks) library. This repo covers virtual-screening workflows specifically. MD/preparation workflows live in a sibling repo `biobb_md_workflows`.

Two workflows, each a CLI command (installed via `pyproject.toml` `[project.scripts]`):
- `cavity_analysis` — clusters an MD trajectory (or accepts pre-computed representative structures) and runs Fpocket cavity detection + filtering on each representative structure.
- `vs_autodock` — high-throughput virtual screening: docks a ligand library (SMILES or SDF) against a receptor pocket using AutoDock Vina, then ranks by binding affinity.

## Commands

Install:
```bash
conda env create -f environment.yml   # installs env + pip-installs this package (`pip: - .`)
conda activate biobb_vs
```

Run a workflow:
```bash
cavity_analysis --structures_path data/receptor/receptor.pdb --filtering_selection "resid 31 or resid 21"
vs_autodock --ligand_lib data/ligands/zinc_200_425_001_reduced.sdf \
--structure_path data/receptor/receptor.pdb \
--pocket_selection "resid 37 or resid 49 or resid 112" \ 
--box_offset 5 \
--cpus 4 \
--exhaustiveness 8 \
```
`--help` on either command lists all flags.

There is no automated test suite (no pytest). `tests/{cavity_analysis,vs_autodock}/run.sl` are SLURM sbatch scripts that exercise each workflow end-to-end against fixtures in `data/`. To validate a change, run the workflow directly (or via `run.sl` after fixing the conda env path inside it) and inspect `output/log.out` and the generated summary/ranking files.

## Architecture

Each workflow lives in `biobb_vs_workflows/<workflow>/<workflow>.py` and follows the same shape:

1. **`main()`** parses argparse CLI flags and calls the workflow function of the same name (also exported from the package `__init__.py`).
2. **`config_contents(...)`** builds a YAML config *as an f-string* (not a template file) embedding the resolved CLI args, one section per BioBB step (`stepN_toolname`).
3. **`create_config_file()`** writes that YAML to `output/config.yml`.
4. `settings.ConfReader` (biobb_common) reads it back and splits it into `global_paths` / `global_prop` dicts keyed by step name (e.g. `global_paths["step2_box"]`).
5. Each BioBB step is invoked as `tool_func(**global_paths[step_name], properties=global_prop[step_name])`. Steps chain via `dependency/stepX/output_key` path strings that `ConfReader` resolves to the actual prior-step output path.
6. `conf.get_prop_dic(prefix=name)` / `get_paths_dic(prefix=name)` re-derive the same step templates under a per-item subfolder (`output/<name>/stepN_.../...`) — this is how both workflows fan out over collections (one subfolder per ligand, or per cluster/model) while reusing a single YAML template.
7. Global properties `restart` / `remove_tmp` (set in the YAML) make steps idempotent/skippable on re-run — the CLI's own `--restart` flag controls this.

**`vs_autodock.py`**: branches on ligand library extension. `.sdf` path iterates ligands via `openbabel.pybel`; `.smi` path parses lines with `read_ligand_lib()`. Each ligand gets its own `output/<ligand_name>/` subfolder (protonation/conversion via `babel_convert` → dock via `autodock_vina_run`, wrapped in try/except so one failing ligand doesn't abort the run). Ranking is derived by regex-parsing `REMARK VINA RESULT` lines out of the output pdbqt (`get_affinity`/`get_ranking`), written to `scores.csv`. Unless `--debug`, per-ligand subfolders are deleted after scoring (`clean_output`) — only `scores.csv`, `receptor.pdb`, `ligand_library.txt`, and (if `--keep_poses`) a `poses/` folder survive.

**`cavity_analysis.py`**: if `--structures_path` is not given, clusters `--traj_path`/`--top_path` with `gmx_cluster` (AMBER trajectories are first converted via `cpptraj_convert`); otherwise treats each PDB under `structures_path` as an already-clustered representative model. Has a retry loop (`MAX_CLUSTER_RETRIES`, `CUTOFF_INCREASE_FACTOR`) that increases the clustering cutoff and re-clusters if the combined centroid PDB exceeds `MAX_CLUSTER_ATOMS`, since `extract_model` (step5) chokes/slows on very large atom counts. For every model: `fpocket_run` → `fpocket_filter` (score/druggability/volume thresholds) → `filter_residue_com` (repo-local filter, not a BioBB tool — uses MDAnalysis to keep only pockets whose center of mass is within `--distance_threshold` of `--filtering_selection`). Produces `summary_by_volume.yml`, `summary_by_drug_score.yml`, `summary_by_score.yml` sorted per-model rankings.

**`pymol_cavity_analysis.py`**: standalone PyMOL visualization helper for `cavity_analysis` output (loads models + highlights filtered pockets/residues). Not registered as a console script — run directly inside a PyMOL Python environment.

**`vs_autodock_mp_legacy.py`**: older multiprocessing-based version of the docking workflow (parallelizes across ligands, sources ligands from DrugBank/`biobb_io` instead of a user-supplied library). Not wired into `pyproject.toml` scripts — kept for reference only, not the maintained entrypoint.

## Plans

- At the end of each plan, give me a list of unresolved questions to answer, if any. Make the questions extremely concise. Sacrifice grammar for the sake of concision.
