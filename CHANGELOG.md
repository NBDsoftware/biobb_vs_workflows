# Changelog

## Release Version 0.0.2 ##

This release adds full docs + an end-to-end tutorial, makes both workflows extract a clean
protein from input structures by default, exposes fpocket tuning params, and fixes a couple
of install/config bugs.

### ⚠️ Breaking changes
- **vs_autodock, cavity_analysis** (`--structures_path`): input structures are now
  protein-extracted by default (water/ligands/ions/cofactors stripped) before pocket
  detection/docking. Use `--skip_extraction` to keep the old behavior.
- Dependency cleanup: dropped `biobb_pdb_tools` and `propka`, unpinned `openbabel`, removed
  the `biopython<=1.81` cap, and constrained Python to `>=3.11,<3.12` (vina has no 3.12
  build yet). Recreate your environment from `environment.yml`/`pyproject.toml`.

### 🐛 Bug Fixes
- **cavity_analysis**: the `filter_residue_com` step was silently skipped even when
  `--filtering_selection` was given; it now runs automatically whenever a selection is
  provided.
- Fixed a missing-comma/dependency typo in `pyproject.toml` that broke installation.

### 🚀 Features & Improvements
- **cavity_analysis**: new `--min_radius`/`--max_radius`/`--num_spheres` flags exposing
  fpocket's `-m`/`-M`/`-i` alpha-sphere tuning to control pocket size/shape filtering.
- **cavity_analysis**: new `--debug` flag to keep temporary files.
- **vs_autodock**: `--debug` now also disables cleanup of intermediate BioBB step files
  (previously only affected per-ligand subfolder cleanup).
- **All workflows**: shared `biobb_vs_workflows.common.to_yaml()` helper for consistent
  YAML rendering (`None`/bool/list) across both config templates.
- New Sphinx docs site (published via GitHub Pages) and a full end-to-end tutorial notebook
  (ligand preparation with gypsum-dl → cavity analysis → docking), runnable on Colab or a
  local Jupyter, plus a standalone ligand-preparation notebook and overhauled READMEs.
- Added `data/` test fixtures (receptor, complexes, ligand libraries, cluster PDBs) and
  `tests/{cavity_analysis,vs_autodock}/run.sl` SLURM scripts exercising both workflows
  end-to-end.
