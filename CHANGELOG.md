# Changelog

## Unreleased ##

### 📦 Dependencies
- Bumped the whole BioBB stack to 5.3.x (`biobb_vs` 5.3.0, `biobb_common` 5.3.1,
  `biobb_structure_utils` 5.3.0, `biobb_gromacs` 5.3.1, `biobb_io` 5.3.0,
  `biobb_chemistry` 5.3.0, `biobb_analysis` 5.3.0). BioBB recipes pin `biobb_common`
  exactly, so the set has to move together. Recreate your environment.
- **`biobb_analysis` is now installed via pip, not conda — do not move it back.** Its
  conda recipe requires `gromacs>=2026`, which pulls `libhwloc` → `libxml2>=2.14`, while
  `openbabel 3.1.1` (hard-pinned by `biobb_chemistry` 5.3.0, and used directly by
  `vs_autodock`) caps `libxml2<2.14` on every available build. The two are mutually
  exclusive, so an all-conda 5.3.x stack does not solve. The PyPI package declares only
  `biobb-common==5.3.1` and runs fine against gromacs 2025.x — the 5.2.1→5.3.0 changes in
  `gmx_cluster`/`cpptraj_convert` are container-path fixes, nothing needing GROMACS 2026.
- Added `ambertools` explicitly to all three env files; it used to arrive transitively via
  conda `biobb_analysis` and is still needed for the `cpptraj_convert` (AMBER trajectory)
  path in `cavity_analysis`.
- Corrected the stale `python=3.11` comment — vina 1.2.7 does ship py3.12+ builds now. The
  pin stays because 3.11 is the verified configuration, not because 3.12 is unavailable.

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
