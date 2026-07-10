# Virtual screening

High-throughput virtual screening: dock a ligand library (SMILES or SDF) against a receptor pocket with AutoDock Vina, then rank by binding affinity.

## Description

The input is a ligand library (SMILES `.smi` or SDF `.sdf`), a target structure (PDB) and a pocket. The pocket is given either as an Fpocket result (see `cavity_analysis`) or as a residue selection. The workflow builds a docking box around the pocket, prepares the receptor, then docks every ligand and ranks them by binding affinity.

AutoDock Vina predicts a binding affinity for each ligand (in kcal/mol). The workflow keeps the best affinity per ligand and writes the ranking to `scores.csv`. Docking poses of the top ligands can optionally be saved.

- **Step 1**: Selection of the cavity used to build the docking box — a pocket from an input zip (`--input_pockets_zip`, see `cavity_analysis`) or a residue selection (`--pocket_selection`).
- **Step 2**: Creation of the box surrounding the selected cavity or residues. `--box_offset` sets the padding between the outermost atom and the box edge.
- **Step 3**: Addition of H atoms and partial charges to the receptor (`.pdb` → `.pdbqt`). Vina ignores the partial charges, but correct receptor protonation still matters because it decides which atoms are H-bond donors/acceptors.
- **Step 4**: Ligand preparation with OpenBabel. If the library is SMILES, ligands are protonated and given a 3D conformer (`.smi` → `.pdbqt`, generated at pH 7.4). If the library is SDF, ligands are only converted, keeping the input protonation and conformer (`.sdf` → `.pdbqt`).
- **Step 5**: Docking with AutoDock Vina (rigid receptor, flexible ligand).
- **Step 6**: Save poses of the top-scoring ligands (only with `--keep_poses`).

Ligands are docked one at a time. A ligand that fails to convert or dock is skipped and left out of the ranking (a success rate is reported in the log).

## Usage

```bash
conda activate biobb_vs
vs_autodock --ligand_lib data/ligands/zinc_200_425_001_reduced.sdf \
  --structure_path data/receptor/receptor.pdb \
  --pocket_selection "resid 37 or resid 49 or resid 112" \
  --box_offset 5 --cpus 4 --exhaustiveness 8
```

The `config.yml` is auto-generated from the CLI arguments into `--output`. `--restart` resumes from the last completed step when re-run against the same output folder. Run `vs_autodock --help` for the full option list.

## Options


### Inputs

Define the pocket with either `--input_pockets_zip` or `--pocket_selection` (mutually exclusive).

| Flag | Default | Description |
|------|---------|-------------|
| `--ligand_lib` | *required* | Ligand library: SMILES (`.smi`, one `smiles name` per line) or SDF (`.sdf`, one or more ligands). |
| `--structure_path` | *required* | Target structure (PDB); remove unneeded ligands/ions/cofactors first (H added at pH 7). |
| `--input_pockets_zip` | `None` | Fpocket pockets zip file. |
| `--pocket_selection` | `None` | Residue selection (MDAnalysis syntax, e.g. `resid 37 49 112`) defining the pocket. |

### Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--pocket_num` | `1` | Pocket number to use from `--input_pockets_zip`. |
| `--box_offset` | `5.0` | Extra distance (Å) between the outermost residue atom and the box boundary. |
| `--num_top_ligands` | all | Number of top ligands to save in the ranking. |
| `--keep_poses` | `False` | Save docking poses for the top ligands. |
| `--vina_bin` | `vina` | AutoDock Vina binary. |
| `--cpus` | `1` | CPUs per docking. |
| `--exhaustiveness` | `8` | Vina sampling runs (`4` faster, `8` more accurate). |
| `--debug` / `-d` | `False` | Keep intermediate files for debugging. |
| `--restart` | `False` | Restart from the last completed step. |
| `--output` | `working_dir_path` | Output directory. |

## Recommendations

- **Prefer prepared SDF ligands over SMILES.** For SMILES, OpenBabel (`obabel`) perceives bonds from the generated 3D coordinates and protonates for **pH 7.4** using tabulated per-group pKa rules. This is heuristic. If you already have well-prepared 3D, protonated ligands, pass them as SDF so they are docked as-is.
- **Tune exhaustiveness to the library size.** It trades accuracy for speed. For large libraries, start with a low value to screen fast, then re-dock the best-scoring ligands with a higher value.
- **Clean the receptor first.** Remove ligands, ions, and cofactors you do not need. Hydrogens are added automatically at pH 7. The receptor is treated as rigid.
- **Keep the box small.** A smaller box makes the search easier and faster; Vina cannot place the ligand outside the box. `--box_offset` adds padding around the pocket residues (default 5 Å); a warning is printed above 5 Å.
- **Validate before screening.** Dock a known binder or the native ligand first and check the pose before running the full library.
- **Use the scores to rank, not to measure.** Vina affinities are approximate. Docking is non-deterministic, so scores and poses change slightly between runs.

## Output

Unless `--debug`, per-ligand subfolders are deleted after scoring. Surviving outputs:

- `scores.csv` — ranking of successfully docked ligands. Columns: `Rank, Affinity, Index, Identifier`. Ranked by affinity (most negative first). Limited to `--num_top_ligands` when given.
- `receptor.pdb` — copy of the input receptor.
- `ligand_library.txt` — absolute path to the ligand library used.
- `poses/` (only with `--keep_poses`) — one PDB per top ligand (`<name>_poses.pdb`) with its docking poses.

With `--debug`, every per-ligand subfolder is kept with all intermediate files (prepared ligand, box, docking output).

## Limitations

- **No parallelization between ligands.** Ligands are docked one after another; `--cpus` only parallelizes a single docking.
- **Rigid receptor.** Only the ligand is flexible. Side chains cannot move and flexible side chains are not supported.
- **Approximate scoring.** Predicted affinity is not the experimental binding energy. Accuracy varies by target; evaluate against known actives.
- **Fixed protonation.** Receptor H are added at pH 7 (auto mode) and SMILES ligands are protonated at pH 7.4. Neither is configurable from the command line, and no tautomer/stereoisomer enumeration is done.
