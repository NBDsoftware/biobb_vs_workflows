# Virtual screening

High-throughput virtual screening: dock a ligand library (SMILES or SDF) against a receptor pocket with AutoDock Vina, then rank by binding affinity.

## Description

<!-- TODO: refine content -->

The input is a ligand library (SMILES `.smi` or SDF `.sdf`), a target structure (PDB) and either a pocket from an Fpocket analysis or a residue selection. The workflow docks each ligand to the target and ranks them by affinity, writing the top scorers to `scores.csv` (rank, name, identifier, affinity). Poses of the top ligands can optionally be kept.

**Pipeline steps:**

- **Step 1**: Selection of the cavity used to build the docking box — a pocket from an input zip (`--input_pockets_zip`, see `cavity_analysis`) or a residue selection (`--pocket_selection`).
- **Step 2**: Creation of the box surrounding the selected cavity or residues.
- **Step 3**: Addition of H atoms to the receptor (`.pdb` → `.pdbqt`). Vina ignores the pdbqt charges, but correct receptor protonation matters.
- **Step 4**: Ligand preparation with OpenBabel. `.smi` → protonated 3D conformer `.pdbqt` (`--gen3d -p <ph>`); `.sdf` → `.pdbqt` keeping the input protonation/conformer.
- **Step 5**: Docking with AutoDock Vina (rigid receptor, partial ligand flexibility).
- **Step 6**: Save poses of the top-scoring ligands (only with `--keep_poses`).

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

Command-line arguments take priority over the auto-generated `config.yml`.

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

<!-- TODO: refine content -->

## Output

<!-- TODO: refine content -->

Unless `--debug`, per-ligand subfolders are deleted after scoring. Surviving outputs: `scores.csv`, `receptor.pdb`, `ligand_library.txt`, and (with `--keep_poses`) a `poses/` folder.

## Limitations

<!-- TODO: refine content -->
