# Virtual screening

High-throughput virtual screening: dock a ligand library (SMILES or SDF) against a receptor pocket with AutoDock Vina or gnina, then rank by score.

## Description

The input is a ligand library (SMILES `.smi` or SDF `.sdf`), a target structure (PDB) and a pocket. The pocket is given either as an Fpocket result (see `cavity_analysis`) or as a residue selection. The workflow builds a docking box around the pocket, prepares the receptor, then docks every ligand and ranks them.

Two docking engines are available through `--docking_engine`: **`vina`** (default, AutoDock Vina) and **`gnina`** (Vina + CNN rescoring, needs a gnina binary and a `biobb_vs` built from source — see **Installing gnina** below). See the [docking engines reference](../../docs/docking_engines/index.md) for how each one scores poses and how to choose between them.

The workflow keeps the best pose per ligand and writes the ranking to `scores.csv`. Docking poses of the top ligands can optionally be saved.

- **Step 0**: Extraction of the protein from the receptor structure (drops waters/ligands/ions) with `extract_molecule`. Runs by default; disable with `--skip_extraction` to keep a bound cofactor/ligand/ion.
- **Step 1**: Selection of the cavity used to build the docking box — a pocket from an input zip (`--input_pockets_zip`, see `cavity_analysis`) or a residue selection (`--pocket_selection`).
- **Step 2**: Creation of the box surrounding the selected cavity or residues. `--box_offset` sets the padding between the outermost atom and the box edge.
- **Step 3**: Addition of H atoms and partial charges to the receptor (`.pdb` → `.pdbqt`). Neither engine uses the partial charges, but correct receptor protonation still matters because it decides which atoms are H-bond donors/acceptors. gnina takes this `.pdbqt` as given instead of re-perceiving the receptor itself.
- **Step 4**: Ligand preparation with OpenBabel. If the library is SMILES, ligands are protonated and given a 3D conformer at pH 7.4. If the library is SDF, the input protonation and conformer are kept as they are. Ligands are converted to `.pdbqt` for vina; gnina reads the `.sdf` directly, so no conversion happens.
- **Step 5**: Docking (rigid receptor, flexible ligand) with the selected engine.
- **Step 6**: Save poses of the top-scoring ligands (only with `--keep_poses`).

Ligands are docked one at a time. A ligand that fails to prepare or dock is skipped and left out of the ranking (a success rate is reported in the log).

## Usage

```bash
conda activate biobb_vs
virtual_screening --ligand_lib data/ligands/zinc_200_425_001_reduced.sdf \
  --structure_path data/receptor/receptor.pdb \
  --pocket_selection "resid 37 or resid 49 or resid 112" \
  --box_offset 5 --cpus 4 --exhaustiveness 8
```

With gnina, adding CNN rescoring and a fixed seed:

```bash
virtual_screening --ligand_lib data/ligands/zinc_200_425_001_reduced.sdf \
  --structure_path data/receptor/receptor.pdb \
  --pocket_selection "resid 37 or resid 49 or resid 112" \
  --docking_engine gnina --gnina_cnn fast --gnina_seed 42 \
  --box_offset 5 --cpus 4 --exhaustiveness 8
```

The `config.yml` is auto-generated from the CLI arguments into `--output`. `--restart` resumes from the last completed step when re-run against the same output folder. Run `virtual_screening --help` for the full option list.

## Installing gnina

gnina is **not** conda-installable and is not pulled in by `environment.yml`. The `gnina` engine needs two things the default install does not provide:

1. A **gnina binary** — download a [release binary](https://github.com/gnina/gnina/releases) (`chmod +x` it; release assets arrive mode `644`) or use the [Docker image](https://hub.docker.com/u/gnina), then point `--gnina_bin` at it.
2. A **`biobb_vs` built from source** — the conda-forge package ships no gnina module. Without it the workflow exits with a clear error, and the `vina` engine keeps working.

A GPU is strongly recommended. gnina is Linux + NVIDIA only; `--gnina_no_gpu` forces CPU execution, but note that on a machine with no CUDA libraries at all the release binary will not even load, since it is linked against CUDA/cuDNN — that is a separate problem from `--gnina_no_gpu`.

## Options


### Inputs

Define the pocket with either `--input_pockets_zip` or `--pocket_selection` (mutually exclusive).

| Flag | Default | Description |
|------|---------|-------------|
| `--ligand_lib` | *required* | Ligand library: SMILES (`.smi`, one `smiles name` per line) or SDF (`.sdf`, one or more ligands). |
| `--structure_path` | *required* | Target structure (PDB). Protein is auto-extracted (waters/ligands/ions dropped) unless `--skip_extraction`; H added at pH 7. |
| `--input_pockets_zip` | `None` | Fpocket pockets zip file. |
| `--pocket_selection` | `None` | Residue selection (MDAnalysis syntax, e.g. `resid 37 49 112`) defining the pocket. |

### Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--pocket_num` | `1` | Pocket number to use from `--input_pockets_zip`. |
| `--box_offset` | `5.0` | Extra distance (Å) between the outermost residue atom and the box boundary. |
| `--num_top_ligands` | all | Number of top ligands to save in the ranking. |
| `--keep_poses` | `False` | Save docking poses for the top ligands. |
| `--skip_extraction` | `False` | Skip protein extraction from the receptor (keep cofactors/ligands/ions). |
| `--cpus` | `1` | CPUs per docking. Keep it ≤ `--exhaustiveness`. |
| `--exhaustiveness` | `8` | Sampling runs (`4` faster, `8` more accurate). |
| `--debug` / `-d` | `False` | Keep intermediate files for debugging. |
| `--restart` | `False` | Restart from the last completed step. |
| `--output` | `working_dir_path` | Output directory. |

### Docking engine

| Flag | Default | Description |
|------|---------|-------------|
| `--docking_engine` | `vina` | Docking engine: `vina` or `gnina`. |
| `--vina_bin` | `vina` | AutoDock Vina binary. |
| `--gnina_bin` | `gnina` | gnina binary (see the **Installing gnina** section). |
| `--gnina_cnn_scoring` | gnina's own (`rescore`) | Where gnina uses the CNN: `none`, `rescore`, `refinement`, `metrorescore`, `metrorefine`, `all`. The dominant cost knob. |
| `--gnina_cnn` | gnina's own (3-model ensemble) | CNN model, or a `PREFIX_ensemble` name. `fast` is a single model, ~4× faster. |
| `--gnina_scoring` | gnina's own | Empirical scoring function: `vina`, `vinardo`, `ad4_scoring`, `dkoes_*`. |
| `--gnina_rank_by` | `CNNaffinity` | Score to rank ligands by: `CNNaffinity`, `CNNscore`, `minimizedAffinity` or `CNN_VS` (`CNNaffinity * CNNscore`, gnina's own screening pipeline default). Falls back to `minimizedAffinity` with `--gnina_cnn_scoring none`. |
| `--gnina_pose_sort_order` | matches `--gnina_rank_by` | gnina property (`CNNscore`, `CNNaffinity` or `Energy`) that selects which poses gnina keeps, before this workflow ranks them. Defaults to whatever matches `--gnina_rank_by`, so the poses ranked over are the ones gnina itself would keep under that score. |
| `--gnina_seed` | none | Random seed. Docking is stochastic, set it for reproducible runs. |
| `--gnina_no_gpu` | `False` | Force gnina onto the CPU even when a GPU is available. |

**gnina cost.** The CNN settings dominate runtime far more than `--exhaustiveness` does — see the [gnina reference](../../docs/docking_engines/gnina.md) for relative costs and GPU vs CPU numbers.

## Recommendations

- **Prefer prepared SDF ligands over SMILES.** For SMILES, OpenBabel (`obabel`) perceives bonds from the generated 3D coordinates and protonates for **pH 7.4** using tabulated per-group pKa rules. This is heuristic. If you already have well-prepared 3D, protonated ligands, pass them as SDF so they are docked as-is — this matters more with gnina (see its [limitations](../../docs/docking_engines/gnina.md)). Neither engine samples ring conformations or stereochemistry — whatever you provide is what gets docked.
- **With gnina, rank by `CNNaffinity`, not `CNNscore`** (the default does this) — see the [gnina reference](../../docs/docking_engines/gnina.md) for why.
- **Leave `--gnina_pose_sort_order` alone unless reproducing a benchmark.** It defaults to match `--gnina_rank_by`, so gnina keeps the poses that score best under the same metric this workflow ranks by — see the [gnina reference](../../docs/docking_engines/gnina.md) for why a mismatch quietly changes *which* poses are ranked, not just their order.
- **Try more than one scoring function.** `--gnina_scoring vinardo` often does better than the default in virtual screening. Ranking the same library with two engines/functions and combining by *rank* (not by raw score, whose scales are not comparable) is a stronger signal than any single method's top hits.
- **Set `--gnina_seed`.** Docking is a stochastic Monte Carlo search, so an unseeded run is not reproducible.
- **Tune exhaustiveness.** It trades accuracy for speed. For large libraries, start with a low value to screen fast, then re-dock the best-scoring ligands with a higher value. **Match the number of cpu cores used with the exhaustiveness value** as each Monte Carlo chain will run on a different core (e.g. optimal run would be 4 cpu cores with exhaustiveness value 4, suboptimal run would be running on 4 cpu cores with exhaustiveness 6).
- **Receptor cleaning is automatic.** By default only the protein is kept (waters, ligands, and ions are stripped with `extract_molecule`). Pass `--skip_extraction` to keep a bound cofactor/ligand/ion you need. Hydrogens are added automatically at pH 7. The receptor is treated as rigid.
- **Keep the box small.** A smaller box makes the search easier and faster; neither engine can place the ligand outside the box. `--box_offset` adds padding around the pocket residues (default 5 Å); a warning is printed above 5 Å.
- **Validate before screening.** Dock a known binder or the native ligand first and check the pose before running the full library.
- **Use the scores to rank, not to measure.** These are docking scores, not binding free energies. Docking is non-deterministic, so scores and poses change slightly between runs unless a seed is set.

## Output

Unless `--debug`, per-ligand subfolders are deleted after scoring. Surviving outputs:

- `scores.csv` — ranking of successfully docked ligands, limited to `--num_top_ligands` when given. The columns depend on the engine:

  | Engine | Columns |
  |--------|---------|
  | `vina` | `Rank, Affinity, Index, Identifier` — ranked by affinity, most negative first |
  | `gnina` | `Rank, minimizedAffinity, CNNaffinity, CNNscore, CNN_VS, Index, Identifier` — ranked by `--gnina_rank_by` |
  | `gnina --gnina_cnn_scoring none` | `Rank, minimizedAffinity, Index, Identifier` |

  `Affinity` (vina) and `minimizedAffinity` (gnina) are the same quantity in the same units (kcal/mol, lower is better), so the two engines stay comparable. Note that with gnina the rows are ordered by `--gnina_rank_by`, which defaults to `CNNaffinity` — so **`minimizedAffinity` is not monotonic with `Rank`**, and a ligand with a poor empirical affinity can rank highly if the CNN likes it.

- `receptor.pdb` — copy of the input receptor.
- `ligand_library.txt` — absolute path to the ligand library used.
- `poses/` (only with `--keep_poses`) — one file per top ligand with its docking poses: `<name>_poses.pdb` for vina, `<name>_poses.sdf` for gnina. The gnina SDF keeps every score as an SD data field, which a conversion to PDB would discard.

With `--debug`, every per-ligand subfolder is kept with all intermediate files (prepared ligand, box, docking output).

## Limitations

- **No parallelization between ligands.** Ligands are docked one after another; `--cpus` only parallelizes a single docking. This costs more with gnina, which reloads its CNN models on every invocation — batching the library into a single gnina call would amortize that, and is not implemented yet.
- **Rigid receptor.** Only the ligand is flexible. Side chains cannot move and flexible side-chain docking is not exposed, even though gnina supports it.
- **Approximate scoring.** Predicted affinity is not the experimental binding energy. Accuracy varies by target; evaluate against known actives.
- **Fixed protonation.** Receptor H are added at pH 7 (auto mode) and SMILES ligands are protonated at pH 7.4. Neither is configurable from the command line, and no tautomer/stereoisomer enumeration is done.
- **gnina has extra scoring limitations** (CNN grid size, `--gnina_cnn` affecting which poses survive, not just their order) — see the [gnina reference](../../docs/docking_engines/gnina.md).