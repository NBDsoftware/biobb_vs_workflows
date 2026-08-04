# Command-line reference

**Audience:** everyone.

This is the grouped reference with real defaults and the behaviours `--help` does not mention. The
authoritative option list is the binary itself:

```bash
gnina --help          # the documented options
gnina --help_hidden   # plus the internal/testing group
gnina --version
```

Definitions live in [gninasrc/main/main.cpp](../../gninasrc/main/main.cpp), using
Boost.Program_options. Any option can also be supplied from a file via `--config`.

## Reading the defaults column

Some options declare a default that `--help` shows as `(=x)`; others are initialized elsewhere in the
code and show no default in `--help` even though they have one. Both are listed below as **effective**
defaults — what you get when the flag is absent. Where they differ, it is noted.

Also note three options are `value<bool>`, **not switches** — they require an argument:

```bash
--addH 0              # correct
--addH                # error: missing argument
```

Those are `--addH`, `--stripH` and `--autobox_extend`. Everything else that looks like a flag is one.

And three options accept **multiple values**: `-l/--ligand`, `--cnn` and `--cnn_model`.

## Input

| Option | Default | Purpose |
|---|---|---|
| `-r, --receptor` | — | Rigid part of the receptor |
| `-l, --ligand` | — | Ligand(s). Repeatable, and each file may hold many ligands |
| `--flex` | — | Flexible side chains as a PDBQT file |
| `--flexres` | — | Flexible side chains as a comma-separated `chain:resid` list (insertion code optional) |
| `--flexdist_ligand` | — | Reference file for distance-based flex selection |
| `--flexdist` | −1 (off) | Flex every residue within this distance of `--flexdist_ligand` |
| `--flex_limit` | −1 (off) | **Hard** cap: refuse to dock if exceeded |
| `--flex_max` | −1 (off) | **Soft** cap: keep only the N closest residues |

`--flex_limit` and `--flex_max` are mutually exclusive, and both are silently ignored when `--flexres`
is given. See [flexible-docking.md](flexible-docking.md).

## Search space

Required, one way or another.

| Option | Default | Purpose |
|---|---|---|
| `--center_x/y/z`, `--size_x/y/z` | 0 | Explicit box definition |
| `--autobox_ligand` | — | Define the box from a reference structure's bounding coordinates |
| `--autobox_add` | **4** | Padding added on all six sides. No declared default, so `--help` shows none |
| `--autobox_extend` | **true** | Extend any box side shorter than the ligand's longest interatomic distance. Takes an argument |
| `--no_lig` | false | No ligand — for sampling or minimizing flexible residues only |

See [inputs-and-search-space.md](inputs-and-search-space.md), including why `--autobox_add` cuts both
ways.

## Covalent docking

| Option | Default | Purpose |
|---|---|---|
| `--covalent_rec_atom` | — | Receptor atom, as `chain:resnum:atom_name` or `x,y,z` |
| `--covalent_lig_atom_pattern` | — | SMARTS for the ligand atom forming the bond. **All matches are docked** |
| `--covalent_lig_atom_position` | — | Optional explicit `x,y,z` placement, instead of OpenBabel's `GetNewBondVector` |
| `--covalent_fix_lig_atom_position` | false | Fix that atom at the given position rather than only seeding from it |
| `--covalent_bond_order` | **1** | Bond order of the covalent bond |
| `--covalent_optimize_lig` | false | UFF-optimize the complex. Changes ligand bond lengths and angles |

Setting `--covalent_rec_atom` with CNN scoring enabled produces a warning recommending
`--cnn_scoring none`. Take it — see [covalent-docking.md](covalent-docking.md).

## Scoring and minimization

| Option | Default | Purpose |
|---|---|---|
| `--scoring` | `vina` | `ad4_scoring`, `default` (= `vina`), `dkoes_fast`, `dkoes_scoring`, `dkoes_scoring_old`, `vina`, `vinardo` |
| `--custom_scoring` | — | Custom scoring-function file |
| `--custom_atoms` | — | Custom atom-type parameter file (10 fields per type) |
| `--score_only` | false | Score the given pose. No search, no output file |
| `--minimize` | false | Energy-minimize the given pose |
| `--local_only` | false | Local search within the autobox. Usually you want `--minimize` |
| `--randomize_only` | false | Generate random, clash-avoiding poses |
| `--num_mc_steps` | 0 → heuristic | Fixed MC steps per chain. Heuristic is `105 × (50 + movable_atoms + 10 × ndof)` |
| `--max_mc_steps` | 0 (no cap) | Cap the heuristic rather than replacing it |
| `--num_mc_saved` | **50** | Top poses retained per MC chain. No declared default |
| `--temperature` | 0 → **1.2** | Metropolis acceptance temperature |
| `--minimize_iters` | 0 → `(25 + movable_atoms)/3` | Minimization iterations |
| `--accurate_line`, `--simple_ascent`, `--minimize_early_term`, `--minimize_single_full` | false | Minimizer tuning |
| `--approximation` | `linear` | `linear`, `spline`, `exact` — plus an undocumented `gpu`. Case-sensitive |
| `--factor` | 32 | Granularity of the potential approximation |
| `--force_cap` | 1000 | Max allowed force; gentler minimization of clashing structures |
| `--user_grid` | — | AutoDock map-based user grid data |
| `--user_grid_lambda` | −1.0 | Weight for the user grid |
| `--print_terms` | false | Dump all available terms |
| `--print_atom_types` | false | Dump all atom types |

**`--minimize` overrides several of these** to values suited to real convergence:
`--minimize_iters` → 10000, `--approximation` → `spline`, `--factor` → 10, `--force_cap` → 10.

`--local_only` and `--minimize` **silently skip** ligands whose bounding extent exceeds 100 Å.

## CNN scoring

| Option | Default | Purpose |
|---|---|---|
| `--cnn_scoring` | `rescore` | `none`, `rescore`, `refinement`, `metrorescore`, `metrorefine`, `all` |
| `--cnn` | 3-model ensemble | Built-in model name(s), or `PREFIX_ensemble`. Repeatable |
| `--cnn_model` | — | TorchScript `.pt` file(s). Repeatable, and combinable with `--cnn` |
| `--cnn_rotation` | **0** | Evaluate multiple pose rotations, max 24. Measured to make no difference |
| `--cnn_mix_emp_force` | false | Merge CNN and empirical forces during refinement |
| `--cnn_mix_emp_energy` | false | Merge CNN and empirical energies |
| `--cnn_empirical_weight` | 1.0 | Weight of the empirical part in those mixes |
| `--cnn_center_x/y/z` | from the box | CNN grid centre, if different from the docking box |
| `--cnn_verbose` | false | Verbose CNN output |

Undocumented `--cnn_scoring` aliases: `no` → `none`, `docking` → `all`, and anything starting with
`refine` or `min` → `refinement`.

Two pseudo-names for `--cnn`, valid only when given alone: **`fast`** (one distilled model) and
**`default1.0`** (the 1.0 five-model ensemble). See [scoring-cnn.md](scoring-cnn.md) for the full model
list and what the names mean.

**`--cnn_scoring none` silently forces `--pose_sort_order Energy`.**

## Output

| Option | Default | Purpose |
|---|---|---|
| `-o, --out` | — | Output file. Format from the extension; `.sdf.gz` recommended |
| `--out_flex` | — | Output file for flexible residues |
| `--full_flex_output` | false | Write the whole structure to `--out_flex`, not just the flexible residues |
| `--log` | — | Optional log file |
| `--atom_terms` | — | Write per-atom interaction terms to a file |
| `--atom_term_data` | false | Embed per-atom terms in the output SD data instead |
| `--pose_sort_order` | `CNNscore` | `CNNscore`, `CNNaffinity`, `Energy` — plus an undocumented `vina` alias for `Energy` |

Which score fields each format carries is in [output-and-ranking.md](output-and-ranking.md); the short
version is that SDF carries all of them and PDBQT does not.

## Misc

| Option | Default | Purpose |
|---|---|---|
| `--cpu` | auto-detect | Number of CPU cores. Set this explicitly |
| `--seed` | PID + time | Random seed. Set this for reproducibility |
| `--exhaustiveness` | **8** | Number of independent MC chains |
| `--num_modes` | **9** | Max binding modes to output |
| `--min_rmsd_filter` | **1.0** | RMSD threshold for the output diversity filter |
| `-q, --quiet` | false | Suppress output messages |
| `--addH` | **true** | Add ligand hydrogens. Takes an argument |
| `--stripH` | **false** | Strip ligand hydrogens. Takes an argument |
| `--device` | 0 | GPU device index — **ignored by the PyTorch backend**; use `CUDA_VISIBLE_DEVICES` |
| `--no_gpu` | false | Force CPU-only execution |
| `--config` | — | Read any of the above from a file |
| `--help`, `--help_hidden`, `--version` | — | Usage, usage with hidden options, version |

The merged pose pool is `max(--num_modes, --num_mc_saved)`, so raising `--num_modes` past 50 enlarges
the internal pool too.

## Hidden options

Shown only by `--help_hidden`, described in the source as "Hidden options for internal testing". Not
for routine use, but useful to know exist.

| Option | Default | Purpose |
|---|---|---|
| `--verbosity` | 1 | Output verbosity. Raise it to see the weighted term list in use |
| `--flex_hydrogens` | false | Enable torsions affecting only hydrogens (e.g. OH). The help text calls this "stupid but provides compatibility with Vina" |
| `--outputmin` | — | Write `minout.sdf` of the minimization trajectory, with the given amount of interpolation |
| `--cnn_gradient_check` | false | Internal gradient checks |
| `--gpu_docking` | false | GPU acceleration for **non-CNN** scoring. Warns "experimental and not recommended" |

`--gpu_docking` and `--approximation gpu` together are the non-CNN GPU path. It is genuinely
experimental; the warning is not boilerplate.

## Verifying this page

This table drifts. To check it against your build:

```bash
# every option name the binary accepts
gnina --help --help_hidden 2>&1 | grep -oE '\-\-[a-zA-Z_0-9]+' | sort -u
```

Compare against the options named above. If they disagree, the binary is right.
