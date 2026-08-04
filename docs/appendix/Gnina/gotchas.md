# Gotchas cheat sheet

**Audience:** everyone. Skim once before your first run; re-read after your first surprise. Each item
links to the page with the full explanation.

## Inputs

- **Always remove the crystal ligand** from the receptor file before docking — gnina treats every atom
  you give it as receptor, with no way to tell "this was the bound ligand." →
  [inputs-and-search-space.md](inputs-and-search-space.md)
- **The receptor is accepted exactly as given** — alternate conformations, waters, cofactors, metal
  ions, all of it. Curate the structure yourself. →
  [inputs-and-search-space.md](inputs-and-search-space.md)
- **Check protonation.** No-hydrogen input goes through OpenBabel's inference; supply a PDBQT for
  exact control, since gnina passes PDBQT through unmodified. →
  [inputs-and-search-space.md](inputs-and-search-space.md)
- **Ligand input must be a real 3D conformer**, not a 2D depiction — bond lengths and angles are never
  sampled and must already be correct. → [inputs-and-search-space.md](inputs-and-search-space.md)
- **Ring pucker and stereochemistry are never sampled** — generate and dock multiple conformers
  yourself if you need to explore them. → [inputs-and-search-space.md](inputs-and-search-space.md)

## Running

- **Set `--seed`** for reproducibility, and **set `-o`** — without it you only get scores printed to
  the terminal, not poses. → [running-and-sampling.md](running-and-sampling.md)
- **Gzip your SDF output** (`.sdf.gz`) — fully supported, and the only format that carries every score
  field, not just size savings. → [output-and-ranking.md](output-and-ranking.md)
- **Set `--cpu` explicitly.** For docking, keep it ≤ `--exhaustiveness`; for `--minimize`/`--local_only`
  the relationship is different — `--cpu` is per-ligand parallelism there and `--exhaustiveness` does
  not apply. → [running-and-sampling.md](running-and-sampling.md)
- **Batch ligands into one multi-ligand file** instead of one job per ligand, for any large-scale run.
  → [virtual-screening.md](virtual-screening.md)
- **`--addH`, `--stripH`, `--autobox_extend` require an argument** — they are not switches. → [cli-reference.md](cli-reference.md)
- **`--local_only`/`--minimize` silently skip ligands** wider than 100 Å. → [running-and-sampling.md](running-and-sampling.md)

## Scoring and output

- **`--cnn_scoring none` silently forces `--pose_sort_order Energy`**, regardless of what you passed.
  → [output-and-ranking.md](output-and-ranking.md)
- **Changing `--pose_sort_order` (or the CNN model) can change *which* poses you get**, not just their
  order — sorting happens before the diversity filter and the `--num_modes` cutoff.
  → [output-and-ranking.md](output-and-ranking.md)
- **`CNNaffinity_variance` only appears for ensembles** — a single-model run (e.g. `--cnn fast`) omits
  it entirely. → [output-and-ranking.md](output-and-ranking.md)
- **PDBQT output omits `CNN_VS` and the variance** — use SDF if you need them. →
  [output-and-ranking.md](output-and-ranking.md)
- **There is no single agreed ranking metric for virtual screening** — CNNaffinity, CNNscore and
  CNN_VS each have a credible source recommending them. → [virtual-screening.md](virtual-screening.md)
- **`--device` is ignored by the PyTorch backend** — use `CUDA_VISIBLE_DEVICES`. →
  [install.md](install.md)

## Flexible and covalent docking

- **Default to rigid-receptor docking**; flex only specific, structurally-justified side chains. →
  [flexible-docking.md](flexible-docking.md)
- **Flexible residues enlarge the autobox**, on top of adding degrees of freedom. →
  [flexible-docking.md](flexible-docking.md)
- **`--flex_limit`/`--flex_max` are ignored when `--flexres` is given.** →
  [flexible-docking.md](flexible-docking.md)
- **Covalent docking: score with `--cnn_scoring none`.** The CNN was not trained on covalent
  complexes, and gnina warns you about exactly this. → [covalent-docking.md](covalent-docking.md)
- **Every SMARTS match is docked** in covalent mode — a loose pattern silently multiplies your output
  and can crowd out the attachment point you meant to test. →
  [covalent-docking.md](covalent-docking.md)

## Installing and building

- **On a GPU-less machine, `--no_gpu` is not enough on its own** — the release binary is linked
  `BIND_NOW` against CUDA/cuDNN and will not load at all unless those libraries are present. →
  [install.md § Running with no NVIDIA GPU](install.md#running-with-no-nvidia-gpu)
- **`-DCMAKE_CUDA_ARCHITECTURES=all` is overridden** by a non-cache `set()` in `CMakeLists.txt` that
  forces `all-major`. → [developing.md](developing.md)
- **pygnina's `minimize` is an unexposed `abort()` stub** — only `set_receptor` and `score` work. →
  [developing.md](developing.md)
