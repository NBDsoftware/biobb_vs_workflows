# gnina

gnina is a fork of AutoDock Vina (by way of smina) that adds convolutional
neural network (CNN) rescoring and GPU acceleration on top of the same
pose-sampling algorithm.

## How it works

1. **Sample** poses the same way Vina does (Monte Carlo search + local
   optimization).
2. **Rescore** the resulting poses with a CNN trained on protein-ligand
   structures, producing two extra scores per pose:
   - **`CNNscore`** (0-1) — how likely the pose itself is correct.
   - **`CNNaffinity`** (pK units, **higher is better**) — predicted binding
     affinity.
   - **`CNN_VS`** = `CNNaffinity × CNNscore`, a combined screening score.
3. **Rank and filter** poses — sorting happens *before* the redundancy
   filter, so changing the ranking metric or CNN model can change which
   poses survive, not just their order.

The empirical part of the score (`minimizedAffinity`) uses the same kind of
scoring function as Vina (`vina`, or the reparameterized `vinardo`, both
selectable via `--gnina_scoring`) and is in the same units (kcal/mol, lower
is better), so it stays comparable across engines even when the CNN score
doesn't.

## Choosing a ranking metric

The default is `CNNaffinity`
 (sort by predicted affinity, which is what ranks compounds in a screen). `CNNscore` answers a different question
(sort by network pose score, which answers whether a pose is right). This workflow ranks by `CNNaffinity` by default
(`--gnina_rank_by`), which is a reasonable default for screening.

## Cost

CNN settings dominate gnina's runtime far more than `--exhaustiveness` does.
Rough relative per-ligand cost on CPU, from gnina's own measurements:

| Setting | Relative cost |
|---------|---------------|
| `--gnina_cnn_scoring none` | ~1× (no CNN scores at all) |
| `--gnina_cnn fast` | ~3× |
| default 3-model ensemble | ~10× |
| `--gnina_cnn_scoring refinement` | ~100×, and does **not** improve pose prediction over the default `rescore` mode |

A GPU changes these numbers by roughly an order of magnitude. A common
strategy: screen a large library with `none` or `fast`, then re-dock the
best hits with the default ensemble.

## GPU and installation

A GPU is strongly recommended — gnina is Linux + NVIDIA only. `--gnina_no_gpu`
forces CPU execution, but that's a separate concern from having no CUDA
libraries at all: gnina's release binary is linked against CUDA/cuDNN and
**won't even load** on a machine without them. See the "Installing gnina"
section of the [virtual screening workflow](../workflows/virtual_screening.md)
for the actual setup steps.

## Limitations

- The CNN's input grid spans about 24 Å, so ligands larger than roughly 20 Å
  across will start to see artifacts in their CNN scores.
- Only the ligand is flexible; gnina supports flexible side chains but this
  workflow doesn't expose that option.
- CNN scores were not trained on covalent complexes or unusual geometries —
  they degrade with a poor input conformer more than Vina's empirical terms
  do, so well-prepared 3D ligands matter more here than with `vina`.

## References

- McNutt AT, Francoeur P, Aggarwal R, et al. *GNINA 1.0: molecular docking
  with deep learning.* J Cheminform. 2021;13:43.
  [doi:10.1186/s13321-021-00522-2](https://doi.org/10.1186/s13321-021-00522-2)
- McNutt AT, et al. *GNINA 1.3: the next increment in molecular docking with
  deep learning.* J Cheminform. 2025;17:28.
  [doi:10.1186/s13321-025-00973-x](https://doi.org/10.1186/s13321-025-00973-x)

For the CLI flags exposed by this repo (`--gnina_bin`, `--gnina_cnn_scoring`,
`--gnina_cnn`, `--gnina_scoring`, `--gnina_rank_by`, ...), see the
[virtual screening workflow](../workflows/virtual_screening.md).
