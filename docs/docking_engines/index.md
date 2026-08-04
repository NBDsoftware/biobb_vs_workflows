# Docking engines

The `virtual_screening` workflow supports two docking engines, selected with
`--docking_engine`. Both use the same pose-sampling algorithm and the same
docking box, so their results are directly comparable.

| Engine | Scoring | Units (lower/higher better) | Relative speed |
|--------|---------|------------------------------|-----------------|
| [`vina`](vina.md) (default) | Empirical function | kcal/mol, lower is better | Fast |
| [`gnina`](gnina.md) | Empirical function + CNN rescoring | `CNNaffinity` (pK), higher is better | Slower, GPU recommended |

Use `vina` for a fast baseline with no extra setup. Use `gnina` when you have
a GPU (or can tolerate slower CPU runs) and want a learned pose/affinity
score on top.

See the [virtual screening workflow](../workflows/virtual_screening.md) for
the CLI flags, and the pages above for the scientific background on each
engine's scoring function.

```{toctree}
:hidden:

vina
gnina
```
