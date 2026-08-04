# Performance: accuracy and runtime

**Audience:** users deciding whether, and how, to use gnina. Every number below is labelled with its
dataset and its source — the redocking/cross-docking/whole-protein distinction from
[concepts.md](concepts.md) changes the numbers by a large factor, so a figure without that label is not
usable.

Metric: **TopN**, the percentage of cases with a pose within 2 Å RMSD of the true pose in the top N
outputs. Top1 is the strict version and the one usually quoted.

## Redocking and cross-docking, defined pocket

`[GNINA1.0]`, Top1, pocket defined by `--autobox_ligand`:

| Task | Vina | gnina (1.0 default ensemble) |
|---|---|---|
| Redocking | 58% | **73%** |
| Cross-docking | 27% | **37%** |

## Whole-protein docking

`[GNINA1.0]`, Top1, box = entire receptor, exhaustiveness 8:

| Task | Vina | gnina (1.0 default ensemble) |
|---|---|---|
| Redocking | 31% | **38%** |
| Cross-docking | 12% | **16%** |

Both drop sharply from the defined-pocket numbers above — expected, given the much larger search space.
See [inputs-and-search-space.md](inputs-and-search-space.md) for why exhaustiveness behaves differently
here.

## 1.0 → 1.3: the default ensemble

`[GNINA1.3]`, Wierbowski cross-docking dataset and Posebusters/Astex redocking sets:

| Metric | gnina 1.0 | gnina 1.3 |
|---|---|---|
| Cross-docking Top1 | 37% | **40%** |
| Redocking Top1 (Posebusters) | **69%** | 67% |
| CPU-only time per complex | 30 s | **23 s** |

The redocking dip is real and the paper explains it: CrossDocked2020 v1.3 filtered out problematic
redocked poses that the v1.0 models had been (over-)fit to, and cross-docking was deliberately
prioritized as the task that matters prospectively. `[GNINA1.3]` calls this "a sensible strategy" given
that redocking is a synthetic benchmark relative to real prospective use. A regression on the easier,
less realistic task in exchange for a gain on the harder, realistic one is a reasonable trade — but it
is one to be aware of if you were comparing against 1.0 numbers on a pure redocking benchmark.

## `--cnn fast`: the distilled single model

`[GNINA1.3]`, same datasets:

| Metric | 1.0 default ensemble (5 models) | `--cnn fast` (1 model) |
|---|---|---|
| Cross-docking Top1 | 37% | ~36% |
| Redocking Top1 | 69% | 64% |
| CPU-only time per complex | 30 s | **16 s** |
| vs. Vina empirical alone | — | +1.3 s |
| vs. `fast` on GPU | — | <1 s slower |

`fast` sacrifices a small amount of accuracy — larger on redocking than cross-docking — for a large
runtime win, and lands within noise of the CPU cost of empirical scoring alone. It sits, along with the
1.3 default ensemble, on the Pareto frontier of the accuracy/cost trade-off `[GNINA1.3]`; see
[scoring-cnn.md](scoring-cnn.md) for what it actually is (a knowledge-distilled model) and when to reach
for it.

## Caffe → PyTorch, independent of model choice

Same models, same weights, different deep-learning backend `[GNINA1.3]`:

| | Caffe (1.0) | PyTorch (1.3) |
|---|---|---|
| CPU-only, per complex, 4 cores | 129 s | **~30 s** |

No change in pose-prediction accuracy — this is purely an inference-engine speedup, mostly from better
multiprocessing support. Reported with 4 cores requested; the paper notes the benefit is likely larger
still on many-core systems.

## Virtual screening (DUD-E)

`[GNINA1.3]`, median across 102 targets, ranked by CNNscore (see
[virtual-screening.md](virtual-screening.md) for why the ranking metric itself is contested):

| Version | AUC | nEF1% |
|---|---|---|
| gnina 1.0 | 0.75 | 0.25 |
| gnina 1.3 | **0.78** | **0.27** |
| `--cnn fast` (1.3) | ~0.78 | worse than default |

1.3 improves on 1.0 for 68 of 102 targets. gnina was not trained on DUD-E, so this is a genuine
out-of-domain evaluation, though DUD-E's known dataset biases still complicate cross-method comparison.

Earlier, method-comparison evaluation `[VS1.0]`, `[Workshop2021]`: gnina's default scoring and Vinardo
both beat plain Vina on median AUC and top-1% enrichment; a consensus of the two did better still.

## Covalent docking

`[GNINA1.3]`, Scarpino 207-complex benchmark, success = top pose within 2 Å RMSD. Full detail and the
scoring-function recommendation in [covalent-docking.md](covalent-docking.md):

| Setting | Success |
|---|---|
| Covalent mode, generated conformer, no positioning info | 36.2% |
| Covalent mode, experimental conformer + specified position | 66.6% |
| Covalent mode off, CNN scoring | 27.5% |
| Covalent mode off, Vina scoring | 15.8% |

## CNNscore calibration

`[GNINA1.0]`, fraction of top poses within 2 Å RMSD, by CNNscore threshold:

| | Redocking | Cross-docking |
|---|---|---|
| Fraction of top poses within 2 Å, given CNNscore > 0.8 | ~79% | ~56% |
| Fraction of top poses that clear a CNNscore of 0.8 at all | 87% | ~15% |

The second row is the one to remember day to day: a low CNNscore on a real, non-cognate target is
the expected outcome, not evidence of a broken run. See [scoring-cnn.md](scoring-cnn.md).

## Locally measured CPU timings

The numbers above are all from the papers' benchmark hardware (dual 16-core Xeon 5218, RTX 2080 Ti,
CUDA 10.2/cuDNN 7.6.5 for the 1.0 paper). The table below is a **single local measurement**, included
because it is concrete rather than because it generalizes — treat it as one data point on one machine,
not a benchmark.

Machine: 12-core laptop, no NVIDIA GPU, `--cpu 8 --exhaustiveness 8 --seed 0`, redocking
`test/gnina/data/184l_lig.sdf`, all 9 poses written. Full setup in
[local/cpu-only-workstation.md](local/cpu-only-workstation.md).

| Mode | Wall clock | Peak RSS |
|---|---|---|
| `--cnn_scoring none` (Vina empirical only, no torch) | 2.3 s | 0.42 GB |
| `--cnn fast` (single distilled model) | 6.1 s | 0.45 GB |
| Default (`--cnn_scoring rescore`, 3-model ensemble) | 23.7 s | 0.56 GB |

Consistent with the published `--cnn fast` positioning: cheapest first is `--cnn_scoring none` for
large batches, `--cnn fast` as the middle ground when a CNN score is wanted per pose, the default
ensemble for one-off runs. Avoid `--cnn_scoring refinement` (roughly 10× `rescore`) and `all` entirely
on CPU.
