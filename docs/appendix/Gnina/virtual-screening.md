# Virtual screening

**Audience:** users. Read [performance.md](performance.md) alongside this.

## What to expect

gnina performs well on the DUD-E benchmark, and importantly it was **not trained on DUD-E**, so this is
a genuine out-of-domain evaluation. DUD-E has known biases that inflate results for methods trained on
it; gnina is not one of them, though the biases still complicate any comparison drawn on this benchmark.
`[GNINA1.3]`

Median across the 102 DUD-E targets `[GNINA1.3]`:

| Version | AUC | nEF1% |
|---|---|---|
| gnina 1.3 (default ensemble) | 0.78 | 0.27 |
| gnina 1.0 (default ensemble) | 0.75 | 0.25 |

1.3 improves on 1.0 for 68 of the 102 targets. The single `--cnn fast` model has comparable AUC to 1.0
but worse enrichment — a real trade-off, not a free lunch.

Earlier evaluation `[Workshop2021]`, `[VS1.0]` found gnina's default scoring and Vinardo both beating
plain Vina on median AUC and top-1% enrichment, with a consensus of the two better still.

## Which metric should you rank by?

This is worth flagging explicitly because **three credible sources give three different answers**, and
the guide is not going to pretend otherwise:

| Source | Metric | Rationale given |
|---|---|---|
| `[Workshop2021]` | **`CNNaffinity`** | Pose score answers "is this the right pose", not "is this a good binder". For ranking *compounds* you want the affinity head. |
| `[GNINA1.3]` | **`CNNscore`** | The published DUD-E numbers above are computed by ranking on the pose score. |
| `scripts/deepdock.py` | **`CNN_VS`** | The repository's own screening pipeline defaults `--target_metric` to `CNN_VS`, i.e. `CNNaffinity × CNNscore`. |

The existence of `CNN_VS` as a built-in output field (see
[output-and-ranking.md](output-and-ranking.md)) is what makes the disagreement legible: it is the
product of the two, so a compound needs both a plausible pose and a good predicted affinity to rank
highly. That is a defensible reading of what each head is for, and it is what the code's own pipeline
picked.

Practical guidance, given that state of affairs:

- If you are reproducing published numbers, use the metric that publication used.
- If you are running a real screen, `CNN_VS` is a reasonable default and is what the in-repo pipeline
  uses. Note that `--pose_sort_order` does **not** accept it, so ranking by `CNN_VS` means reading the
  SD field and sorting yourself.
- Whichever you choose, look at more than one. Which brings us to consensus.

## Consensus scoring

> **Recommendation** `[Workshop2021]` — whenever you have more than one scoring approach available, rank
> compounds by each separately and combine by **summing ranks**, not raw scores.

Ranks, not scores, because different tools and functions have incompatible scales and units and you
would otherwise have to reconcile them. A compound ranked highly by multiple differently-parameterized
methods is a stronger signal than any single method's top hit.

A worked illustration of why this matters, and why a single summary number misleads: on one DUD-E
target, Vinardo gave noticeably better early enrichment than CNN affinity — even though CNN affinity was
better on average across the whole benchmark. `[Workshop2021]` Look at full ROC curves, not just AUC.

## Running at scale

Practical recommendations, several of them opinionated `[Workshop2021]`:

**Pre-filter your library.** Remove highly flexible ligands especially. Docking cost scales with degrees
of freedom (the step-count heuristic is literally linear in them; see
[running-and-sampling.md](running-and-sampling.md)), and there are diminishing returns on spending that
compute on ligands unlikely to dock well anyway.

**Set `--cpu` explicitly on every job.** A commonly-seen failure on shared clusters is launching many
jobs that each auto-detect and try to claim *all* available cores. For docking, keep
`--cpu` ≤ `--exhaustiveness`; for `--minimize` batches the rule is different — see
[running-and-sampling.md](running-and-sampling.md).

**Batch ligands into one input file.** Do not launch one gnina process per ligand. Every invocation
re-reads and reprocesses the receptor through OpenBabel and pays process startup; feed a single
multi-ligand file so that cost amortizes across the batch. Pre-converting your receptor to PDBQT once
skips the OpenBabel receptor processing entirely on every run.

**Sharing a GPU across concurrent jobs is fine** — by default the GPU is only touched during final CNN
rescoring, not the whole run. But each concurrent process needs its own slice of VRAM, so if GPU memory
is the limiting factor you may not run as many concurrent jobs as you have cores.

**Pick the cheap CNN deliberately.** `--cnn fast` exists for this. `[GNINA1.3]` is specific about the
intended pattern: many single-threaded jobs with `--cnn fast`, followed by a rescreen of the top hits
with the default 1.3 ensemble to cut false positives. It also notes the flip side — when you have ample
compute or GPUs, the speedup does not justify the extra complexity of a hierarchical screen. Numbers in
[performance.md](performance.md).

**Consider a smaller library instead.** `[Workshop2021]`'s bottom line was blunt: *"my final
recommendation with high-throughput screening is to not do it"* — brute-force docking an entire large
library is expensive and there are often better uses of the same compute. The suggested alternative is a
cheap pre-filter followed by docking a shortlist. That is what the next section is about.

## Pharmacophore pre-filtering with Pharmit

[Pharmit](https://pharmit.csb.pitt.edu) is a much cheaper complement to brute-force docking. You define a
3D pharmacophore query — a spatial arrangement of hydrogen-bond donors and acceptors, hydrophobic points
and charged features, optionally with an exclusion shape to reject molecules overlapping the receptor —
and it searches pre-computed conformer libraries, including large commercial aggregators (tens of
millions of compounds) and user-contributed libraries, in seconds to minutes.

It also has a minimize step of its own: rigid-receptor, flexible-ligand minimization against the query
pose using the Vina scoring function. That gives you both a score and, usefully, **how far each hit moved
from its pharmacophore-aligned pose** — a large movement suggests the compound does not actually want to
make the interactions you asked for, even if its post-minimization score looks fine.

The trade-off: it is fast and lets you encode expert insight directly, but you need a *good* query, which
is hard without a bound reference ligand showing you real interaction geometry, and a single query only
tests one hypothesis about how the ligand binds.

The combined workflow `[Workshop2021]`:

1. Refine a pharmacophore query interactively until you have roughly 1,000–10,000 hits — few enough to be
   manageable, many enough not to be overly specific.
2. Download them.
3. Rescore with gnina (`--minimize`, and/or CNN scoring) to bring in the CNN, which Pharmit has no access
   to.

Because Pharmit already applies a Vina-based minimize/score step, it is worth switching the scoring
function for your gnina pass — Vinardo, say — to see whether it changes the picture. That is exactly the
single-target case mentioned above where Vinardo beat CNN affinity on early enrichment.

## Iterative screening with `deepdock.py`

[scripts/deepdock.py](../../scripts/deepdock.py) is an iterative, large-scale screening pipeline built on
gnina and a SPRINT model. It reads SMILES or Parquet input, generates 3D conformers, docks in batches,
and selects the next iteration's candidates from a combination of SPRINT-model predictions and
fingerprint similarity, writing Parquet output. Parallelism is via Dask (`dask.bag`,
`dask.dataframe`, `dask_jobqueue`, `dask.distributed`) with a cluster YAML — see
[scripts/slurm.yaml](../../scripts/slurm.yaml) — plus `multiprocessing`.

It is a subcommand CLI:

```
deepdock.py {prepare, initial_batch, next_batch, dock_batch, select_batch, analyze, topn, all} \
            --dir OUTDIR -r RECEPTOR --autobox_ligand LIG [...]
```

`all` skips already-completed steps; every other subcommand overwrites its outputs. Notable options:
`--target_metric` (default `CNN_VS`), `--batch_size` (default 100,000), `--num_batches` (default 5),
`--molweight_cutoff` (default 1200), `--cluster`, `--sprint_checkpoint`, `--target_sequence`,
`--smiles-only` (reduced Parquet with only `db`/`name`/`smiles`, skipping fingerprint and SPRINT
generation), `--max_workers`, `--iolimit`, and `--gnina_executable`.

Beyond gnina itself it needs `dask`, `distributed`, `dask_jobqueue`, `torch`, `rdkit`, `xgboost`,
`ultrafast`, `pyarrow`, `pandas`, `numpy`, `scipy`, `matplotlib`, `pyyaml` and `biopython`.

> **Status:** this script is the most actively developed part of the repository, and it is documented
> only by its own `--help` — [scripts/README.md](../../scripts/README.md) covers `makeflex.py` only.
> Treat the option list above as a snapshot and check `--help` before relying on it.

## Under the hood

- `CNN_VS` computation: `result_info::write` in [result_info.cpp](../../gninasrc/lib/result_info.cpp)
- Sort-order handling: `main.cpp`, in the per-ligand output block
- Receptor reprocessing happens per invocation in `main.cpp`'s setup path, which is what batching avoids
