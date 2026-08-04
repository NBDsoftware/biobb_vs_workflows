# Concepts

**Audience:** users. Read this once, early. Everything else in the guide assumes it.

## Lineage: AutoDock Vina → smina → gnina

Understanding the lineage explains most of gnina's design choices, including some that look odd.

**AutoDock Vina** was written from scratch by Oleg Trott at Scripps, published in 2009, last updated
in 2011. Despite the name, it shares no code with the older AutoDock (dock4) — "Vina" is a distinct
program, not a shorthand for "AutoDock". Vina's focus was raw performance: a fast, empirically
parameterized scoring function optimized primarily for pose prediction, with affinity prediction a
secondary goal.

**smina** forked Vina around 2011 to add two features Vina's authors declined to merge upstream:
true energy minimization (converging to a local optimum of the scoring function, rather than Vina's
approximate refinement) and a pluggable scoring-function system. smina's defaults are designed to
reproduce Vina's behaviour as closely as possible, with one deliberate exception: a
hydrogen-optimization bug that smina fixed and Vina never did. smina is dual-licensed Apache/GPLv2 —
GPL because it links OpenBabel — and is still maintained, but its feature set is intentionally kept
stable and minimal. `[Workshop2021]`

**gnina** forked smina to add CNN scoring. Unlike smina, gnina does not aim for backward
compatibility with anything upstream; its goal is to be better than what came before. This is also
why it is a much heavier thing to build and deploy: on top of everything smina needed, it pulls in
CUDA, cuDNN and libtorch. `[Workshop2021]`

Two concrete differences from smina are worth knowing, because they mean gnina and smina do not
produce bit-identical output even with the CNN disabled `[GNINA1.0]`:

- gnina computes in single (32-bit) precision, smina in double (64-bit), because the CNN work has to
  move to the GPU efficiently. Measured effect on redocking: most output poses identical, small
  differences in some.
- gnina has `--autobox_extend`, which smina does not — see
  [inputs-and-search-space.md](inputs-and-search-space.md).

With `--cnn_scoring none` and `--autobox_extend 0`, gnina is otherwise the smina pipeline.

## A score is not an energy

The single most important framing: **call gnina's output a score, not an energy.** Docking scores
are fit to correlate with binding affinity, but they rest on approximations severe enough that
treating them as physical free energies will mislead you.

- **The receptor is normally rigid.** Even with flexible side-chain docking, only side-chain torsions
  move — the backbone never does. See [flexible-docking.md](flexible-docking.md).
- **The ligand is flexible only in its rotatable torsions.** Bond lengths, bond angles, ring
  conformations and stereochemistry are never sampled. Whatever you supply is what comes back.
- **There is no explicit solvent model** in the default (Vina) scoring function. `ad4_scoring` has a
  desolvation term, but that is a different scoring function, not a model of water molecules — see
  [scoring-empirical.md](scoring-empirical.md).

These approximations are exactly what make docking fast enough to be useful at scale. They are also
why "the score is a prediction of an affinity" is a more honest framing than "the score is an
energy". `[Workshop2021]`

The CNN scores carry the same caveat in a different form: `CNNaffinity` is reported in pK units and
is genuinely useful for ranking, but it is a neural network's guess conditioned on a 3D grid, not a
measured dissociation constant. See [scoring-cnn.md](scoring-cnn.md) for what the numbers are
calibrated to mean.

## The pipeline, end to end

Every gnina run follows the same six stages.

1. **Inputs** — a receptor, a ligand or multi-ligand file, and a definition of where to dock (the
   box). See [inputs-and-search-space.md](inputs-and-search-space.md).
2. **Sampling** — independent Monte Carlo chains explore ligand poses and torsions against the rigid
   receptor, each retaining the best poses it has seen. See
   [running-and-sampling.md](running-and-sampling.md).
3. **Refinement** — the surviving poses are locally minimized. By default this uses the empirical
   scoring function, in its exact functional form rather than the grid approximation used during
   sampling.
4. **Scoring** — final scores are computed. By default the CNN enters here and only here.
5. **Ranking and filtering** — poses are sorted by the chosen criterion, then a diversity filter
   drops any pose too close to a better-ranked one. See
   [output-and-ranking.md](output-and-ranking.md).
6. **Output** — the surviving poses are written out, each tagged with its scores.

The Monte Carlo sampling core is unchanged from Trott's original Vina implementation. Everything
gnina and smina added lives in the refinement, scoring and ranking layers wrapped around it.

One consequence worth internalizing early: **sampling and scoring fail differently.** If the correct
pose was never generated, no scoring function can rescue it — that is a sampling problem, addressed
with `--exhaustiveness` and a sensible box. If the correct pose was generated but ranked below a
wrong one, that is a scoring problem, addressed by a better scoring function. gnina's central claim
is that CNN rescoring buys more than the equivalent compute spent on extra sampling — quantified in
[performance.md](performance.md).

## Redocking, cross-docking, whole-protein

Almost every number in [performance.md](performance.md) is qualified by which of these tasks it
describes, and the differences are large. The vocabulary is worth learning before you read any
benchmark, including your own.

**Redocking** takes a ligand out of the complex it was solved in and docks it back into that same
receptor structure. It is easy to set up and easy to score, since the crystal pose is the answer.
It is also the best case by construction: the receptor is already in the right conformation for that
ligand. Useful as a sanity check and for comparing methods, but it systematically overstates what
you will get prospectively. `[GNINA1.0]`

**Cross-docking** docks a ligand into a *different* structure of the same protein — a non-cognate
receptor. This is the realistic case, because in real work the receptor was solved with some other
ligand, or none. Numbers are much lower than redocking, and it is what GNINA 1.3 deliberately
optimized for, accepting a small redocking regression to get it. `[GNINA1.3]`

**Whole-protein docking** dispenses with a known pocket and searches the entire protein surface.
Accuracy drops further again, and the tuning advice inverts — see
[inputs-and-search-space.md](inputs-and-search-space.md). `[GNINA1.0]`

The standard success metric across all three is **TopN**: the percentage of cases where a pose within
2 Å RMSD of the true pose appears in the top N ranked outputs. Top1 is the strict version and the one
usually quoted.

## What each release changed

**gnina 1.0 (2021)** `[GNINA1.0]` established the CNN-scoring pipeline and, importantly, derived its
own defaults rather than inheriting them: `--exhaustiveness 8`, `--autobox_add 4`, `--num_mc_saved 50`,
`--num_modes 9`, `--min_rmsd_filter 1.0`, `--cnn_rotation 0`. It selected a default 5-model CNN
ensemble by greedy forward search, and found that CNN *rescoring* captures essentially all the
benefit — CNN-guided refinement costs an order of magnitude more for no meaningful accuracy gain.
That negative result is why the cheap `rescore` mode is the default today.

**gnina 1.3 (2025)** `[GNINA1.3]` is three changes:

- **Caffe replaced by PyTorch.** No behavioural change for users, but a large CPU speedup — average
  docking time without a GPU fell from 129 s to about 30 s per complex, largely from better
  multiprocessing. It is also why you can now hand gnina an arbitrary TorchScript model with
  `--cnn_model`; see [developing.md](developing.md).
- **Models retrained on CrossDocked2020 v1.3**, which fixed ligand/receptor misalignment and
  incorrect bond typing in earlier versions. Cross-docking ranking improved across the board;
  redocking dipped slightly, because the dataset update filtered out problematic redocked poses.
- **Knowledge distillation**, which compresses an ensemble's ranking power into a single model. This
  is what `--cnn fast` is, and it is the reason a CNN-scored screen on CPU is now affordable. See
  [scoring-cnn.md](scoring-cnn.md).
- Plus a genuinely new capability: **covalent docking**. See
  [covalent-docking.md](covalent-docking.md).

> *Changed since 1.0:* the default CNN ensemble is 3 models, not 5. The 1.0 ensemble is still
> available as `--cnn default1.0`.

Independent evaluations consistently find gnina outperforming Vina and reaching performance similar
to commercial tools, and it has been used prospectively — notably in CACHE Challenge #1 `[CACHE1]`
and in a 7-million-compound screen. `[GNINA1.3]`

## Under the hood

- CLI definitions and the top-level pipeline: [gninasrc/main/main.cpp](../../gninasrc/main/main.cpp)
- Monte Carlo search: [monte_carlo.cpp](../../gninasrc/lib/monte_carlo.cpp),
  [mutate.cpp](../../gninasrc/lib/mutate.cpp)
- Minimization: [quasi_newton.cpp](../../gninasrc/lib/quasi_newton.cpp),
  [bfgs.cu](../../gninasrc/lib/bfgs.cu)
- Empirical scoring terms: [everything.h](../../gninasrc/lib/everything.h),
  [builtinscoring.cpp](../../gninasrc/lib/builtinscoring.cpp)
- CNN scoring: [cnn_torch_scorer.cpp](../../gninasrc/lib/cnn_torch_scorer.cpp),
  [torch_model.cpp](../../gninasrc/lib/torch_model.cpp),
  [dl_scorer.cpp](../../gninasrc/lib/dl_scorer.cpp)
