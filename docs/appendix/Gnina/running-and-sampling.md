# Running gnina, and how sampling works

**Audience:** users.

## Quickstart

The minimum that does something useful:

```bash
gnina -r receptor.pdb -l ligand.sdf --autobox_ligand reference_ligand.sdf -o docked.sdf.gz
```

Two additions are worth making reflexive:

- **`--seed <int>`.** Docking is a stochastic Monte Carlo process. Without a fixed seed you get a
  different answer — usually similarly good, not identical — on every run. A seed gives you
  reproducibility on the same machine and build; different systems or library versions may still
  differ. Left unset, the seed is derived from the process ID and the clock.
- **`-o out.sdf.gz`.** Without an output file you only see scores printed to the terminal, not the
  poses. Gzipped SDF is fully supported and is the format that preserves every score field — see
  [output-and-ranking.md](output-and-ranking.md).

## Common invocations

```bash
# Pocket docking
gnina -r rec.pdb -l lig.sdf --autobox_ligand orig.sdf -o docked.sdf.gz

# Whole-protein docking (see inputs-and-search-space.md)
gnina -r rec.pdb -l lig.sdf --autobox_ligand rec.pdb -o whole.sdf.gz --exhaustiveness 64

# Flexible side chains within 3.5 A of a reference ligand (read flexible-docking.md first)
gnina -r rec.pdb -l lig.sdf --autobox_ligand orig.sdf \
      --flexdist_ligand orig.sdf --flexdist 3.5 -o flex.sdf.gz

# Vinardo empirical scoring, no CNN
gnina -r rec.pdb -l lig.sdf --autobox_ligand orig.sdf \
      --scoring vinardo --cnn_scoring none -o vinardo.sdf.gz

# Single distilled CNN model instead of the 3-model ensemble: much cheaper
gnina -r rec.pdb -l lig.sdf --autobox_ligand orig.sdf --cnn fast -o fast.sdf.gz

# Minimize and score ligands already positioned in the pocket, no search
gnina -r rec.pdb -l posed.sdf --minimize -o minimized.sdf.gz
```

## The run modes

gnina is not only a docking program; several flags replace the search entirely, and the distinctions
matter because they change what `--cpu` and `--exhaustiveness` mean.

| Flag | What it does |
|---|---|
| *(none)* | Full docking: Monte Carlo search, refinement, scoring |
| `--score_only` | Score the pose exactly as given. No search, no minimization, no output file |
| `--minimize` | Energy-minimize the pose as given. This is what you want for rescoring pre-posed ligands |
| `--local_only` | Local search within the autobox. Usually you want `--minimize` instead |
| `--randomize_only` | Generate random poses, attempting to avoid clashes. No receptor needed |

`--minimize` silently changes several defaults to values suited to converging properly rather than to
fast approximate refinement: `--minimize_iters` becomes 10000, `--approximation` becomes `spline`,
`--factor` becomes 10, and `--force_cap` drops to 10 so clashing input structures relax gently
instead of exploding.

**`--local_only` and `--minimize` silently skip ligands** whose bounding extent exceeds 100 Å, with a
warning. In a large batch that is quiet pose loss — check your output count against your input count.

## Sampling internals

Each independent Monte Carlo chain samples the ligand's **degrees of freedom**: the 6 rigid-body
motions (translation in x/y/z, rotation as pitch/roll/yaw) plus the ligand's rotatable torsions.
Never bond lengths, bond angles or ring conformations.

Every chain starts from a **fully randomized** pose. The conformation you handed in is discarded
except for the parts that are never sampled, like ring pucker — there is no bias toward your input
geometry.

### What one Monte Carlo step actually does

Worth being precise about, because the usual summary ("perturb a random degree of freedom") is only
two-thirds right. The move is chosen uniformly at random over `2 + n_torsions` entities per ligand,
then `[GNINA1.0]`, `mutate_conf` in [mutate.cpp](../../gninasrc/lib/mutate.cpp):

- **Translation** — a random step inside a sphere of radius ~2 Å. A perturbation.
- **Whole-molecule rotation** — a quaternion increment scaled by the molecule's gyration radius. A
  perturbation.
- **One torsion** — set to a **fresh uniform random angle in [-π, π]**. Not a perturbation: the
  previous value is discarded entirely.

Because the choice is uniform over `2 + n_torsions` options, torsion moves dominate for flexible
ligands as an emergent consequence of the counting, not as a deliberate bias. A ligand with 10
rotatable bonds spends 10/12 of its moves re-randomizing a torsion.

Flexible-residue torsions join the same pool. That is the mechanism by which flexible docking makes
sampling harder: every side-chain torsion you add is another entity competing for the same moves. See
[flexible-docking.md](flexible-docking.md).

After the move, a fast local refinement runs using a **soft potential** — softened so that a move that
happens to overlap atoms does not generate an enormous repulsive force that would fling the ligand out
of the pocket. The result is then accepted or rejected by the Metropolis criterion at
`--temperature` (default 1.2).

During sampling the empirical scoring function is evaluated through a **precalculated grid
approximation**: one grid per ligand atom type, with the full ligand scored by interpolating and
summing per-atom values. This is what `--approximation` and `--factor` control. No such approximation
exists for CNN scoring, because a CNN score is not additive over atoms `[GNINA1.0]` — which is exactly
why `--cnn_scoring all` costs orders of magnitude more than the default.

### From chains to output

Each chain retains its best poses as it goes — the top 50 by default (`--num_mc_saved`) — even if the
chain later wanders elsewhere, with a 1 Å RMSD dedup applied during accumulation. At the end, the
chains' pools are merged and the top poses retained.

One correction to a commonly-repeated claim: the merged pool size is
**`max(--num_modes, --num_mc_saved)`**, not `--num_mc_saved` alone. Asking for more output poses than
50 therefore enlarges the internal pool too, which is why very large `--num_modes` values improve even
the *first* nine poses `[GNINA1.0]` — and why they cost more.

What happens to that pool — sorting, the diversity filter, and the `--num_modes` cutoff — is in
[output-and-ranking.md](output-and-ranking.md).

## Tuning knobs

**`--exhaustiveness`** (default 8) — the number of independent Monte Carlo chains, and the main way to
buy more sampling. The chains are fully independent and parallelize perfectly: 8 chains on 8 cores
costs about the same wall-clock time as 1 chain on 1 core.

For a targeted pocket, returns diminish quickly past 8: 8 → 16 doubles the compute for a small gain,
which is why 8 is the default `[GNINA1.0]`. This is **not** true for whole-protein docking — see
[inputs-and-search-space.md](inputs-and-search-space.md).

**`--cpu`** — how many cores gnina may use. Left unset it auto-detects, which on a shared machine is
usually the wrong thing.

> **Recommendation** — set `--cpu` explicitly on every job. `[Workshop2021]`

The relationship with `--exhaustiveness` depends on the mode, and this is easy to get backwards:

- **Docking**: parallelism comes from the chains. If `--cpu` is less than `--exhaustiveness`, chains
  queue and run in batches, and you lose the "more exhaustiveness is nearly free" property. Keep
  `--cpu` ≤ `--exhaustiveness`.
- **`--minimize` / `--local_only`**: there is no search, so `--exhaustiveness` is irrelevant. `--cpu`
  is per-ligand parallelism — each worker thread takes a whole ligand from a bounded queue. Set it to
  the cores you have.

**`--num_mc_steps` / `--max_mc_steps`** — by default the steps per chain come from a heuristic that
scales with the ligand: `105 × (50 + num_movable_atoms + 10 × degrees_of_freedom)`. Pinning it to a
constant is mostly useful for quick-and-dirty runs where you want a large batch of not-terrible random
poses fast rather than a converged search. `--max_mc_steps` caps the heuristic instead of replacing
it, which is the safer knob for bounding worst-case cost across a heterogeneous library.

**`--num_modes`** (default 9) — how many final poses to write, subject to the diversity filter. Not the
same thing as `--num_mc_saved`, though they interact as described above.

**Knobs that were measured not to matter.** `[GNINA1.0]` explored these and found no significant
effect on docking performance, so leave them alone unless you have a specific reason:
`--cnn_rotation` (default 0) and `--min_rmsd_filter` (default 1.0). Likewise
`--cnn_empirical_weight` with `--cnn_mix_emp_energy`/`--cnn_mix_emp_force` did not meaningfully change
results.

## Validating results with `obrms`

Use `obrms` (shipped with OpenBabel) to compute RMSD between poses. It correctly accounts for graph
isomorphism and symmetry — a benzene ring rotated 180° gives an RMSD of zero, as it should, rather than
a spuriously large number from naive atom-index matching. It is the tool both papers use to score
their benchmarks `[GNINA1.0]`.

The CLI differs between major versions: OpenBabel 2 takes a single hyphen for the "first-only" flag,
OpenBabel 3 takes two. Current installation instructions target OpenBabel 3, which is a substantial
improvement generally.

When you interpret your own RMSD numbers, keep the redock/cross-dock distinction from
[concepts.md](concepts.md) in view — a redocking RMSD is not evidence about prospective performance.

## Under the hood

- Monte Carlo driver: [monte_carlo.cpp](../../gninasrc/lib/monte_carlo.cpp); chain parallelism in
  [parallel_mc.cpp](../../gninasrc/lib/parallel_mc.cpp)
- The move set: `count_mutable_entities` and `mutate_conf` in
  [mutate.cpp](../../gninasrc/lib/mutate.cpp)
- Step-count heuristic, pool sizing and the `--minimize` default overrides: `main.cpp`, in the
  `parallel_mc` setup
- Grid approximation of the empirical function:
  [precalculate.h](../../gninasrc/lib/precalculate.h),
  [precalculate_gpu.cu](../../gninasrc/lib/precalculate_gpu.cu)
- Minimization: [quasi_newton.cpp](../../gninasrc/lib/quasi_newton.cpp),
  [bfgs.cu](../../gninasrc/lib/bfgs.cu)
