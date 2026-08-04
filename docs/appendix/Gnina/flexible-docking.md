# Flexible side-chain docking

**Audience:** users.

gnina, like smina and Vina before it, can treat selected side chains as flexible. **The backbone is
always rigid** — only side-chain torsions are sampled, never backbone motion. Flexible residues are
selected once at the start of docking and the selection is never updated during sampling.
`[GNINA1.0]`

## Read this before you use it

> **Recommendation** `[GNINA1.0]` — use rigid-receptor docking by default. Reach for flexible side
> chains only when you have specific structural evidence that a particular residue needs to move, and
> hand-pick it.

This is the published conclusion, not a matter of taste. `[GNINA1.0]` evaluated flexible against rigid
cross-docking across ~7,900 complexes, binning results by how much the binding site actually differs
between the docking target and the cognate receptor:

- For **low** target–cognate side-chain RMSD — a highly similar binding site, so effectively a
  redocking situation — rigid docking is better on average, as you would expect.
- For **higher** RMSD, where you would expect flexibility to pay, the picture is unclear. Flexible
  docking looks equivalent or slightly better on average, but there are few systems with target–cognate
  RMSD above 6 Å, so the apparent improvement is inconclusive.
- Overall RMSD distributions for top poses are similar between the two, with slightly *more* low-RMSD
  systems for rigid docking.

Success is very system-dependent, and flexibility is never free. Three costs compound:

1. **More degrees of freedom.** Every side-chain torsion joins the same pool of Monte Carlo moves as the
   ligand's own torsions (see [running-and-sampling.md](running-and-sampling.md)), so ligand sampling
   gets diluted. The default step-count heuristic scales with degrees of freedom, so the run also just
   costs more.
2. **A larger search box.** Flexible side chains are included in the calculation of the autobox bounds
   `[GNINA1.0]`. Turning on flexibility silently enlarges your search space, on top of adding DOF.
3. **Models trained for rigid docking.** The default CNN models were not trained on flexible-docking
   data. `[GNINA1.0]` is explicit that optimizing defaults and training models for this task was left
   to future work — and as of 1.3 that work has not landed.

> **Recommendation** `[GNINA1.0]`, `[Workshop2021]` — if backbone flexibility genuinely matters for your
> system, dock against an *ensemble* of receptor conformations (NMR models, multiple PDB entries) rather
> than reaching for side-chain flexibility. gnina cannot sample backbone motion at all, so an ensemble
> is the only way to represent it.

## Selecting residues

Four mechanisms, in rough order of how much judgement they ask of you:

**`--flexres <chain:resid,...>`** — an explicit, manually curated list, e.g. `A:300`. An insertion code
may be appended. This requires you to look at the structure and decide which side chains genuinely need
to move. More setup effort, and the recommended approach whenever you have any reason to believe
specific residues matter.

**`--flex <file.pdbqt>`** — flexible side chains defined directly in a PDBQT file. The most explicit
option, and how Vina originally did it.

**`--flexdist_ligand <file>` + `--flexdist <distance>`** — automatically flex every residue with a side
chain atom within the given distance of a reference ligand. Convenient for systematic or batch
evaluations, at the cost of specificity. `[GNINA1.0]` used 3.5 Å for its own evaluation as a reasonable
representation of a binding site — which is also the setting whose results argued for defaulting to
rigid.

**`--flex_limit` / `--flex_max`** — bounds for pipeline use:
- `--flex_limit N` is a *hard* cap: if more than N residues are selected, gnina refuses to dock and
  terminates with a warning. Useful for avoiding bottlenecks in a large screen `[GNINA1.0]`.
- `--flex_max N` is a *soft* cap: keep only the N residues closest to `--flexdist_ligand`.

Two behaviours neither the help text nor the papers mention:

- The two are **mutually exclusive** — passing both is an error.
- Both are **silently ignored when `--flexres` is given** (with a warning). If you hand-pick residues,
  you are trusted with the count.

## Output

Since most of the receptor is unchanged, gnina writes out **only the flexible side chains** by default
rather than duplicating the whole protein for every pose. That is much more efficient, but it means the
flexible output is not a usable structure on its own.

- `--out_flex <file>` — where the flexible residues go.
- `--full_flex_output` — write the entire structure instead, if you would rather have that.

To reassemble a complete receptor afterwards:

```bash
python scripts/makeflex.py RIGID.pdb FLEXIBLE.pdb OUT.pdb
```

It handles multi-model `FLEXIBLE.pdb` files (multiple `MODEL`/`ENDMDL` records). See
[scripts/README.md](../../scripts/README.md).

## A caveat on connectivity

`[GNINA1.0]` found that disulfide bonds between cysteine residues are **allowed to break during
sampling**, with a software warning. In their 7,970-system evaluation this produced different
connectivity in the output for four systems, plus five more with unexplained connectivity changes — all
nine discarded from analysis. If you are flexing cysteines, check the connectivity of your output.

## Under the hood

- Flexible-residue selection, caps and the box interaction: `main.cpp`, in the flex setup block
- Side-chain torsions enter the same move pool via `count_mutable_entities` in
  [mutate.cpp](../../gninasrc/lib/mutate.cpp)
- PDBQT flexible-residue parsing: [parse_pdbqt.cpp](../../gninasrc/lib/parse_pdbqt.cpp),
  [PDBQTUtilities.cpp](../../gninasrc/lib/PDBQTUtilities.cpp)
