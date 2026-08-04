# Empirical scoring functions

**Audience:** users and developers.

The empirical scoring functions are weighted sums of mostly-pairwise potential terms. The term
implementations live in [everything.h](../../gninasrc/lib/everything.h); the named functions and their
fixed weights are assembled in [builtinscoring.cpp](../../gninasrc/lib/builtinscoring.cpp).

Two flags let you inspect the machinery: `--print_terms` and `--print_atom_types` dump the available
terms and atom types, and raising `--verbosity` prints the weighted term list actually in use.

## The built-in functions, term by term

These are the real weights from the source, not a paraphrase, because the differences between the
functions are structural rather than cosmetic.

**`--scoring vina`** (the default; `--scoring default` is an alias for it):

| Weight | Term |
|---|---|
| −0.035579 | `gauss(o=0,_w=0.5,_c=8)` |
| −0.005156 | `gauss(o=3,_w=2,_c=8)` |
| 0.840245 | `repulsion(o=0,_c=8)` |
| −0.035069 | `hydrophobic(g=0.5,_b=1.5,_c=8)` |
| −0.587439 | `non_dir_h_bond(g=-0.7,_b=0,_c=8)` |
| ~1.923 | `num_tors_div` |

Read as shapes rather than numbers, that is: a van-der-Waals-like attractive well built from two
gaussians plus a repulsion wall, zeroed where the two atoms' radii exactly touch; a more favourable
version of the curve when both atoms are hydrophobic; and a hydrogen-bond term in which the atoms are
actually *preferred* to overlap slightly, so the donor/acceptor distance shrinks. That last one has a
discontinuity built in.

`num_tors_div` is the one non-pairwise term: it penalizes the ligand for having more rotatable bonds,
standing in for the conformational entropy lost on binding.

**What Vina does not have.** No electrostatic terms and no partial charges. Only whether an atom is a
donor, an acceptor, or neither matters — a formally +2 atom that happens to be a donor scores
identically to any other donor. There is no metal term either: metals are modelled as
hydrogen-bond donors, and you can verify this directly in the atom-type table, where every metal type
(Mg, Mn, Zn, Ca, Fe, and the catch-all `GenericMetal`) carries `donor = true`.

Note that this limitation is **specific to Vina**, not to gnina — see `ad4_scoring` below.

**`--scoring vinardo`** — an independent reparameterization by a different group `[GNINA1.0]`, and
*not* simply Vina with different weights:

| Weight | Term |
|---|---|
| −0.045 | `gauss(o=0,_w=0.8,_c=8)` |
| 0.80 | `repulsion(o=0,_c=8)` |
| −0.035 | `hydrophobic(g=0.0,_b=2.5,_c=8)` |
| −0.60 | `non_dir_h_bond(g=-0.6,_b=0,_c=8)` |
| 0.0 | `num_tors_div` |

Vinardo **drops Vina's second gaussian entirely** and widens the first (`w=0.5` → `w=0.8`), so it is
four pairwise terms rather than five. It also reshapes the hydrophobic term substantially
(`g=0.5,b=1.5` → `g=0.0,b=2.5`) and ships its own complete atom-parameter table with different radii
per type. gnina adds a Boron type that Vinardo never parameterized, using defaults.

Note the `num_tors_div` weight: both functions compute it as `5 × w / 0.1 − 1`, and Vinardo's `w` of
0.02 makes that **exactly zero**. Vinardo therefore applies no rotatable-bond penalty at all, where
Vina's works out to about 1.923. If you are comparing scores between the two, that is a large
structural difference for flexible ligands, not a rounding detail.

Whether it beats Vina is system-dependent: "in my hands it's a little mixed which one's better"
`[Workshop2021]`. For virtual screening specifically, DUD-E evaluation found Vinardo ahead of default
Vina on average.

> **Recommendation** `[Workshop2021]` — try Vinardo alongside the default, and if you have both, combine
> them by consensus rather than picking one. See [virtual-screening.md](virtual-screening.md).

**`--scoring ad4_scoring`** — an implementation of the AutoDock 4 scoring function, and the reason the
"no electrostatics, no solvent" statement has to be scoped to Vina:

| Weight | Term |
|---|---|
| 0.1560 | `vdw(i=6,_j=12,_s=0,_^=100,_c=8)` |
| 0.0974 | `non_dir_h_bond_lj(o=-0.7,_^=100,_c=8)` |
| 0.1159 | `ad4_solvation(d-sigma=3.5,_s/q=0.01097,_c=8)` |
| 0.1465 | `electrostatic(i=1,_^=100,_c=8)` |
| 0.2744 | `num_tors_add` |

It has both an **electrostatic** term and a **desolvation** term. It is a 6-12 Lennard-Jones potential
rather than Vina's gaussians. If you need charge-aware scoring, this is the built-in that has it.

**`--scoring dkoes_scoring`, `dkoes_scoring_old`, `dkoes_fast`** — kept mainly for
paper-reproducibility. They are small functions built on `vdw(i=4,j=8)`, an h-bond term, a torsion
penalty and a constant offset, with `dkoes_scoring` adding desolvation. Not recommended for new work.

## Custom scoring functions

You can assemble your own function from the same term library, as long as every term is already
implemented in [everything.h](../../gninasrc/lib/everything.h) — the file format selects and weights
existing terms, it does not define new ones.

The format ([custom_terms.cpp](../../gninasrc/lib/custom_terms.cpp)) is deliberately simple:

- One term per line, **weight first, then the parameterized term name**, whitespace-separated.
- Everything after the term name is ignored, which makes it a convenient place for inline commentary.
- Lines whose first character is `#` are skipped, as are blank and whitespace-only lines.
- An unrecognized term name is an error (`Unknown term`), as is an unparseable parameter
  (`Could not convert parameters`).

[examples/kitchensink.score](../../examples/kitchensink.score) shows the full range of available terms
and their parameters, and uses the trailing-comment property to document them:

```
1.0  gauss(o=0,_w=0.5,_c=8)		o is offset, w is width of gaussian
1.0  electrostatic(i=1,_^=100,_c=8)	i is the exponent of the distance, see everything.h
```

Load one with `--custom_scoring <file>`, plus `--custom_atoms <file>` if you also need custom
atom-type parameters (ten fields per type; see
[inputs-and-search-space.md](inputs-and-search-space.md)).

### Dumping a per-term breakdown

Combined with `--score_only`, a custom scoring file is the practical way to get per-term values for a
pose. gnina prints a `#`-prefixed header row of term names followed by one row per scored structure,
which is effectively a whitespace-separated table:

```
## Name gauss(o=0,_w=0.5,_c=8) gauss(o=3,_w=2,_c=8) repulsion(o=0,_c=8) ...
## 184L_I4B_A_401 45.89304 721.11029 0.55266 ...
```

Grep for `^##` and you have a CSV. `[Workshop2021]`

### The hacks, and what replaced them

`[Workshop2021]` described two uses of custom terms as "genuinely hacky but effective":

- **Custom metal coordination** — if you have a single iron in the receptor and want a specific ligand
  nitrogen to coordinate with it, add a strongly attractive nitrogen–iron term. Sampling then strongly
  prefers that geometry. Still a reasonable trick.
- **Soft covalent docking** — add a strongly attractive gaussian between, say, a chlorine and a
  sulfur atom type, then arrange the system so the only chlorine in the ligand and the only sulfur in
  the receptor are the pair you want pulled together. This does not create a bond; it just makes the
  proximity so favourable that sampling insists on it. There is no physically real chlorine/sulfur
  interaction being modelled — you are borrowing an existing term to fake a local preference.

  > *Changed since the workshop:* gnina 1.3 has real covalent docking. Use the `--covalent_*` options
  > instead; see [covalent-docking.md](covalent-docking.md). The hack is only of interest now if you
  > want a soft bias rather than an actual bond.

## Building a target-specific empirical model

If you have a labelled dataset for one target — actives and decoys, or a congeneric series with known
affinities — the term library gives you interpretable features to train on:

1. Score your labelled examples with a function that has every term active at weight 1.0, e.g.
   [examples/kitchensink.score](../../examples/kitchensink.score) via `--custom_scoring` with
   `--score_only`.
2. Extract the per-term breakdown as above into a table of numeric features per example.
3. Train a conventional classifier or regressor — logistic regression is enough — on those features
   against your labels.

`[Workshop2021]` demonstrated this on DUD-E-style active/decoy labels for a single target and got a
much better cross-validated AUC than any general-purpose scoring function or consensus. The caveat
that came with it is the important part: **results that good from a small, simple model should make
you suspicious rather than confident.**

- Such a model tends to learn incidental properties that separate *that particular* active/decoy set
  — rigidity, for example, or, surprisingly and without a fully satisfying explanation, a preference
  for molecules with technically unfavourable electrostatic and hydrogen-bond geometry (possibly
  because "good" molecules can tolerate a clash that "bad" molecules cannot). Validate it hard —
  scaffold-based holdouts, comparison against a plain similarity search to the training set — before
  trusting it on a new library.
- The genuine advantage over training a CNN for the same purpose is **interpretability**: the features
  are named physical terms, so you can inspect the weights and sanity-check what was learned. That is
  exactly how the "likes electrostatic clashes" red flag above was found.

> **Recommendation** `[Workshop2021]` — never use a model trained this way to guide docking or pose
> optimization, i.e. never plug it in as `--custom_scoring` for search or minimization. Its training
> data is (molecule, label) pairs with no information about *bad poses of the same molecule*. A
> scoring function used during search needs to distinguish good from bad poses of one ligand, and an
> active-versus-decoy classifier has no such signal. Use it only to re-rank poses generated by a
> proper pose-aware function.

That generalizes into a useful sanity check for any structure-based scoring model: it should be
**pose-sensitive.** A genuinely good active placed in an obviously wrong pose should score badly. If it
does not, the model is not doing structure-based prediction.

## Under the hood

- Term implementations and their parameter grammar: [everything.h](../../gninasrc/lib/everything.h)
- Named functions, weights and the `vinardo_data` atom table:
  [builtinscoring.cpp](../../gninasrc/lib/builtinscoring.cpp)
- Custom-file parsing, `custom_terms::add`: [custom_terms.cpp](../../gninasrc/lib/custom_terms.cpp)
- Grid precalculation used during sampling: [precalculate.h](../../gninasrc/lib/precalculate.h)
