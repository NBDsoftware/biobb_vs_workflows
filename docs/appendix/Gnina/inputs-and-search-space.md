# Inputs and the search space

**Audience:** users. This is where most docking runs go wrong, and none of it produces an error
message.

## You always have to look at your data

Whatever you hand gnina as "the receptor" *is* the receptor, full stop. It does not sanity-check your
structure. If there is a zinc ion, an alternate-conformation lysine sitting on top of another
residue, or a leftover cofactor, gnina accepts it exactly as given — it has no way to know you did not
mean to include it. `[Workshop2021]`

### Receptor preparation

- Extract the protein atoms from your PDB file and go through them deliberately: missing atoms,
  alternate residue conformations, cofactors, binding-site waters. Decide to keep or discard each,
  rather than discovering later that the decision was made for you.
- **Remove the crystal ligand before docking.** This is the classic mistake. gnina does not
  distinguish "protein" from "ligand still sitting in the receptor file" — every atom in the receptor
  file becomes a rigid receptor atom, so the crystal ligand ends up permanently occupying the pocket
  you are trying to dock into.
- Prefer a receptor already in a bound (holo) conformation where you have the choice.
- **Protonation matters**, because hydrogen-bond donor/acceptor status depends on it, and the Vina
  scoring function has nothing else to go on (see [scoring-empirical.md](scoring-empirical.md)). By
  default, if you supply a receptor with no hydrogens, OpenBabel infers protonation. OpenBabel only
  ever *adds* hydrogens in chemically sensible places — it will not invent unrelated ones — but if you
  want exact control over hydrogen placement, supply the receptor as a **PDBQT file**: gnina takes
  PDBQT exactly as given and does not run it through OpenBabel at all.

### Ligand preparation

- Any format OpenBabel can read is accepted, including gzipped variants. But the input **must be a
  genuine 3D conformation**, not a 2D depiction laid out for drawing. gnina will happily dock a flat
  structure — bond lengths and angles all wrong — and may even report a good-looking score, but the
  pose is meaningless.
- You do *not* need the exactly right starting conformation, because rotatable torsions are sampled
  anyway (and in fact fully re-randomized; see
  [running-and-sampling.md](running-and-sampling.md)). What you do need is reasonable bond lengths,
  bond angles and protonation state, because none of those are ever sampled.
- **Ring conformations and stereochemistry are never sampled.** Whatever ring pucker or stereoisomer
  you provide is what comes back. To explore boat-versus-chair, or multiple stereoisomers, generate
  and dock multiple input conformers yourself. Docking several conformers per molecule is a
  documented way to improve results generally, not just for rings. `[GNINA1.3]`
- **Ligands much larger than about 20 Å across hit a CNN limitation**, because the scoring grid is
  23.5 Å on a side. See [scoring-cnn.md](scoring-cnn.md).

## Defining the box

There are three ways to tell gnina where to dock.

**1. Explicit coordinates** — `--center_x/y/z` plus `--size_x/y/z`.

**2. `--autobox_ligand <file>`** — draws a box around the bounding coordinates of the atoms in
`<file>`, then pads it symmetrically by `--autobox_add` (default 4 Å on all six sides).

**3. Flexible-residue definitions** contribute to the box too — see
[flexible-docking.md](flexible-docking.md), because that has a consequence people miss.

The autobox reference does not have to be a real ligand. It just needs atoms with Cartesian
coordinates. Alpha-sphere output from a pocket-detection tool like fpocket or MDpocket — literally a
cloud of unbonded "atoms" marking a cavity — works fine, as does a file of binding-site residue
coordinates. `[Workshop2021]`

### `--autobox_extend`

On by default. If any side of the auto-generated box is shorter than the longest distance between any
two atoms in the ligand, that side is extended to that distance `[GNINA1.0]`. The point is that poses
outside the box incur a penalty, so a long thin box would constrain the ligand's rotation as an
artefact of box geometry rather than of the binding site. Extending guarantees the ligand can rotate
freely in place.

This is one of the two things that make gnina differ from smina even with the CNN off; pass
`--autobox_extend 0` to reproduce smina's box exactly. Note it takes an *argument* — it is a
`value<bool>`, not a switch.

### `--autobox_add` is a real trade-off, in opposite directions

Flagged in the workshop as a trick question, and worth getting right `[Workshop2021]` `[GNINA1.0]`.

Increasing the box size does **not** slow gnina down — it does not change the number of sampled
degrees of freedom. But it does change docking *quality*, and which way depends on your task:

- **Redocking**: a smaller box is better. The crystal ligand's footprint is by construction the right
  answer, so a tight box is a helpful constraint that shrinks the search space. `[GNINA1.0]` measured
  redocking accuracy *decreasing* as `--autobox_add` grows.
- **Cross-docking**: a smaller box can hurt. The right box for the actual binding site may be
  considerably larger than the reference ligand used to define it. If your pocket can accommodate
  ligands bigger than your reference, a tight `--autobox_add` becomes a wrong constraint.

The papers are blunt about what this means for benchmarking: a low `--autobox_add` *artificially*
improves redocking numbers by unrealistically constraining the search to a space you only know
because you already have the answer. `[GNINA1.0]`

The default of 4 Å was chosen as the compromise. There is no universally correct value — ask whether
your case is closer to redocking (favour smaller) or to prospective work with potentially larger
ligands than your reference (favour larger).

## Whole-protein docking

Whole-protein docking is not a separate mode; it is a box choice. Point `--autobox_ligand` at the
entire receptor instead of a small reference ligand:

```bash
gnina -r rec.pdb -l lig.sdf --autobox_ligand rec.pdb -o whole_docked.sdf.gz --exhaustiveness 64
```

This is a legitimate way to dock without prior knowledge of the pocket, and it performs respectably —
but expect a large accuracy drop: Top1 falls from 73% to 38% for redocking and from 37% to 16% for
cross-docking relative to a defined pocket `[GNINA1.0]`. With the whole surface in play, the ligand
tends to settle in local minima far from the real site, and most of a protein surface is not
hospitable to binding, so once sampling finds *a* pocket it rarely leaves.

Two things behave differently from targeted docking:

**Exhaustiveness keeps paying.** Targeted docking shows clear diminishing returns past the default of
8. Whole-protein docking does not — the search space is much larger, so cranking `--exhaustiveness`
well beyond the default keeps helping. Concretely, going 8 → 16 raises Top9 from 48% to 58%
(redocking) and 22% to 29% (cross-docking) for whole-protein, versus 87% → 88% and 53% → 55% for a
defined pocket. `[GNINA1.0]`

> **Recommendation** `[GNINA1.0]` — for whole-protein docking, set exhaustiveness as high as your time
> budget allows.

**Better scoring beats more sampling.** The default CNN ensemble at exhaustiveness 16 (Top1 43%
redocking, 20% cross-docking) outperforms Vina scoring at exhaustiveness 64 (38% and 16%) — four times
the sampling, still behind. `[GNINA1.0]`

## Under the hood

- Box construction and `--autobox_extend`: `main.cpp`, around the autobox handling
- PDBQT passthrough and receptor parsing:
  [PDBQTUtilities.cpp](../../gninasrc/lib/PDBQTUtilities.cpp),
  [parse_pdbqt.cpp](../../gninasrc/lib/parse_pdbqt.cpp)
- `--custom_atoms` parsing (`setup_atomconstants_from_file` in `main.cpp`) takes ten fields per atom
  type: `ad_radius`, `ad_depth`, `ad_solvation`, `ad_volume`, `covalent_radius`, `xs_radius`,
  `xs_hydrophobe`, `xs_donor`, `xs_acceptor`, `ad_heteroatom`
