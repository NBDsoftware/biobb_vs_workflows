# Covalent docking

**Audience:** users. New in gnina 1.3 `[GNINA1.3]`.

## What gnina does and does not do

gnina does **not** model a chemical reaction. Like most covalent docking tools, it expects you to supply
the ligand **already in its bound covalent form**, including whatever chemical modification the reaction
would have produced — an epoxide already opened, a Michael acceptor already added across, and so on.
You then tell gnina which ligand atom and which receptor atom are bonded, and it positions and samples
the resulting construct.

Given a receptor/ligand atom pairing, gnina:

1. Repositions the ligand so the covalent ligand atom sits within bonding distance of the receptor atom.
2. Creates the bond, at `--covalent_bond_order` (default 1).
3. Treats the **residue-plus-ligand as a single flexible residue**. Its internal torsions are sampled and
   minimized during Monte Carlo search and refinement, but **no rigid-body transformation is applied** —
   the construct is anchored to the protein, so translation and rotation of the ligand as a unit are
   gone. This is the key mechanical difference from ordinary docking.
4. For CNN scoring, keeps ligand atoms typed as ligand atoms, so the receptor/ligand distinction the
   network relies on is preserved.

Positioning uses OpenBabel's `GetNewBondVector` heuristic on the receptor atom (after reducing its
hydrogen count) to find a sensible placement, then `OBBuilder::Connect` to rotate and translate the
ligand into a reasonable bonding geometry.

## Specifying the atoms

```bash
gnina -r rec.pdb -l covalent_form_of_lig.sdf \
      --autobox_ligand orig.sdf \
      --covalent_rec_atom A:145:SG \
      --covalent_lig_atom_pattern '[$(C=O)]' \
      --cnn_scoring none \
      -o covalent_docked.sdf.gz
```

| Option | Meaning |
|---|---|
| `--covalent_rec_atom` | The receptor atom, as `chain:resnum:atom_name` — or as `x,y,z` Cartesian coordinates |
| `--covalent_lig_atom_pattern` | A **SMARTS expression** matching the ligand atom that forms the bond |
| `--covalent_lig_atom_position` | Optional explicit `x,y,z` placement for the covalent ligand atom, instead of the OpenBabel heuristic |
| `--covalent_fix_lig_atom_position` | With the above, *fix* the atom at that position rather than only using it to build the initial structure |
| `--covalent_bond_order` | Bond order, default 1 |
| `--covalent_optimize_lig` | UFF-optimize the residue+ligand construct. Note this changes the ligand's bond lengths and angles |

### Every match is docked — keep the SMARTS sharp

If the SMARTS matches multiple ligand atoms, **all pairings of ligand and receptor atoms are
evaluated**, and the number of output poses expands correspondingly. `[GNINA1.3]`

This interacts with `--num_modes` in a way that is easy to get wrong. Suppose your SMARTS matches three
carbonyl carbons and you leave `--num_modes` at its default of 9. gnina generates poses for all three
attachment points, then sorts the *combined* pool, applies the 1 Å diversity filter, and truncates to 9.
Nothing guarantees those 9 are spread across the three attachment points — if one attachment scores well,
it can plausibly occupy most or all of the output, and you will never see the alternatives you thought
you were enumerating.

If you genuinely want to compare attachment points, either write a SMARTS that matches exactly one atom
and run once per candidate, or raise `--num_modes` well above the default and check the poses' geometry.
The same reasoning applies to `--min_rmsd_filter`: poses at *different* attachment points are far apart
in RMSD, so the filter will not conflate them, but poses at the same attachment point compete normally.

## Score it with Vina, not the CNN

> **Recommendation** `[GNINA1.3]` — run covalent docking with `--cnn_scoring none`.

gnina tells you this itself. Leave CNN scoring on with `--covalent_rec_atom` set and it warns:

```
WARNING: CNN scoring not yet calibrated for covalent docking.  Recommend running with --cnn_scoring none
```

The reason is straightforward domain-of-applicability: **the CNN models were never trained on covalent
complexes.** On the Scarpino 207-complex covalent redocking benchmark, Vina scoring performs
significantly better than the CNN. `[GNINA1.3]`

This is a nice illustration of a general point, and the paper draws it explicitly: a model that
outperforms a physics-based function inside its training domain can lose badly outside it. Note the
comparison cuts the other way on the *same* benchmark when covalent docking is switched off — CNN 27.5%
versus Vina 15.8%. It is specifically the covalent setting that the CNN has not seen.

## What accuracy to expect

Success rate on the Scarpino benchmark, measured as the fraction of targets whose top-ranked pose is
within 2 Å RMSD `[GNINA1.3]`:

| Setting | Success |
|---|---|
| Default: a generated ligand conformer, no positioning information | **36.2%** |
| Best case: experimental ligand conformer, covalent atom position specified | **66.6%** |
| Covalent docking *disabled*, CNN scoring | 27.5% |
| Covalent docking *disabled*, Vina scoring | 15.8% |

The spread between 36% and 67% is the value of prior information, and intermediate settings land in
between. The bottom two rows are the important comparison for deciding whether to bother: **enabling
covalent mode matters far more than which scoring function you use.**

With Vina scoring, gnina 1.3's covalent docking is competitive with, but does not outperform, the state
of the art. `[GNINA1.3]`

## Before 1.3: the soft-covalent hack

Earlier gnina and smina users faked covalent docking with a custom scoring term — a strongly attractive
gaussian between two atom types chosen so that only the intended pair could match. That never created a
bond; it just made the proximity overwhelmingly favourable. `[Workshop2021]`

It is still described in [scoring-empirical.md](scoring-empirical.md), and remains of some interest if
you want a soft *bias* toward proximity rather than an actual bond — but for real covalent docking, use
the options on this page.

## Under the hood

- Covalent setup, atom resolution and construct building:
  [covinfo.cpp](../../gninasrc/lib/covinfo.cpp), [covinfo.h](../../gninasrc/lib/covinfo.h)
- Option definitions and the CNN warning: `main.cpp`, in the covalent option group and the
  post-parse validation block
