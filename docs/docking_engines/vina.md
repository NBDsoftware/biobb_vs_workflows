# AutoDock Vina

AutoDock Vina is a molecular docking program: given a receptor and a ligand,
it searches for binding poses and scores them with an empirical function.
It's the default docking engine for the `virtual_screening` workflow.

## How it works

- **Search.** Vina generates candidate poses with an Iterated Local Search —
  a Monte Carlo-like global search where each step is refined by a local
  BFGS optimization. `--exhaustiveness` controls how many independent
  searches run; more searches trade speed for a better chance of finding the
  true best pose.
- **Scoring.** Poses are ranked by a weighted sum of pairwise atomic terms —
  two Gaussian attraction terms, a repulsion term, a hydrophobic term, and a
  directional hydrogen-bond term — plus a penalty proportional to the
  ligand's rotatable bonds. The weights were fit against the PDBbind dataset
  of experimental protein-ligand structures and affinities. There's no
  explicit electrostatics or solvation term.
- **Score.** The result is a predicted binding affinity in kcal/mol —
  **lower (more negative) is better**. It's a relative ranking score, not a
  physical binding free energy.

## Practical notes

- The scoring function only looks at atoms within an 8 Å cutoff, and treats
  the receptor as rigid (only the ligand's torsions are optimized).
- Vina is deterministic given a fixed random seed; without one, scores and
  poses vary slightly between runs.

## Reference

Trott O, Olson AJ. *AutoDock Vina: improving the speed and accuracy of
docking with a new scoring function, efficient optimization, and
multithreading.* J Comput Chem. 2010;31(2):455-461.
[doi:10.1002/jcc.21334](https://doi.org/10.1002/jcc.21334)

For the CLI flags exposed by this repo (`--vina_bin`, `--exhaustiveness`,
`--cpus`, ...), see the [virtual screening workflow](../workflows/virtual_screening.md).
