# Output and ranking

**Audience:** users and developers.

## What gnina writes

The terminal table shows four columns per pose:

```
mode |  affinity  |  intramol  |    CNN     |   CNN
     | (kcal/mol) | (kcal/mol) | pose score | affinity
-----+------------+------------+------------+----------
    1       -6.76        0.57       0.8746      5.077
```

The output *file* carries more than that. Written to SDF (`result_info::write` in
[result_info.cpp](../../gninasrc/lib/result_info.cpp)):

| SD tag | Meaning | Written when |
|---|---|---|
| `minimizedAffinity` | The empirical score, kcal/mol. This is the `affinity` column | always |
| `minimizedRMSD` | RMSD from the input pose to the output pose | `--minimize` / `--local_only` |
| `CNNscore` | CNN pose score, 0–1 | CNN scoring enabled |
| `CNNaffinity` | CNN predicted affinity, pK | CNN scoring enabled |
| `CNN_VS` | **`CNNaffinity × CNNscore`** | CNN scoring enabled |
| `CNNaffinity_variance` | Spread across the ensemble | **ensembles only** — absent for single-model runs |
| `atomic_interaction_terms` | Per-atom interaction term values | `--atom_term_data` |

Three of these are worth calling out because they are easy to miss:

**`CNN_VS` is a built-in virtual-screening score.** It is simply the product of the pose score and the
predicted affinity, which downweights a confident-looking affinity attached to a dubious pose. It is
directly relevant to the "which metric should I rank by" question in
[virtual-screening.md](virtual-screening.md) — and note it is *not* an option for
`--pose_sort_order`, so using it means post-processing your output.

**`CNNaffinity_variance` only appears for ensembles.** A single-model run has no spread, so the field
is omitted rather than written as zero. If you are parsing output and expecting it unconditionally,
`--cnn fast` will surprise you. The terminal log also labels this `CNNvariance` while the SD tag is
`CNNaffinity_variance` — same quantity, two spellings.

**The empirical score is recomputed after CNN scoring**, deliberately, so that
`minimizedAffinity` is always a Vina-family energy even when the CNN did the ranking. The code comment
is explicit: "we want vina energies not CNN".

### Output format matters

Not all formats carry all fields.

- **`.sdf` and `.sdf.gz`** take a native fast path that writes every field above. This is the format to
  use, and the reason `.sdf.gz` is the standing recommendation — it is not only about disk space.
- **`.pdbqt`** takes its own native path, writing `REMARK minimizedAffinity`, `minimizedRMSD`,
  `CNNscore` and `CNNaffinity`. It **omits `CNN_VS` and the variance.**
- **Anything else** is converted by OpenBabel, which sets a reduced field set.

The format is chosen purely from the `-o` extension. An unrecognized extension is an error
(`Invalid format`), not a silent fallback.

## Ranking, then filtering

The order of operations here is the thing that surprises people, so it is worth stating exactly. After
the Monte Carlo chains finish and their pools are merged:

1. **Refine** each pose (locally minimize).
2. **Score** — compute CNN score and affinity, then recompute the empirical energy.
3. **Sort** the whole pool by `--pose_sort_order`.
4. **Filter for diversity** — walk the sorted list and drop any pose within `--min_rmsd_filter`
   (default 1 Å) RMSD of a better-ranked pose already kept.
5. **Truncate** to `--num_modes` (default 9) and write.

`--pose_sort_order` accepts `CNNscore` (default), `CNNaffinity`, or `Energy` — plus an undocumented
`vina` alias for `Energy`.

### Changing the sort order changes *which* poses you get

This is the subtlety worth internalizing: **sorting happens before filtering and truncation.** So
changing `--pose_sort_order` does not re-order the same nine output poses. It re-sorts the full internal
pool, which changes which poses survive the diversity filter, which changes which nine reach the
cutoff. A different sort order can produce an entirely different, even non-overlapping, set of output
poses.

The same reasoning explains something that otherwise looks like a bug: **switching CNN models changes
the output set, not just the order**, even at a fixed seed. A different CNN assigns different scores,
so a different sort order emerges, so different poses survive filtering. Comparing `--cnn fast`
against the default ensemble on a fixed seed will show this.

Two related notes:

- **`--cnn_scoring none` silently forces `Energy`** regardless of what you passed to
  `--pose_sort_order`. There is no CNN score to sort by, so this is correct — but it means the two flags
  are coupled.
- There is a **second, hardcoded 1 Å RMSD dedup** applied inside the Monte Carlo chains while poses
  accumulate, independent of `--min_rmsd_filter`. So the pool that reaches step 3 above is already
  somewhat diverse; `--min_rmsd_filter` controls only the final pass.

`[GNINA1.0]` measured `--min_rmsd_filter` at 0.5, 1.0 and 1.5 and found no significant effect on docking
performance — the CNN ranks the poses it sees accurately either way. Leave it at the default unless you
specifically want more or less pose diversity in your output.

## Getting per-atom detail

`--atom_terms <file>` writes per-atom interaction term values to a file; `--atom_term_data` embeds
them in the output SD data instead. Useful for understanding *why* a pose scored as it did, and for
building interpretable per-target models — see
[scoring-empirical.md](scoring-empirical.md).

For per-term (rather than per-atom) breakdowns, use `--score_only` with a custom scoring file, also
covered in [scoring-empirical.md](scoring-empirical.md).

## Under the hood

- Field writing and the per-format paths: `result_info::write` in
  [result_info.cpp](../../gninasrc/lib/result_info.cpp)
- Sorting, `remove_redundant`, and the `--num_modes` truncation: `main.cpp`, in the per-ligand output
  block
- The recognized SD tag list for round-tripping: [PDBQTUtilities.cpp](../../gninasrc/lib/PDBQTUtilities.cpp)
