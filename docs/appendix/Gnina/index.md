# The gnina guide

**Audience:** everyone. Start here.

gnina (pronounced "NEE-na") is a molecular docking program. Give it a receptor structure, a ligand,
and a definition of where to dock, and it samples ligand poses and scores them — aiming both to
predict the correct binding pose and to rank compounds by likely affinity. It adds
convolutional-neural-network (CNN) scoring and GPU acceleration on top of a sampling engine
inherited from AutoDock Vina by way of smina.

This guide is the "why" and "how" layer. [README.md](../../README.md) owns installation one-liners,
the raw generated `--help` output, and the citation list; this directory owns the reasoning,
the scientific context, and the behaviour that `--help` does not tell you. No fact is deliberately
duplicated between the two except a single quickstart command.

## Reading paths

**First time with gnina** — result before theory, because docking output is hard to reason about in
the abstract:

1. [install.md](install.md) — get a working binary
2. [running-and-sampling.md](running-and-sampling.md) § Quickstart — one command, one output file
3. [concepts.md](concepts.md) — what the pipeline actually does, and what a score is not
4. [inputs-and-search-space.md](inputs-and-search-space.md) — the part that most often goes wrong
5. [output-and-ranking.md](output-and-ranking.md) — read what you just produced
6. [scoring-cnn.md](scoring-cnn.md) § How to read a CNNscore
7. [gotchas.md](gotchas.md) — skim once, re-read after your first surprise

**"Should I use gnina, and can I trust this number?"**
[performance.md](performance.md) → [concepts.md](concepts.md) →
[scoring-cnn.md](scoring-cnn.md) § How to read a CNNscore.

**"I'm screening a library"**
[virtual-screening.md](virtual-screening.md) → [performance.md](performance.md) →
[output-and-ranking.md](output-and-ranking.md) → [scoring-cnn.md](scoring-cnn.md) (for `--cnn fast`).

**Developers** — start at [developing.md](developing.md), which opens with the code map and the two
main call chains before anything else.

## Contents

| Page | What's in it |
|---|---|
| [concepts.md](concepts.md) | Vina → smina → gnina lineage; why a score is not an energy; the six-stage pipeline; redocking vs cross-docking vs whole-protein; what 1.0 and 1.3 each changed |
| [install.md](install.md) | Release binary, Docker, source build, GPU requirements, smoke test, and running with no NVIDIA GPU |
| [inputs-and-search-space.md](inputs-and-search-space.md) | Receptor and ligand preparation, protonation, defining the box, `--autobox_add`, whole-protein docking |
| [running-and-sampling.md](running-and-sampling.md) | Quickstart, the run modes, Monte Carlo internals, the tuning knobs that matter, validating with `obrms` |
| [scoring-empirical.md](scoring-empirical.md) | The term library, how Vina/Vinardo/AD4 actually differ, custom scoring functions, target-specific models |
| [scoring-cnn.md](scoring-cnn.md) | The 3D-grid representation, model families, ensembles, knowledge distillation, `--cnn_scoring` modes, interpreting the scores |
| [output-and-ranking.md](output-and-ranking.md) | Every field gnina writes, per-format differences, and how sorting interacts with the diversity filter |
| [performance.md](performance.md) | Published accuracy and runtime, by version, task and scoring choice |
| [flexible-docking.md](flexible-docking.md) | Flexible side chains, why rigid is the right default, reassembling output |
| [covalent-docking.md](covalent-docking.md) | Covalent docking: input requirements, atom selection, and why to score it with Vina |
| [virtual-screening.md](virtual-screening.md) | Screening at scale, which metric to rank by, Pharmit prefiltering, the `deepdock.py` pipeline |
| [cli-reference.md](cli-reference.md) | Every option, grouped, with real defaults and behavioural footguns |
| [gotchas.md](gotchas.md) | One-line cheat sheet, each item linking to the page that owns it |
| [developing.md](developing.md) | Code map, call chains, what's dead, build internals, adding a CNN model, tests and CI |

## Conventions

Reading this guide is easier if you know which claims come from where, because the three kinds of
source carry very different weight.

**Unmarked statements describe the current code.** They were checked against this repository and the
installed binary. If one is wrong, the code is the authority and the guide has a bug.

**Published results carry a trailing citation key.** For example: cross-docking Top1 improves from
27% to 37% `[GNINA1.0]`. Keys resolve in the table below. A number without a key is either measured
locally (and labelled as such) or a bug.

**Opinions and recommendations are marked as such**, so you can tell a benchmark from a preference:

> **Recommendation** `[GNINA1.0]` — use rigid-receptor docking by default; reach for flexible side
> chains only with specific structural evidence.

**Things that changed** get a one-line note rather than being silently overwritten, because plenty of
gnina material on the web predates the current release:

> *Changed since 1.0:* the default CNN ensemble is 3 models, not 5.

**Source references.** Files are linked by identity — [everything.h](../../gninasrc/lib/everything.h).
Specific lines appear as unlinked text plus the symbol name — "the `--no_gpu` handling in
`main.cpp:1057`" — because line numbers rot within weeks in an active repository while symbol names
stay greppable.

### Sources

| Key | Source |
|---|---|
| `[GNINA1.0]` | McNutt et al., *GNINA 1.0: molecular docking with deep learning*, J Cheminform 13:43 (2021). Reprint: [sources/GNINA_1/](sources/GNINA_1/) |
| `[GNINA1.3]` | McNutt et al., *GNINA 1.3: the next increment in molecular docking with deep learning*, J Cheminform 17:28 (2025). Reprint: [sources/GNINA_1.3/](sources/GNINA_1.3/) |
| `[VS1.0]` | Sunseri & Koes, *Virtual screening with Gnina 1.0*, Molecules 26:7369 (2021) |
| `[CNN2017]` | Ragoza et al., *Protein–ligand scoring with convolutional neural networks*, J Chem Inf Model 57:942 (2017) — the primary methods citation for CNN scoring |
| `[CrossDocked]` | Francoeur et al., *Three-dimensional convolutional neural networks and a cross-docked data set for structure-based drug design*, J Chem Inf Model 60:4200 (2020) — the CrossDocked2020 dataset |
| `[KD2024]` | McNutt, Li, Francoeur, Koes, *Condensing molecular docking CNNs via knowledge distillation*, ChemRxiv (2024) — training detail behind the distilled models |
| `[CACHE1]` | Dunn et al., *CACHE Challenge #1: Docking with GNINA is all you need*, J Chem Inf Model (2024) — a prospective application |
| `[Workshop2021]` | CECAM / Liverpool ChiroChem open-source-software workshop, David Koes presenting, with Q&A from Drew McNutt, Rocco Meli and others. Video and slides linked from [README.md](../../README.md#help); notebook and slides in [docs/rsc_workshop2021/](../rsc_workshop2021/) |

`[Workshop2021]` is a recorded talk, not a peer-reviewed result. It is the only source for several
practical recommendations in this guide, and it is now several releases old — where a paper has since
published the same comparison, the guide cites the paper.

## Getting help

- Slack: invite link in [README.md](../../README.md#help)
- GitHub issues: bugs and feature requests
- Example Colab notebook: linked from [README.md](../../README.md#help)

## Verified against

- Repository: branch `master`, commit `6fe1ce2b`
- Binary: `gnina v1.3.3 master:6fe1ce2`, checked via `--help`, `--help_hidden`, `--score_only` and a
  full docking run
- Papers: the two reprints in [sources/](sources/)

When you change behaviour that this guide describes, update the page and this line together.
[cli-reference.md](cli-reference.md) and the model list in [scoring-cnn.md](scoring-cnn.md) are the
two places that go stale fastest; both say how to regenerate themselves.
