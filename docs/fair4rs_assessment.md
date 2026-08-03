# FAIR4RS & Robustness Assessment

A review of how well `biobb_vs_workflows` meets **FAIR4RS** and
data-pipeline-engineering principles, plus a roadmap to close the gaps.

## Where the project stands today

Legend: 🟢 settled · 🟡 partial · 🔴 missing

### FAIR4RS

| Principle | Reality in the repo | Status |
|---|---|---|
| **Findable** | Public GitHub repo (`NBDsoftware/biobb_vs_workflows`), git tags (`0.0.1`–`0.0.2`), GitHub Pages docs, Colab tutorial. No `__version__` exposed or logged anywhere. No `CITATION.cff`. Zenodo DOI planned. | 🟡 |
| **Accessible** | Repo + tags + `LICENSE` present. License is **CC BY-NC-SA 4.0** — academic-open / commercially-restricted, with a separate commercial license available on request. | 🟢 |
| **Interoperable** | Standard formats in/out (PDB, SDF/SMILES ligand libs, YAML config), thin argparse CLIs, importable workflow functions. | 🟢 |
| **Reusable** | `--help` + README/docs exist, but reproducibility isn't in place yet — see Reproducibility below. | 🟡 |

### Robustness

| Principle | Reality in the repo | Status |
|---|---|---|
| **Reproducibility** | Python version and `openbabel` are pinned exactly; every other dependency (`biobb_*`, `mdanalysis`, `biopython`) is a loose `>=` lower bound. | 🟡 |
| **Provenance** | Each run writes `config.yml` (resolved paths baked in) and `log.out` + per-step logs — each BioBB step logs its own module version. No run manifest, input checksums, git commit hash, CLI argv, or `biobb_vs_workflows` package version are not recorded. | 🟡 |
| **Environment portability** | Conda `environment.yml` exists. No container image yet. | 🟡 |
| **Modularity** | Genuinely modular at the BioBB-step level. The only shared code between the two workflows is `common.to_yaml()`; `create_config_file()` and the config-scaffolding pattern are byte-identical copy-paste between workflows. | 🟡 |
| **Validation** | Step outputs are validated **internally by the biobbs**. Each workflow has its own ad-hoc `check_arguments()` (existence checks, mutual exclusivity, `box_offset` bounds), but no shared validation module, and argparse `choices=` only on the `virtual_screening` engine flags. Two spots log an error then silently `return`: the unsupported-ligand-library-extension check in `virtual_screening.py` and `get_clusters_population()` in `cavity_analysis.py`. | 🟡 |

## Roadmap

### 1. Improve clarity of errors
- Turn the two silent `log.error(...); return` exits into `raise`/`sys.exit`: the unsupported-library-extension branch in `virtual_screening.py` and `get_clusters_population()` in `cavity_analysis.py`.

### 2. Version exposure & provenance
- Expose `__version__` (via `importlib.metadata`) and log it at the start of every run.
- One shared helper writing `output/run_manifest.json` — git commit (if
  resolvable), full `sys.argv`, SHA-256 of each input file, UTC timestamp, resolved
  `config.yml` path, and a `conda env export` / `pip freeze` snapshot.

### 3. Validation
- argparse `choices=` for enum-like flags (e.g. ligand-library format).
- Extract the two workflows' duplicated `check_arguments()` logic into `common/`, alongside the existing `to_yaml()` helper.

### 4. Testing & automation
- Extend CI beyond the existing docs-build workflow: add `ruff` lint + `pip install .`.
- One tiny end-to-end smoke test on a small system (reuse `tests/*/run.sl`'s logic against the `data/` fixtures, wired into CI instead of SLURM-only).
- Dependency automation (Dependabot/Renovate) watching the loose `biobb_*` version floors.

### 5. FAIR4RS metadata & portability
- Create a `CITATION.cff` (doesn't exist yet).
- Mint a **Zenodo DOI** for easy academic citation: (1) org admin
  authorizes Zenodo for the NBDsoftware org, (2) enable the repo in Zenodo, (3) cut a GitHub Release
  → DOI auto-minted, (4) paste the concept-DOI into `CITATION.cff` and a README badge.
- License: **settled** for now (CC BY-NC-SA 4.0 dual-license) — worth revisiting whether a more code-conventional license fits better long-term.
- Build an Apptainer/Singularity or Docker image per release for one-command portable execution.
