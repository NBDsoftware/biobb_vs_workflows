# Developing gnina

**Audience:** developers.

## Start here

Four things before the repository layout table, because they are what actually orients you.

**1. `main()` and every CLI flag are in one file.**
[gninasrc/main/main.cpp](../../gninasrc/main/main.cpp) defines the Boost.Program_options option
groups and drives the whole pipeline from parsing through output. It is the map to every flag in
[cli-reference.md](cli-reference.md) — if a flag's behaviour is unclear, this is where to look first.

**2. Two call chains matter.**

- **Scoring:** `main.cpp` → [cnn_torch_scorer.cpp](../../gninasrc/lib/cnn_torch_scorer.cpp) (model
  selection and ensemble resolution) → [torch_model.cpp](../../gninasrc/lib/torch_model.cpp) (grid
  construction, metadata parsing) → [dl_scorer.cpp](../../gninasrc/lib/dl_scorer.cpp) (device setup).
- **Sampling:** [monte_carlo.cpp](../../gninasrc/lib/monte_carlo.cpp) →
  [mutate.cpp](../../gninasrc/lib/mutate.cpp) (the move set) →
  [quasi_newton.cpp](../../gninasrc/lib/quasi_newton.cpp) /
  [bfgs.cu](../../gninasrc/lib/bfgs.cu) (minimization).

**3. What's live versus dead.**
[cnn_data.cpp](../../gninasrc/lib/cnn_data.cpp) looks like the CNN model registry and is not — it is
not in `LIB_SRCS`, nothing references `builtin_cnn_models` or `default_model_name`, and its 22-model
list and `crossdock_default2018` default are Caffe-era leftovers. The real registry is generated at
build time; see [Adding a CNN model](#adding-a-cnn-model) below.
[covinfo.cpp](../../gninasrc/lib/covinfo.cpp) / `.h` is live and is where covalent docking setup lives.
Also obsolete: [docs/tutorial/layers/](../tutorial/layers/) (Caffe layer reference, 59 files) and
[include/caffe/util/nccl.hpp](../../include/caffe/util/nccl.hpp), the only file under `include/`.

**4. Build, test, then the model-bundling contract**, both covered below.

## Building from source

```bash
mkdir build && cd build && cmake .. && make -j8 && make install
```

Requirements and the reasoning behind them are in [install.md](install.md). Two build-system details
worth knowing that are not obvious from the README:

- **`project(gnina C CXX CUDA)` is at [CMakeLists.txt:2](../../CMakeLists.txt)**, not further down where
  you might expect a language declaration. `find_package(CUDA 12.0 REQUIRED)` is around line 38. There
  is no `CPU_ONLY` option anywhere in the project.
- **`CMAKE_CUDA_ARCHITECTURES` is hard-set.** `CMakeLists.txt` does a non-cache
  `set(CMAKE_CUDA_ARCHITECTURES "all-major")`, which **overrides** a `-DCMAKE_CUDA_ARCHITECTURES=all`
  passed on the command line — exactly the value the README suggests for heterogeneous clusters. If you
  need `all` rather than `all-major`, edit the `set()` call. There is also a `CUDA_VERSION_MAJOR >= 13`
  branch selecting `TORCH_CUDA_ARCH_LIST`.
- **`GNINA_FORCE_EXTERNAL_LIBS`** disables the libtorch/libmolgrid auto-fetch (`FetchContent` /
  `ExternalProject_Add`) and requires pre-installed copies instead. Useful in a controlled build
  environment where fetching from `download.pytorch.org` at configure time is undesirable.
- `docker/compat/Dockerfile` produces the maximum-compatibility semi-static binary. It is *semi*-static:
  only the torch, NCCL and c10 static archives are linked `--whole-archive`; cuDNN and most CUDA
  libraries (`cudart`, `cusparse`, `cufft`, `cublas`, `cublasLt`, `cusolver`, `curand`, `nvToolsExt`)
  remain ordinary dynamic links. See [install.md](install.md) for what that implies for CPU-only hosts.

## Adding a CNN model

Custom TorchScript models are loaded with `--cnn_model <file>`, which can be combined with built-in
`--cnn` names. To package one, `torch.jit.trace` it and attach JSON metadata describing the input grid
and atom typing — the full recipe is in
[gninasrc/lib/models/README.md](../../gninasrc/lib/models/README.md):

```python
d = {
    'resolution': 0.5,
    'dimension': 23.5,
    'recmap': '...',   # newline/space-separated receptor atom type names
    'ligmap': '...',   # same, for the ligand
}
extra = {'metadata': json.dumps(d)}
z = torch.zeros((1, 28, 48, 48, 48))
script = torch.jit.trace(model, z)
script.save('my_model.pt', _extra_files=extra)
```

Load it with `--cnn_model my_model.pt`. See [scoring-cnn.md](scoring-cnn.md) for what `resolution` and
`dimension` mean and why they are per-model rather than global.

**Bundling a model into the binary itself** (rather than loading it at runtime) means adding it to the
`.pt` list in [gninasrc/CMakeLists.txt](../../gninasrc/CMakeLists.txt). At configure time,
[gninasrc/make_model_cpp.py](../../gninasrc/make_model_cpp.py) generates `torch_models.cpp`, embedding
each listed `.pt` file and deriving its `--cnn` name from the filename (`.` → `_`). This script, plus
the CMake list, **is** the model registry — not `cnn_data.cpp`.

## Training new CNN models

- Training scripts and pretrained PyTorch models: [gnina-torch](https://github.com/RMeli/gnina-torch)
- Legacy Caffe training scripts: [gnina/scripts](https://github.com/gnina/scripts)
- Sample pretrained models: [gnina/models](https://github.com/gnina/models)
- Training dataset: [CrossDocked2020](https://github.com/gnina/models/tree/master/data/CrossDocked2020)
  `[CrossDocked]`. The DUD-E docked poses used in the original gnina paper are also available, but
  training virtual-screening models directly on DUD-E is
  [not recommended](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0220113) due to
  known dataset biases — see [virtual-screening.md](virtual-screening.md).

Transfer learning — starting from a pretrained model and fine-tuning on a smaller target-specific
dataset — is reasonable to try and has been published successfully by others, but is not guaranteed:
one documented attempt at fine-tuning toward a Factor Xa-like dataset did not work well
`[Workshop2021]`. If you try it, build in a way to actually evaluate whether it helped before trusting
it. Knowledge distillation (see [scoring-cnn.md](scoring-cnn.md)) is the training technique behind the
`fast` and `*_KD*` models, if you want a worked example of the loss formulation `[KD2024]`.

The hardest part of any of this is getting your data into the model in the first place: training
scripts expect a simple text manifest, one labelled receptor/ligand example per line.

## Repository layout

| Path | What's there |
|---|---|
| `gninasrc/main/` | `main.cpp` — the executable entry point and every CLI option definition |
| `gninasrc/lib/` | Core docking library: Monte Carlo (`monte_carlo.cpp`), minimization (`quasi_newton.cpp`, `bfgs.cu`), empirical terms (`everything.h`, `builtinscoring.cpp`, `custom_terms.cpp`), CNN/torch integration (`cnn_torch_scorer.cpp`, `torch_model.cpp`, `dl_scorer.cpp`), covalent docking (`covinfo.cpp/.h`), option struct (`user_opts.cpp/.h`), GPU kernels (`grid_gpu.cu`, `model.cu`, `conf_gpu.cu`, `precalculate_gpu.cu`), molecule/PDBQT I/O, `models/` (bundled `.pt` weights + typing README), vendored LLVM-style option parser (`CommandLine2/`) |
| `gninasrc/gninagrid/` | Standalone grid-generation tool |
| `gninasrc/gninatyper/` | Atom-typing utility |
| `gninasrc/gninavis/` | CNN visualization tool (needs RDKit) |
| `gninasrc/gninaserver/` | Query server, with a Python `client.py` |
| `gninasrc/fromgnina/`, `gninasrc/tognina/` | Format conversion utilities |
| `gninasrc/pygnina/` | Boost.Python bindings (`bindings.cpp`, `setup.py`) — see [status](#python-bindings-pygnina) below |
| `gninasrc/make_model_cpp.py` | Generates `torch_models.cpp`, the real CNN model registry, at build time |
| `docker/` | `ubuntu-22.04/`, `ubuntu-24.04/` (standard builds), `compat/` (maximum-compatibility semi-static build) |
| `docs/` | This guide (`source/`), `rsc_workshop2021/` (workshop slides/notebook), `tutorial/` (obsolete Caffe layer reference), `install_apt_debian.md` |
| `examples/` | `kitchensink.score` — a worked custom scoring function |
| `scripts/` | `makeflex.py` (reassemble a receptor after flexible docking), `deepdock.py` (iterative screening pipeline — see [virtual-screening.md](virtual-screening.md)), `slurm.yaml`, `makemodel.ipynb`, `split_caffe_proto.py` (Caffe-era), `README.md` (covers `makeflex.py` only) |
| `test/` | CTest-driven suite: `test/gnina/`, `test/gninagrid/`, `test/pygnina/` |
| `Eigen/` | Vendored, patched Eigen 3.2 snapshot |
| `cmake/` | `version.cmake` (embeds git tag/rev/branch into the binary, exposed via `--version`), `Modules/` (Find scripts for libmolgrid, NCCL, RDKit) |
| `LICENSE.APACHE`, `LICENSE.GNU` | Dual license texts |
| `CITATION.cff` | Machine-readable citation metadata |

Two directories exist locally but are **untracked** — `docs/source/` (this guide) and `gnina_bins/`
(downloaded/built release binaries, several GB). Be careful with `git add -A` around either.

## Python bindings (pygnina)

[gninasrc/pygnina/](../../gninasrc/pygnina/) provides **Boost.Python** bindings (not pybind11),
exposing a `GNINA` class and a `result_info` class, built via its own `CMakeLists.txt`/`setup.py`.

> **Status: incomplete.** Only `set_receptor` and `score` are exposed. `GNINA::minimize` exists in C++,
> but its `std::istream` overload is a literal `abort()` and neither overload is bound with `.def(...)`.
> `result_info` also carries `cnnvariance` internally, but only `energy`, `cnnscore`, `cnnaffinity` and
> `write` are exposed. If you need minimize-from-Python or ensemble variance, they are not there yet —
> either shell out to the `gnina` binary or extend the bindings.

Working API, matching [test/pygnina/test_score.py](../../test/pygnina/test_score.py):

```python
from pygnina import GNINA
g = GNINA()
g.set_receptor(receptor_path_or_string, "pdb")
result = g.score(ligand_path_or_string, "sdf")
result.energy(); result.cnnscore(); result.cnnaffinity()
```

Both `set_receptor` and `score` accept either a path/string or a Python file-like object.

## Other executables

- **`gninagrid`** — standalone grid generation.
- **`gninatyper`** — atom-typing utility.
- **`gninavis`** — CNN visualization; requires RDKit, skipped gracefully if absent.
- **`fromgnina` / `tognina`** — format conversion utilities.
- **`gninaserver`** — a query server with a bundled Python client.
- **`scripts/makeflex.py`** — see [flexible-docking.md](flexible-docking.md).
- **`scripts/deepdock.py`** — see [virtual-screening.md](virtual-screening.md).

## Testing and CI

The suite is CTest-driven (`add_subdirectory(test)` in the top-level `CMakeLists.txt`):
`test/gnina/` (`test_gnina.py`, `test_cnn.py`, `test_flex.py`, `test_min.py`, plus
`correctness.py`/`speed.py` and C++/CUDA unit tests), `test/gninagrid/` (grid-generation fixtures),
`test/pygnina/` (Python bindings).

Two GitHub Actions workflows run on a **self-hosted GPU runner**, on every push/PR to `master`:

- **`CI.yml`** (also weekly, `cron: 0 0 * * 0`): `cmake .. && make -j8 && make install`, then
  `export CTEST_OUTPUT_ON_FAILURE=1 && make test`. **Not `ctest` directly** — it goes through the
  Makefile target.
- **`Coverage.yml`**: the same build with `-DBUILD_COVERAGE=1`, a hardcoded
  `PYTHONPATH=.../build/external/lib/python3.8/site-packages/`, then
  `ctest --timeout 7200 && ctest -T Coverage` under `continue-on-error: true`, uploaded to Codecov.

**Do not run `make test`/`ctest` on a CPU-only or unconfigured host** —
[test/CMakeLists.txt](../../test/CMakeLists.txt) needs a configured CUDA build tree.
