# Installing gnina

**Audience:** everyone. [README.md](../../README.md#installation) has the exact commands; this page
covers the reasoning behind the choices, and the one situation the README does not address (no
NVIDIA GPU).

gnina is Linux-only. There is no NVIDIA GPU support on macOS, so on a Mac the practical option is
smina instead — CPU-only, empirical scoring, no CNN.

## Pick a route

> **Recommendation** — use a pre-built binary or a Docker image unless you specifically need a build
> tuned to your own system. `[Workshop2021]`

| Route | When |
|---|---|
| [GitHub release binary](https://github.com/gnina/gnina/releases/latest) | Default choice. Works under WSL2. Large download; see below. |
| [Docker image](https://hub.docker.com/u/gnina) | You already run containers, or you need reproducibility. Dockerfiles in [docker/](../../docker/): `ubuntu-22.04`, `ubuntu-24.04`, `compat`. Note the images are large — tens of GB. |
| From source | You need a specific toolchain, you are developing gnina, or you want a small binary. See [developing.md](developing.md). |
| Binary + CUDA wheels | You have no NVIDIA GPU at all. See [below](#running-with-no-nvidia-gpu). |

**Why the release binary is so large.** Rather than requiring users to install matching versions of
CUDA, cuDNN, libtorch, OpenBabel and Boost, the release binaries statically bundle almost everything
except the barest minimum of system libraries, built against old-enough baselines to run on any
reasonably current Linux. CUDA accounts for most of the size. A from-source build linked dynamically
against your own system libraries is roughly 40 MB; the bundled binary is around half a gigabyte.
`[Workshop2021]`

[docker/compat/Dockerfile](../../docker/compat/Dockerfile) is the maximum-compatibility build that
produces this: it builds OpenBabel, MKL-DNN (oneDNN), PyTorch and libmolgrid from source, then hand-writes
a link line producing `gnina.static`. It is a *semi*-static build — the torch, NCCL and c10 static
archives are `--whole-archive`d, while cuDNN and most CUDA libraries remain ordinary dynamic links.
That distinction matters for the no-GPU case below.

Release assets arrive with mode `644`. You need `chmod +x` before running one.

## GPU requirements

An NVIDIA GPU is required for *fast* CNN scoring, not for correctness — everything works on CPU, just
slower. At least 4 GB of VRAM is recommended, more is better.

- `--no_gpu` forces CPU-only execution even when a GPU is present.
- `--device` is **ignored by the PyTorch backend.** If you pass a non-zero value gnina warns and tells
  you to use `CUDA_VISIBLE_DEVICES` instead. The flag still feeds the experimental non-CNN GPU path
  (see [cli-reference.md](cli-reference.md)), but for selecting which GPU does CNN scoring,
  `CUDA_VISIBLE_DEVICES` is the mechanism.
- With no GPU detected and `--no_gpu` not given, gnina warns and continues:

  ```
  WARNING: No GPU detected. CNN scoring will be slow.
  Recommend running with single model (--cnn fast)
  or without cnn scoring (--cnn_scoring=none).
  ```

  That advice is worth taking; see [performance.md](performance.md) for the numbers behind it.

## Building from source, in brief

Full commands are in [README.md](../../README.md#installation). The shape of it:

- System packages via apt: build-essential, git, wget, cmake, Boost (all components), Eigen3, glog,
  protobuf, HDF5, ATLAS, RDKit, jsoncpp, Python 3 with numpy/pytest/pip
- CUDA ≥ 12.0, following NVIDIA's own install instructions, with `nvcc` on `PATH`
- OpenBabel 3 built from source. Versions ≤ 3.1.1 have bond-order-determination bugs; either the
  `dkoes/openbabel` fork or current upstream `openbabel/openbabel` works
- `cmake .. && make && make install` from a `build/` directory inside the checkout

CMake auto-fetches **libtorch** (via `FetchContent`) and **libmolgrid** (via `ExternalProject_Add`) if
it does not find them installed — neither is a git submodule, and `.gitmodules` is empty. Set
`GNINA_FORCE_EXTERNAL_LIBS` to require pre-installed copies instead of fetching. RDKit is optional and
only needed for the `gninavis` visualization tool; the build skips it with a warning if absent.

**One caveat about GPU architectures.** README suggests `-DCMAKE_CUDA_ARCHITECTURES=all` when building
for a cluster with heterogeneous GPUs, but [CMakeLists.txt](../../CMakeLists.txt) unconditionally does
a non-cache `set(CMAKE_CUDA_ARCHITECTURES "all-major")`, which overrides that command-line value. If
you genuinely need `all`, edit the `set()` rather than passing `-D`. See
[developing.md](developing.md) for the rest of the build internals.

## Smoke test

```bash
gnina --version
gnina -r test/gnina/data/184l_rec.pdb -l test/gnina/data/184l_lig.sdf \
      --autobox_ligand test/gnina/data/184l_lig.sdf --seed 0 -o /tmp/out.sdf
```

You should get a table of nine poses with affinity, CNN pose score and CNN affinity columns, and a
`/tmp/out.sdf` containing them. If the scores look plausible but `/tmp/out.sdf` is missing, you
forgot `-o` — see [running-and-sampling.md](running-and-sampling.md).

## Running with no NVIDIA GPU

This is the case the README does not cover, and the failure mode is confusing enough to be worth
spelling out.

### Compute and linkage are two different problems

`--no_gpu` solves the first one only.

**Compute is fine without a GPU.** `--no_gpu` forces the CPU path, and gnina degrades rather than
aborting even without the flag: it only warns, and `initializeCUDA` in
[dl_scorer.cpp](../../gninasrc/lib/dl_scorer.cpp) deliberately returns silently when `cudaSetDevice`
fails — the code comment is literally "be silent if GPU not present". gnina's own test suite runs
with `--no_gpu` throughout `test/gnina/test_cnn.py`.

**Linkage is the actual blocker.** The release binaries are *not* self-contained. They `DT_NEEDED`
seven CUDA/cuDNN userspace libraries — on a current asset:

```
libcudnn.so.9  libcudart.so.12  libcublas.so.12  libcublasLt.so.12
libcusolver.so.11  libcufft.so.11  libcusparse.so.12
```

Without them the process dies at `ld.so` time, before `main()`, no matter which flags you pass.
`readelf -d` also shows `FLAGS_1: NOW PIE` — the binary is linked **`BIND_NOW`**, so every undefined
symbol resolves eagerly at load. **Empty stub `.so` files therefore do not work**; the real libraries
are required.

No *driver* is needed. `libcuda.so.1` is absent from the `NEEDED` list, because `libcudart` dlopens
the driver lazily and on a CPU-only run never gets that far.

Check your own asset rather than trusting the list above, since it varies between releases:

```bash
readelf -d ./gnina.bin | grep NEEDED
```

Older assets additionally need `libnvToolsExt.so.1`.

### Why not Docker or a source build

Neither is a reasonable route on a GPU-less workstation. The official images are tens of GB. A source
build is impossible without the CUDA toolkit regardless of whether a GPU exists:
[CMakeLists.txt](../../CMakeLists.txt) declares `project(gnina C CXX CUDA)` on line 2 and does
`find_package(CUDA 12.0 REQUIRED)`, and there is **no `CPU_ONLY` option anywhere in the project** —
the only CMake options are `BUILD_COVERAGE` and `GNINA_FORCE_EXTERNAL_LIBS`. The `CPU_ONLY` mentions
in [docs/install_apt_debian.md](../install_apt_debian.md) are leftovers from upstream Caffe's
documentation and do not apply to gnina.

So: take a release binary and supply the CUDA libraries from PyPI wheels.

### The recipe

```bash
PREFIX=$HOME/opt/gnina          # wherever you want this to live
mkdir -p "$PREFIX"
# ... download the release asset to $PREFIX/gnina.bin, then:
chmod +x "$PREFIX/gnina.bin"

# 1. See what your asset actually needs.
readelf -d "$PREFIX/gnina.bin" | grep NEEDED

# 2. Install matching NVIDIA wheels into a throwaway venv.
python3 -m venv "$PREFIX/cudaenv"
"$PREFIX/cudaenv/bin/pip" install --no-cache-dir \
  nvidia-cuda-runtime-cu12 nvidia-cublas-cu12 nvidia-cudnn-cu12 \
  nvidia-cufft-cu12 nvidia-cusolver-cu12 nvidia-cusparse-cu12 \
  nvidia-nvjitlink-cu12

# 3. Collect every shipped .so into one directory.
mkdir -p "$PREFIX/cudalibs"
find "$PREFIX/cudaenv" -path '*/nvidia/*/lib/*.so*' \
  -exec ln -sf {} "$PREFIX/cudalibs/" \;

# 4. Confirm nothing is missing. This must print nothing.
LD_LIBRARY_PATH="$PREFIX/cudalibs" ldd "$PREFIX/gnina.bin" | grep 'not found'
```

Notes on the pieces:

- **Pin to the CUDA minor line the binary was built against.** The wheel names above are unpinned for
  readability; in practice pin them (e.g. `nvidia-cudnn-cu12==9.8.0.87`) so a later wheel release does
  not silently change the ABI under you. The `--version` output and the release notes tell you which
  CUDA line an asset targets.
- **`nvidia-nvjitlink-cu12` is not optional**, even though nothing in gnina's `NEEDED` list mentions
  it: the wheel builds of `libcusolver` and `libcusparse` `DT_NEEDED` it.
- `--no-cache-dir` avoids leaving a second copy of the ~2 GB download in `~/.cache/pip`.

Then wrap it, so `LD_LIBRARY_PATH` and `--no_gpu` are not your problem on every invocation:

```bash
#!/usr/bin/env bash
export LD_LIBRARY_PATH="$PREFIX/cudalibs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
case " $* " in
  *" --no_gpu "*) exec "$PREFIX/gnina.bin" "$@" ;;
  *)              exec "$PREFIX/gnina.bin" --no_gpu "$@" ;;
esac
```

The `case` matters: `--no_gpu` is a Boost.Program_options switch, and passing it twice is an error
("multiple occurrences"), so inject it only when it is absent.

A worked instance of this, with concrete paths and pinned versions, is in
[local/cpu-only-workstation.md](local/cpu-only-workstation.md).

### Using gnina on CPU

Everything works, including CNN scoring; the cost is wall-clock time. Cheapest first:

- `--cnn_scoring none` for large batches where you will rescore a shortlist later
- `--cnn fast` as the middle ground when you want a CNN score per pose — this is exactly the case it
  was built for `[GNINA1.3]`
- the default 3-model ensemble for one-off runs

Avoid `--cnn_scoring refinement` (roughly 10× `rescore`) and `all` entirely on CPU. Set `--cpu`
explicitly rather than letting gnina claim every core. Measured timings are in
[performance.md](performance.md).

One thing to expect: pose *ranking* differs between `--cnn fast` and the default ensemble even at a
fixed seed. That is correct behaviour, not a bug — `--pose_sort_order` defaults to `CNNscore` and
re-sorts the pool *before* the diversity filter, so a different CNN surfaces a different final set of
poses, not merely the same set in a different order. See
[output-and-ranking.md](output-and-ranking.md).

Do **not** try `make test` or `ctest` on a setup like this: [test/CMakeLists.txt](../../test/CMakeLists.txt)
needs a configured CUDA build tree, which is exactly what this route avoids.
