# Site-specific example: gnina on a CPU-only laptop

**Not maintained, and not part of the guide.** This is one worked instance of the general recipe in
[../install.md § Running with no NVIDIA GPU](../install.md#running-with-no-nvidia-gpu), kept so the
exact versions do not have to be rediscovered. Read the general recipe first; use this only as a
reference for what a working setup looks like.

Host: Ubuntu 24.04, Intel Iris Xe graphics only, no NVIDIA GPU, no driver, no CUDA toolkit.

## Layout

```
~/opt/gnina/gnina.bin   -> gnina_bins/gnina.cuda12.8.static
~/opt/gnina/cudaenv/    venv holding only the nvidia-*-cu12 wheels
~/opt/gnina/cudalibs/   symlinks to every .so under cudaenv/.../nvidia/*/lib
~/.local/bin/gnina      launcher: sets LD_LIBRARY_PATH, injects --no_gpu, execs gnina.bin
```

## Which binary

`gnina_bins/gnina.cuda12.8.static` is preferred over the older `gnina.1.3.2` asset alongside it: it
reports `v1.3.3` and needs one fewer library (no `libnvToolsExt.so.1`).

Note on provenance: this binary reports `master:6fe1ce2 Built Jun 30 2026`, and `6fe1ce2` is a commit
on *this fork's* master (Jun 29 2026) — so it was built from this repository, not downloaded from
upstream gnina releases. That also means it reflects this fork's current code, which is why it is
usable for checking documented behaviour.

`gnina_bins/` itself is ~3.4 GB of untracked binaries living inside the repository working tree. It is
not in `.gitignore`; be careful with `git add -A`.

## Pinned wheel versions

```bash
python3 -m venv ~/opt/gnina/cudaenv
~/opt/gnina/cudaenv/bin/pip install --no-cache-dir \
  nvidia-cuda-runtime-cu12==12.8.90 nvidia-cublas-cu12==12.8.4.1 \
  nvidia-cudnn-cu12==9.8.0.87 nvidia-cufft-cu12==11.3.3.83 \
  nvidia-cusolver-cu12==11.7.3.90 nvidia-cusparse-cu12==12.5.8.93 \
  nvidia-nvjitlink-cu12==12.8.93
```

These are the 12.8 line the binary was built against.

## PATH caveat

Ubuntu's default `~/.profile` prepends `~/.local/bin` only `if [ -d "$HOME/.local/bin" ]`, evaluated at
login. If that directory did not exist before this setup, `gnina` resolves only in shells started after
the next login — otherwise `export PATH="$HOME/.local/bin:$PATH"`.

## Confirming the install is intact

The first command must print nothing:

```bash
LD_LIBRARY_PATH=~/opt/gnina/cudalibs ldd ~/opt/gnina/gnina.bin | grep "not found"
gnina --version
gnina -r test/gnina/data/184l_rec.pdb -l test/gnina/data/184l_lig.sdf \
      --autobox_ligand test/gnina/data/184l_lig.sdf --seed 0 --cpu 8 -o /tmp/out.sdf
```

Measured timings on this machine are in [../performance.md](../performance.md), under the
locally-measured section.
