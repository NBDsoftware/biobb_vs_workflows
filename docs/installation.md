# Installation

Requirements: `git`, `conda`.

```bash
git clone https://github.com/NBDsoftware/biobb_vs_workflows.git
cd biobb_vs_workflows
conda env create -f environment.yml
conda activate biobb_vs
```

This creates the `biobb_vs` conda environment and installs the package, exposing
the two workflow commands (`cavity_analysis`, `vs_autodock`).
