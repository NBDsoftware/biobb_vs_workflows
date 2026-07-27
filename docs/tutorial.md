# Tutorial

A hands-on notebook that runs a complete **virtual screening** with
the workflows in this repository.

| Step | Command | What it does |
|------|---------|--------------|
| 1 | `cavity_analysis` | Extract a protein monomer of interest and perform a cavity analysis with Fpocket |
| 2 | `vs_autodock` | Perform a virtual screening on a selected pocket |

The tutorial can be run **either on Google Colab (no local install) or on a local Jupyter**.

## Run on Google Colab

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/NBDsoftware/biobb_vs_workflows/blob/master/notebooks/notebook_tutorial.ipynb)

Click the badge above and run the first setup cell. It builds the `biobb_vs_tutorial` conda
environment from `notebooks/colab_environment.yml` — nothing to install on your
machine. The conda solve is the slow part; run it once and wait.

## Run locally

Create the environment from `notebooks/local_environment.yml`, then launch Jupyter from it:

```bash
git clone https://github.com/NBDsoftware/biobb_vs_workflows.git
cd biobb_vs_workflows
conda env create -f notebooks/local_environment.yml
conda activate biobb_vs_tutorial
jupyter notebook notebooks/notebook_tutorial.ipynb
```

The notebook lives at
[`notebooks/notebook_tutorial.ipynb`](https://github.com/NBDsoftware/biobb_vs_workflows/blob/master/notebooks/notebook_tutorial.ipynb)
in the repository.
