# Bioexcel Building Block workflows — virtual screening

BioBB workflows are ready-to-use pipelines built using BioExcel Building Blocks (BioBB) to perform common tasks in biomolecular simulations and modeling. This software has been developed for the [European BioExcel](http://bioexcel.eu/), funded by the European Commission (EU Horizon Europe [101093290](https://cordis.europa.eu/project/id/101093290)).

This repo covers the virtual-screening workflows. MD/preparation workflows (`md_gromacs`, `ligand_parameterization`, `protein_preparation`, `traj_postprocessing`) live in [biobb_md_workflows](https://github.com/NBDsoftware/biobb_md_workflows).

## Workflows

**cavity_analysis**: clusters an MD trajectory and runs a cavity analysis using Fpocket on the representative structures.

**vs_autodock**: high-throughput virtual screening of a selected pocket using a ligand library (SMILES/SDF) and AutoDock Vina.

## Installation

Requirements: `git`, `conda`

```bash
git clone https://github.com/NBDsoftware/biobb_vs_workflows.git
cd biobb_vs_workflows
conda env create -f environment.yml
conda activate biobb_vs
```

## Usage

Once installed, each workflow is available as a CLI command:

| Command | Workflow |
|---------|----------|
| `cavity_analysis` | cavity_analysis |
| `vs_autodock` | vs_autodock |

```bash
cavity_analysis --help
```

## Licensing

This project is offered under a dual-license model, intended to make the software freely available for academic and non-commercial use while preventing its use for profit.

### 1. Academic and Non-Commercial Use

For academic, research, and other non-commercial purposes, this software is licensed under the **Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)**.

Under this license, you are free to:
*   **Share** — copy and redistribute the material in any medium or format.
*   **Adapt** — remix, transform, and build upon the material.

As long as you follow the license terms:
*   **Attribution** — You must give appropriate credit.
*   **Non-Commercial** — You may not use the material for commercial purposes.
*   **Share-Alike** — If you remix, transform, or build upon the material, you must distribute your contributions under the same license as the original.

A full copy of the license is available in the [LICENSE](LICENSE) file in this repository.

### 2. Commercial Use

**Use of this software for commercial purposes is not permitted under the CC BY-NC-SA 4.0 license.**

If you wish to use this software in a commercial product, for-profit service, or any other commercial context, you must obtain a separate commercial license.

Please contact **it@nostrumbiodiscovery.com** to inquire about purchasing a commercial license.

![](https://bioexcel.eu/wp-content/uploads/2019/04/Bioexcell_logo_1080px_transp.png "Bioexcel")
