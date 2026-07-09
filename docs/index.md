# BioBB VS Workflows

Ready-to-use command-line pipelines for structure-based virtual screening, built
on top of [BioExcel Building Blocks (BioBB)](https://mmb.irbbarcelona.org/biobb/).
Developed for the [European BioExcel](http://bioexcel.eu/) project, funded by the
European Commission (EU Horizon Europe
[101093290](https://cordis.europa.eu/project/id/101093290)).

This repository covers cavity detection and high-throughput docking. Related
pipelines live in separate repositories: MD and preparation workflows in
[biobb_md_workflows](https://github.com/NBDsoftware/biobb_md_workflows).

```{image} _static/bioexcel_logo.png
:alt: BioExcel
:width: 250px
:target: https://bioexcel.eu/
:align: center
:class: only-light
```
```{image} _static/bioexcel_logo_white.png
:alt: BioExcel
:width: 250px
:target: https://bioexcel.eu/
:align: center
:class: only-dark

## Workflows

| Command | Purpose |
|---------|---------|
| [`cavity_analysis`](workflows/cavity_analysis.md) | Cluster an MD trajectory and detect + filter cavities with Fpocket on the representative structures. |
| [`vs_autodock`](workflows/vs_autodock.md) | High-throughput virtual screening of a ligand library against a pocket with AutoDock Vina. |

```{toctree}
:maxdepth: 2
:hidden:

installation
workflows/cavity_analysis
workflows/vs_autodock
```

## Licensing

Offered under a dual-license model: free for academic and non-commercial use under
**CC BY-NC-SA 4.0**; a separate commercial license is required for for-profit use
(contact `it@nostrumbiodiscovery.com`). See the `LICENSE` file in the repository.
