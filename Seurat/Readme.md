## Seurat Workflow
_Author:_ Karina Díaz Barba

_Code File:_ `Seurat_10x.r`

This is the readme to understand how to perfome the code to analyse the PBMC 10X data and _create a Seurat object_. QC parameters are applied.

#### General comments
- R version (4.2.1)
- _Seurat version 4_ was used.
- Necessary libraries to run the analysis:
    - library(dplyr)
    - library(Seurat)
    - library(patchwork)
- Libraries to plot (optional):
    - library(ggplot2)
    - library(gridExtra)

#### Input files from the PBMC dataset:
- `barcodes.tsv.gz`
- `features.tsv.gz`
- `matrix.mtx.gz`

#### Output:
- A seurat object `pbmc_Seurat_Object.rds` with the information of the _nFeature_RNA_, _nCount_RNA_ and _percent.mt_ of every cell in the cell. `Input for in Azimuth`.
- Several violin plots graphically representing the cells removed after quality control.

- _Note:_ There are extra analyses for more data exploration but the code is commented because we do not use this in this project, as we have other tools to analyse the data. As an optinal step you can run those lines too and created several plots. 

