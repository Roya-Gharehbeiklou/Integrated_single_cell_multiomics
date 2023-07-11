## ArchR Workflow
_Author:_ Karina Díaz Barba

_Code File:_ `ArchR.r`

This is the readme to understand how to perfomed the code to analyse the PBMC 10X data and create a `ArchR object` with the data 
filtered based on specific QC parameters and the FRIP values for every cell.

#### General comments
- R version (4.2.1)
- _ArchR version v1.0.2_ was used.
- Necessary libraries to run the analysis:
    - library(ArchR)
    - library(IRanges)
    - library(GenomicRanges)
- Reference genome:
    - library(`BSgenome.Hsapiens.UCSC.hg38`)

#### Notes: 
- The path of the location of the libraries needed to run ArchR should be defined  as `.libPaths`. As specified in the code. 
- The First time running ArchR might be necessary to install extra libraries. If that happens run the next command:
`ArchR::installExtraPackages()`
- `Iranges` library needs to be loaded separately. As specified in the code.
- The function of _seqnames_ from _GenomicRanges_ should be loaded as: `seqnames <- GenomicRanges::seqnames`. As specified in the code.
- The number of threads should be set as one `addArchRThreads(threads = 1)` . Theres a bug in the program. If this is not set like this, the  doublets inference will not calculated and there will be a problem creating the arrow files. Read: `https://github.com/GreenleafLab/ArchR/issues/218`  for more information. 


#### Input files from the PBMC dataset:
- `pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz`

#### Output:
- `Save-ArchR-Project.rds` ArchR object with information of different Matrices.
- `Save-ArchR-Project_subSet_QC_Frip.rds`: Subset of the ArchR object. Filtered by FRIP>=0.5.
- `ATAC_se.rds`: Summarized Experiment with ATAC count matrix.
-  Plots of the QC parameters performed. 

