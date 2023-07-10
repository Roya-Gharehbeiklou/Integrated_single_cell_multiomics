## FigR workflow 

This folder contain the files that are used during the gene regulatory network construction using `FigR`.

_Libraries used in each file are listed at the top of the file. Loading libraries in a different order could produce dependency issues_

Please take a close look at the files listed below which should be ran in the same order as listed. File paths are embedded in the scripts:
1. `Preprocessing.R` - Loads `Seurat` and `ArchR` preprocessed data and generates correct input data for `FigR`. Also generates correct input format for data integration performed by `Portal`.
  - Input:
    - `pbmc_Seurat_Azimuth_for_figR.rds`
    - `pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5`
    - `Save-ArchR-Project_subSet_QC_Frip.rds`
  - Output, `Portal` input:
    - `gene_scores_ATACassays.h5`
    - `RNA_count.h5`
  - Output `FigR`:
    - `ATAC_se.rds`
    - `RNA_mat.rds`
2. `Peak-correlations.R` - Calculates gene-peak correlations between all genes for all cells.
  - Input:
    - `ATAC_se.rds`
    - `RNA_mat.rds`
  - Output:
    - `ciscor.csv`
  
