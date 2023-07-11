## FigR workflow 

This folder contain the files that are used during the gene regulatory network construction using `FigR`.

_Libraries used in each file are listed at the top of the file. Loading libraries in a different order could produce dependency issues_

Please take a close look at the files listed below which should be ran in the same order as listed. File paths are embedded in the scripts:
1. `Preprocessing.R` - Loads `Seurat` and `ArchR` preprocessed data and generates correct input data for `FigR`. Also generates correct input format for data integration performed by `Portal`.
  - Input:
    - `pbmc_Seurat_Azimuth_for_figR.rds`: `Seurat` preprocessed object with Azimuth cell annotations.
    - `pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5`: Filtered barcode matrix containing RNA counts (non-normalized).
    - `Save-ArchR-Project_subSet_QC_Frip.rds` - `ArchR` preprocessed object.
  - Output, `Portal` input:
    - `gene_scores_ATACassays.h5` - Gene scores matrix from `AtchR`.
    - `RNA_count.h5` - RNA Count matrix. 
  - Output `FigR`:
    - `ATAC_se.rds` - Summarized Experiment with ATAC count matrix.
    - `RNA_mat.rds` - RNA normalized count matrix. Equa barcodes among ATAC and RNA object.
2. `Peak-correlations.R` - Calculates gene-peak correlations between all genes for all cells.
  - Input:
    - `ATAC_se.rds`: Summarized Experiment with ATAC count matrix.
    - `RNA_mat.rds`: RNA normalized count matrix. Equa barcodes among ATAC and RNA object.
  - Output:
    - `ciscor.csv`: Spearman correlation across all cells between their peak accessibility counts (mean-centered) and the normalized RNA expression.
3. `Peak-correlations.sh` - Helper script for batch execution of  `Peak-correlations.R`
4. `Create_cisTopic.R` - Build cisTopic model by running several modeles with different number of topics.
  - Input:
      - `ATAC_se.rds`: Summarized Experiment with ATAC count matrix.
  - Output:
      - `best_cisTopicObject.rds`: cisTopicObject containing the best performing model accompanied by metadata.
5. `Create_cisTopic.sh` - Helper script for batch execution of `Create_cisTopic.R`
6. `DORCs-Smoothing.R` - Smooth input data gene regulatory network construction using cisTopicObject generated in `Create_cisTopic.R`
  - Input:
    - `best_cisTopicObject.rds`: cisTopicObject containing the best performing model accompanied by metadata.
    - `ciscor.csv`: Spearman correlation across all cells between their peak accessibility counts (mean-centered) and the normalized RNA expression.
    - `ATAC_se.rds`: Summarized Experiment with ATAC count matrix.
    - `RNA_mat.rds`: RNA normalized count matrix. Equa barcodes among ATAC and RNA object.
  - Output:
    - `dorcMAt_smoothed.rds`: Smoothed DORC matrix object unsing KNN.
    - `RNAMat_smoothed.rds`: Smoothed RNA matrix object using KNN.
7. `FigR_GRN_function.R` - Altered function of `runFigRGRN`. Changes lib paths.
8. `Run_GRN.R` - Created GRN using `FigR_GRN_function.R`.
  - Input:
    - `ATAC_se.rds`: Summarized Experiment with ATAC count matrix.
    - `RNA_mat.rds`: RNA normalized count matrix. Equa barcodes among ATAC and RNA object.
    - `ciscor.csv`: Spearman correlation across all cells between their peak accessibility counts (mean-centered) and the normalized RNA expression.
    - `dorcMAt_smoothed.rds`: Smoothed DORC matrix object unsing KNN.
    - `RNAMat_smoothed.rds`: Smoothed RNA matrix object using KNN.
  - Output:
    - `figRGRN.rds`: GRN of the preprocessed data.
9. `Visualize_results.R` - Generates visualization for the generated GRN.
  - Input:
    - `figRGRN.rds`: GRN of the preprocessed data.
  - Output:
    - Several plots of the GRN.

    

