## SCENIC+  

The code is based on the tutorial for SCENIC+ 
(https://scenicplus.readthedocs.io/en/latest/pbmc_multiome_tutorial.html#Generate-pseudobulk-ATAC-seq-profiles,-call-peaks-and-generate-a-consensus-peak-set)  
SCENIC+ uses scATAC-seq data, scRNA-seq data and cell type annotation in order to reconstruct the GRN.  
The work with the novelle tool for GRN reconstruction started with the scATAC preprocessing. However, it was already performed with the ArchR tool, it was impossible to implement the output from ArchR into SCENIC+'s workflow.  
- scATAC preprocessing (pycisTopic module):  
    - creating pseudobulk profiles per cell type (pseudobulk.py)
    - calling peaks per pseudobulk profile and combining them into consensus peak set (peak.py)
- cisTopic creation (pycisTopic module):  
    - creation of cisTopic object (cistopic.py)
    - topic modeling (modeling.py)
    - analysing models (model_analysis.py)
