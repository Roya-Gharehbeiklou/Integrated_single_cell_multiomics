"""
The module contributes to scATAC-seq preprocessing. It creates cisTopic object.
Output: cisTopic object.
The code is based on the tutorial for SCENIC+ 
(https://scenicplus.readthedocs.io/en/latest/pbmc_multiome_tutorial.html#Generate-pseudobulk-ATAC-seq-profiles,-call-peaks-and-generate-a-consensus-peak-set)
"""

import os
import pickle
import scanpy as sc
import pandas as pd

from pycisTopic.pycisTopic.cistopic_class import *

work_dir = '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Dilya/scenic_results/'

#scRNA
adata = sc.read_h5ad(os.path.join(
    work_dir, '../azimuth_results/pbmc_Seurat_Object_QCfiltered_Azimuth.h5ad'))
#scRNA barcodes
scRNA_bc = adata.obs_names
cell_data = adata.obs
cell_data['sample_id'] = '10x_pbmc'
#set data type of the celltype column to str, otherwise the export_pseudobulk function will complain
cell_data['celltype'] = cell_data['predicted.celltype.l2'].astype(str)
del(adata)

#scATAC
fragments_dict = {'10x_pbmc': os.path.join(
    work_dir, '../../../data/output/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz')}
path_to_regions = {'10x_pbmc':os.path.join(
    work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed')}
path_to_blacklist = 'pycisTopic/blacklist/hg38-blacklist.v2.bed'

#scATAC barcodes
#ATAC barcodes from ArchR preprocessing and QC
atac_bc_archr = pd.read_csv('atac_barcodes.csv')
bc_passing_filters = {'10x_pbmc':[]}
bc_passing_filters['10x_pbmc'] = atac_bc_archr.iloc[:,0].tolist()
if not os.path.exists(os.path.join(work_dir, 'scATAC/quality_control')):
    os.makedirs(os.path.join(work_dir, 'scATAC/quality_control'))

pickle.dump(bc_passing_filters,
            open(os.path.join(work_dir, 'scATAC/quality_control/bc_passing_filters.pkl'), 'wb'))

print(f"{len(bc_passing_filters['10x_pbmc'])} barcodes passed QC stats")
print(f"{len(list(set(bc_passing_filters['10x_pbmc']) & set(scRNA_bc)))} cell barcodes pass both scATAC-seq and scRNA-seq based filtering")
#9763 barcodes passed QC stats
#9763 cell barcodes pass both scATAC-seq and scRNA-seq based filtering

# Creating cisTopic object
key = '10x_pbmc'
#please note, in tutorial QC metrics are also passed to the function
cistopic_obj = create_cistopic_object_from_fragments(
                            path_to_fragments = fragments_dict[key],
                            path_to_regions = path_to_regions[key],
                            path_to_blacklist = path_to_blacklist,
                            valid_bc = list(set(bc_passing_filters[key]) & set(scRNA_bc)),
                            n_cpu = 1,
                            project = key,
                            split_pattern = '-')
cistopic_obj.add_cell_data(cell_data, split_pattern='-')
print(cistopic_obj)
# CistopicObject from project 10x_pbmc with n_cells × n_regions = 9763 × 350257

#saving the object
pickle.dump(cistopic_obj, open(
    os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))
