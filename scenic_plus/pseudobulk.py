"""
The module contributes to scATAC-seq preprocessing. It generates pseudobulk profiles 
per cell type. Two files are generated:
- pseudobulk in bed format
- pseudobulk in BigWig format
The code is based on the tutorial for SCENIC+ 
(https://scenicplus.readthedocs.io/en/latest/pbmc_multiome_tutorial.html#Generate-pseudobulk-ATAC-seq-profiles,-call-peaks-and-generate-a-consensus-peak-set)
"""

import os
import warnings
import pickle
import scanpy as sc
import pyranges as pr
import pandas as pd

from pycisTopic.pycisTopic.pseudobulk_peak_calling import export_pseudobulk

warnings.simplefilter(action='ignore')

work_dir = '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Dilya/scenic_results/'
tmp_dir = 'scratch/'

#load the cell type annotation from scRNA
adata = sc.read_h5ad(os.path.join(
    work_dir, '../azimuth_results/pbmc_Seurat_Object_QCfiltered_Azimuth.h5ad'))
cell_data = adata.obs
cell_data['sample_id'] = '10x_pbmc'
#set data type of the celltype column to str, otherwise the export_pseudobulk function will complain
cell_data['celltype'] = cell_data['predicted.celltype.l2'].astype(str)
del(adata)

#generate pseudobulk profiles
#get chromosome sizes (hg38)
target_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
chromsizes = pd.read_csv(target_url, sep='\t', header=None)
chromsizes.columns = ['Chromosome', 'End']
chromsizes['Start'] = [0] * chromsizes.shape[0]
chromsizes = chromsizes.loc[:,['Chromosome', 'Start', 'End']]
# Exceptionally in this case, to agree with CellRangerARC annotations
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace(
    'v', '.') for x in range(len(chromsizes['Chromosome']))]
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].split('_')[1] if len(
    chromsizes['Chromosome'][x].split('_')) > 1 else chromsizes['Chromosome'][x] for x in range(
        len(chromsizes['Chromosome']))]
chromsizes=pr.PyRanges(chromsizes)

#location of scATAC-seq raw data
fragments_dict = {'10x_pbmc': os.path.join(
    work_dir, '../../../data/output/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz')}

bw_paths, bed_paths = export_pseudobulk(input_data = cell_data,
                 # variable by which to generate pseubulk profiles
                 variable = 'celltype',
                 sample_id_col = 'sample_id',
                 chromsizes = chromsizes,
                 # specify where pseudobulk_bed_files should be stored
                 bed_path = os.path.join(work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bed_files/'),
                 # specify where pseudobulk_bw_files should be stored
                 bigwig_path = os.path.join(work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bw_files/'),
                 # location of fragment files
                 path_to_fragments = fragments_dict,
                 n_cpu = 1,
                 normalize_bigwig = True,
                 remove_duplicates = True,
                 _temp_dir = os.path.join(tmp_dir, 'ray_spill'),
                 split_pattern = '-')

#Saving bed and bigwig files as dictionaries
pickle.dump(bed_paths, open(os.path.join(
    work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl'), 'wb'))
pickle.dump(bw_paths, open(os.path.join(
    work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bed_files/bw_paths.pkl'), 'wb'))
