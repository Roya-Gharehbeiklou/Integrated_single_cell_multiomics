"""
The module contributes to scATAC-seq preprocessing. It calls peaks per 
pseudobulk profile (output of pseudobulk.py), which are then combined 
into one consensus peak set.
Output: consensus peak set.
The code is based on the tutorial for SCENIC+ 
(https://scenicplus.readthedocs.io/en/latest/pbmc_multiome_tutorial.html#Generate-pseudobulk-ATAC-seq-profiles,-call-peaks-and-generate-a-consensus-peak-set)
"""

import os
import warnings
import pickle
import pandas as pd
import pyranges as pr

# to run pseudobulk_peak_calling the following modules should be installed: pyBigWig,
#pyranges, ray, lda (in conda env load module GCC)
from pycisTopic.pycisTopic.pseudobulk_peak_calling import peak_calling
from pycisTopic.pycisTopic.iterative_peak_calling import *

warnings.simplefilter(action='ignore')

#Project directory - used as default root directory
work_dir = '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Dilya/scenic_results/'
#Output directory
outDir = work_dir + 'output/'
if not os.path.exists(outDir):
    os.makedirs(outDir)
#temporary files directory
tmp_dir = 'scratch/'

#Call peaks per pseudobulk profile
bed_path = pickle.load(open(os.path.join(
    work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl'), 'rb'))
bw_path =  pickle.load(open(os.path.join(
    work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bed_files/bw_paths.pkl'), 'rb'))

#MACS2 - calculates p-value for each peak
macs_path = 'macs2'

# Run peak calling
narrow_peaks_dict = peak_calling(macs_path,
                                 bed_path,
                                 os.path.join(work_dir,'scATAC/consensus_peak_calling/MACS/'),
                                 #genome size which can be sequenced
                                 genome_size='hs',
                                 #changed cpu from 8 to 1, because otherwise the function won't run
                                 n_cpu=1,
                                 input_format='BEDPE',
                                 shift=73,
                                 ext_size=146,
                                 keep_dup = 'all',
                                 q_value = 0.05,
                                 _temp_dir = os.path.join(tmp_dir, 'ray_spill'))

#saving the dictionary
pickle.dump(narrow_peaks_dict, open(os.path.join(
    work_dir, 'scATAC/consensus_peak_calling/MACS/narrow_peaks_dict.pkl'), 'wb'))

#get chromosome sizes
target_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
chromsizes = pd.read_csv(target_url, sep='\t', header=None)
chromsizes.columns = ['Chromosome', 'End']
chromsizes['Start'] = [0]*chromsizes.shape[0]
chromsizes = chromsizes.loc[:,['Chromosome', 'Start', 'End']]
# Exceptionally in this case, to agree with CellRangerARC annotations
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace(
    'v', '.') for x in range(len(chromsizes['Chromosome']))]
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].split('_')[1] if len(
    chromsizes['Chromosome'][x].split('_')) > 1 else chromsizes['Chromosome'][x] for x in range(
        len(chromsizes['Chromosome']))]
chromsizes = pr.PyRanges(chromsizes)

#merging peaks into consensus peak set
peak_half_width = 250 #Number of base pairs that each summit will be extended in each direction
path_to_blacklist = os.path.join(work_dir, 'pycisTopic/blacklist/hg38-blacklist.v2.bed')
# Get consensus peaks
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict, peak_half_width, chromsizes=chromsizes, path_to_blacklist=path_to_blacklist)

#saving peak set as bed file
consensus_peaks.to_bed(
    path = os.path.join(work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed'),
    keep=True,
    compression='infer',
    chain=False)
