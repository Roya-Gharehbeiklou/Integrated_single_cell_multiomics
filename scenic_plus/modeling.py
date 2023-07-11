"""
The module runs topic modeling. The modeling makes use of LDA (latent drichlet allocation) 
and has an idea of accounting cells that are very similar, but with slightly different features.
The modeling is computationaly intense.
The code is based on the tutorial for SCENIC+ 
(https://scenicplus.readthedocs.io/en/latest/pbmc_multiome_tutorial.html#Generate-pseudobulk-ATAC-seq-profiles,-call-peaks-and-generate-a-consensus-peak-set)
"""

import os
import pickle

from pycisTopic.pycisTopic.cistopic_class import *

#to disable ray's memory monitor
RAY_memory_monitor_refresh_ms=0

work_dir = '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Dilya/scenic_results/'
tmp_dir = 'scratch/'

#There are 2 types of LDA models: Serial (recommender for small-medium sized data sets)
# and Parallel with Mallet (recommended for large datsets).
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))

models=run_cgs_models(cistopic_obj,
                    #running more than 24 models was computationally challenging
                    n_topics=[10,16,20],
                    #basically, n_cpu should be equal to number of topic chosen
                    n_cpu=3,
                    n_iter=300,
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    save_path=None)

if not os.path.exists(os.path.join(work_dir, 'scATAC/models')):
    os.makedirs(os.path.join(work_dir, 'scATAC/models'))

pickle.dump(models,
            open(os.path.join(work_dir, 'scATAC/models/pbmc_models_LDA.pkl'), 'wb'))
