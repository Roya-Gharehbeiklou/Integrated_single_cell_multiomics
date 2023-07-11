"""
Analysing results from topic modeling and selecting the best number of topics.
Creating UMAPs for the resulting model.
The code is based on the tutorial for SCENIC+ 
(https://scenicplus.readthedocs.io/en/latest/pbmc_multiome_tutorial.html#Generate-pseudobulk-ATAC-seq-profiles,-call-peaks-and-generate-a-consensus-peak-set)
"""

import os
import pickle
from pycisTopic.pycisTopic.lda_models import *

work_dir = '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Dilya/scenic_results/'

models = pickle.load(open(os.path.join(work_dir, 'scATAC/models/pbmc_models_LDA.pkl'), 'rb'))
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))

model = evaluate_models(models,
                       select_model=16,
                       return_model=True,
                       metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                       plot_metrics=False,
			save='models')

cistopic_obj.add_LDA_model(model)
pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))

#Visualization
run_umap(cistopic_obj, target  = 'cell', scale=True)
plot_metadata(cistopic_obj, reduction_name = 'UMAP', variables = ['celltype'], save='cistopic_umap')

#use of cell-topic probabilities to visualize cell type specificity
plot_topic(cistopic_obj, reduction_name = 'UMAP', num_columns = 4, save='topic_umaps')