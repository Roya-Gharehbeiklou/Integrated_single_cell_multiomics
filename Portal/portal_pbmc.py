import os
import portal
import scanpy as sc
import numpy as np
import scipy.sparse as sp
import rds2py
import anndata2ri
import h5py
import anndata
import scipy.io
import pandas as pd

np.random.seed(1234)

sc.settings.verbosity = 3
sc.logging.print_header()


count_matrix_RNA = rds2py.read_rds('count_matrix_RNA.rds')
sp_rna_data = rds2py.as_sparse_matrix(count_matrix_RNA)

h5 = h5py.File('RNA_count.h5', "r")
# h5_rna_data = h5['10Xpbmc/data']
h5_rna_barcodes = h5['10Xpbmc/barcodes']
h5_rna_features = h5['10Xpbmc/gene']
h5_rna_gene_names = h5['10Xpbmc/gene_names']

# rna_data = scipy.sparse.csr_matrix(np.array(h5_rna_data).transpose()).copy()
rna_barcodes = np.array(h5_rna_barcodes).astype(str).copy()
rna_features = np.array(h5_rna_features).astype(str).copy()
rna_label = pd.read_csv('ancellTypes.csv', index_col = 0)

adata_rna = sc.AnnData(X=sp_rna_data.transpose())

adata_rna.obs.index = rna_barcodes
adata_rna.obs["cell_type"] = rna_label["x"].values.astype(str)
adata_rna.obs["data_type"] = "rna"
adata_rna.var.index = rna_features
print(adata_rna)

adata_atac = anndata.read_hdf('gene_scores_ATACassays.h5', key="assay001")
data_atac = h5py.File('ATAC_scores.h5', "r")
h5_atac_barcodes = data_atac['matrix/barcodes']
h5_atac_features = data_atac['matrix/features']
h5_atac_cell_type = data_atac['matrix/celltypes']

atac_barcodes = np.array(h5_atac_barcodes).astype(str).copy()
atac_features = np.array(h5_atac_features).astype(str).copy()
atac_cell_type = np.array(h5_atac_cell_type).astype(str).copy()

adata_atac.obs.index = atac_barcodes
adata_atac.obs["cell_type"] = atac_cell_type
adata_atac.obs["data_type"] = "atac"
adata_atac.var.index = atac_features
print(adata_atac)

meta_rna = adata_rna.obs
meta_atac = adata_atac.obs

meta = pd.concat([meta_rna, meta_atac], axis=0)# Print the shape of the AnnData object


# Create a folder for saving results
result_path = "./result"
if not os.path.exists(result_path):
    os.makedirs(result_path)
    
    
model = portal.model.Model(training_steps=3000, npcs=20, n_latent=30, lambdacos=50.0)
#  At npcs=20, n_latent=30, lambdacos=50.0 and hvg_num=12000, we saw the best overlap
model.preprocess(adata_rna, adata_atac, hvg_num=12000) # perform preprocess and PCA
model.train() # train the model
model.eval() # get integrated latent representation of cells

portal.utils.plot_UMAP(model.latent, meta, colors=["data_type", "cell_type"], save=True, result_path=result_path)

import umap
import matplotlib.pyplot as plt


reducer = umap.UMAP(n_neighbors=30,
                    n_components=2,
                    metric="correlation",
                    n_epochs=None,
                    learning_rate=1.0,
                    min_dist=0.3,
                    spread=1.0,
                    set_op_mix_ratio=1.0,
                    local_connectivity=1,
                    repulsion_strength=1,
                    negative_sample_rate=5,
                    a=None,
                    b=None,
                    random_state=1234,
                    metric_kwds=None,
                    angular_rp_forest=False,
                    verbose=True)

embedding = reducer.fit_transform(model.latent)

n_cells = embedding.shape[0]
size = 120000 / n_cells


import os
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn import preprocessing
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib as mpl
import umap
import anndata

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

rgb_10 = [i for i in get_cmap('Set3').colors]
rgb_20 = [i for i in get_cmap('tab20').colors]
rgb_20b = [i for i in get_cmap('tab20b').colors]

rgb2hex_10 = [mpl.colors.rgb2hex(color) for color in rgb_10]
rgb2hex_20 = [mpl.colors.rgb2hex(color) for color in rgb_20]
rgb2hex_20b = [mpl.colors.rgb2hex(color) for color in rgb_20b]
rgb2hex_20b_new = [rgb2hex_20b[i] for i in [0, 3, 4, 7, 8, 11, 12, 15, 16, 19]]
rgb2hex = rgb2hex_20 + rgb2hex_20b_new


meta.loc[meta.cell_type.values.astype(str) == "unknown", "cell_type"] = "Unknown"
meta.loc[meta.data_type.values.astype(str) == "rna", "data_type"] = "rna-seq"
meta.loc[meta.data_type.values.astype(str) == "atac", "data_type"] = "atac-seq"


n_cells = embedding.shape[0]
if n_cells >= 15000:
    size = 120000 / n_cells
else:
    size = 8

le = preprocessing.LabelEncoder()
le.fit(sorted(set(meta["data_type"])))
label = le.fit_transform(meta["data_type"].values)
colours = ListedColormap(["tab:blue", "tab:orange"])

le2 = preprocessing.LabelEncoder()
le2.fit(sorted(set(meta["cell_type"])))
label2 = le.fit_transform(meta["cell_type"].values)
colours2 = ListedColormap(rgb2hex_20b_new)


f = plt.figure(figsize=(20,10))

ax1 = f.add_subplot(1,2,1)
scatter1 = ax1.scatter(embedding[:, 0][::-1], embedding[:, 1][::-1], s=size, c=label[::-1], cmap=colours, label=meta["data_type"].values[::-1])
ax1.set_title("Method", fontsize=40)
ax1.tick_params(axis='both',bottom=False, top=False, left=False, right=False, labelleft=False, labelbottom=False, grid_alpha=0)
# celltype
ax12 = f.add_subplot(1,2,2)
scatter2 = ax12.scatter(embedding[:, 0], embedding[:, 1], s=size, c=label2, cmap=colours2, label=meta["cell_type"].values)
ax12.set_title("Cell type", fontsize=40)
ax12.tick_params(axis='both',bottom=False, top=False, left=False, right=False, labelleft=False, labelbottom=False, grid_alpha=0)


l1 = f.legend(handles=scatter1.legend_elements()[0], labels=sorted(set(meta["data_type"])), loc="upper left", bbox_to_anchor=(1.0, 0.9), 
              markerscale=3., title_fontsize=30, fontsize=20, frameon=False, ncol=1, title="Method")
l2 = f.legend(handles=scatter2.legend_elements(num=len(sorted(set(meta["cell_type"]))))[0], labels=sorted(set(meta["cell_type"])), loc="upper left", bbox_to_anchor=(1.0, 0.7), 
              markerscale=3., title_fontsize=30, fontsize=20, frameon=False, ncol=1, title="Cell type")
l1._legend_box.align = "left"
l2._legend_box.align = "left"

f.subplots_adjust(hspace=.1, wspace=.1)

plt.savefig('result.png', bbox_inches='tight')

file = h5py.File('output.h5', 'w')

# Create a group for the DataFrame
df_group = file.create_group('meta')

# Convert the DataFrame to a NumPy array and save it within the group
df_group.create_dataset('values', data=meta.values)
df_group.create_dataset('index', data=meta.index.values)

# Create a group for the NumPy array
arr_group = file.create_group('model')

# Save the NumPy array within the group
arr_group.create_dataset('latent', data=model.latent)

# Close the HDF5 file
file.close()