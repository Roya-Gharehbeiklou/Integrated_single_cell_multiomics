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
np.random.seed(1234)

sc.settings.verbosity = 3
sc.logging.print_header()


h5 = h5py.File('RNA_count.h5', "r")
h5_rna_data = h5['10Xpbmc/data']
h5_rna_barcodes = h5['10Xpbmc/barcodes']
h5_rna_features = h5['10Xpbmc/gene']
h5_rna_gene_names = h5['10Xpbmc/gene_names']
import pandas as pd

rna_data = scipy.sparse.csr_matrix(np.array(h5_rna_data).transpose()).copy()
rna_barcodes = np.array(h5_rna_barcodes).astype(str).copy()
rna_features = np.array(h5_rna_features).astype(str).copy()
rna_label = pd.read_csv('ancellTypes.csv', index_col = 0)

adata_rna = anndata.AnnData(rna_data)
print("1--------", adata_rna)

adata_rna.obs.index = rna_barcodes
adata_rna.obs["cell_type"] = rna_label["x"].values.astype(str)
adata_rna.obs["data_type"] = "rna"
adata_rna.var.index = rna_features

print(adata_rna)
adata_atac = anndata.read_hdf('gene_scores_ATACassays.h5', key="assay001")
adata_atac.obs.index = rna_barcodes
adata_atac.obs["cell_type"] = rna_label["x"].values.astype(str)
adata_atac.obs["data_type"] = "atac"

adata_atac.var.index = rna_features

print(adata_atac)

meta_rna = adata_rna.obs
meta_atac = adata_atac.obs

meta = pd.concat([meta_rna, meta_atac], axis=0)# Print the shape of the AnnData object


# Create a folder for saving results
result_path = "./result"
if not os.path.exists(result_path):
    os.makedirs(result_path)
    
    
model = portal.model.Model(training_steps=3000, lambdacos=10.0)
model.preprocess(adata_rna, adata_atac) # perform preprocess and PCA
model.train() # train the model
model.eval() # get integrated latent representation of cells

portal.utils.plot_UMAP(model.latent, meta, colors=["data_type", "cell_type"], save=False, result_path=result_path)


# Perform downstream analyses and visualizations using the integrated data
integrated_data = model.integrated

# Save the integrated and analyzed data
integrated_data.write("integrated_and_analyzed_data.h5ad")

# Set up PDF file for saving figures
pdf_path = "figures.pdf"
pdf = PdfPages(pdf_path)

# UMAP visualization
sc.pp.neighbors(integrated_data)
sc.tl.umap(integrated_data)
sc.pl.umap(integrated_data, color=['RNA', 'ATAC'])
pdf.savefig()

# Cluster analysis
sc.pp.neighbors(integrated_data)
sc.tl.leiden(integrated_data)
sc.pl.umap(integrated_data, color=['leiden'])
pdf.savefig()

# Differential gene expression analysis
sc.tl.rank_genes_groups(integrated_data, groupby='leiden')
sc.pl.rank_genes_groups(integrated_data, n_genes=10, sharey=False)
pdf.savefig()

# Heatmap of gene expression
genes_of_interest = ['Gene1', 'Gene2', 'Gene3']
sc.pl.heatmap(integrated_data, var_names=genes_of_interest, groupby='leiden', cmap='viridis')
pdf.savefig()

# Save the PDF file
pdf.close()

