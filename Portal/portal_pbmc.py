import portal
import scanpy as sc
import scipy.sparse as sp
import numpy as np
import rds2py
import scipy.sparse as ss
import rpy2.robjects as robjects
import h5py
import anndata


count_matrix_RNA = rds2py.read_rds('count_matrix_RNA.rds')
sp_rna_data = rds2py.as_sparse_matrix(count_matrix_RNA)

# Convert the sparse matrix to CSR format
sparse_matrix = sp.csr_matrix(sp_rna_data)

adata_RNA = sc.AnnData(X=sp_rna_data)
print(adata_RNA.obs)


h5 = h5py.File('gene_scores_ATACassays.h5', 'r')
se_data_ATAC = h5["assay001"]
np_data_ATA = se_data_ATAC[()]

adata_ATAC = anndata.AnnData(X=np_data_ATA)

# Print the shape of the AnnData object
print(adata_ATAC.obs)

# Create a Portal model
model = portal.model.Model()
model.preprocess(adata_RNA, adata_ATAC)  # Perform preprocess and PCA
model.train()  # Train the model
model.eval()  # Get integrated latent representation of cells


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

