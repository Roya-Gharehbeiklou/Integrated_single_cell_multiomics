import portal
import scanpy as sc
import scipy.sparse as sp
import rds2py
import anndata2ri
import h5py
import anndata as ad
import rpy2.robjects as ro

ro.r['library']("DropSeq.util")
ro.r['library']("Seurat")
ro.r['library']("SingleCellExperiment")
ro.r.assign("obj.path", "pbmc.RData")
ro.r("load(obj.path)")
ro.r("cnt <- pbmc@assays$RNA@counts")
ro.r("celltype <- as.vector(pbmc$celltype)")
ro.r("rm(pbmc)")
ro.r("obj_sc <- CreateSeuratObject(counts=cnt)")
ro.r("names(celltype) <- colnames(x = obj_sc)")
ro.r("obj_sc <- AddMetaData(object = obj_sc, metadata = celltype, col.name = 'celltype')")

ro.r("sce <- as.SingleCellExperiment(obj_sc)")
anndata2ri.activate()
adata_sc = ro.r('as(sce, "SingleCellExperiment")')

adata_RNA = adata_sc
print(adata_RNA)


adata_ATAC = ad.read_hdf('gene_scores_ATACassays.h5', key="assay001")
print(adata_ATAC)


# Print the shape of the AnnData object

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

