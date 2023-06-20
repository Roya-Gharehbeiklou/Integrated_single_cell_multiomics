import rds2py
import scipy.sparse as ss
import h5sparse
import numpy as np

# RNA counts to h5

rObj_RNA = rds2py.read_rds('count_matrix_RNA.rds')
sp_mat_RNA = rds2py.as_sparse_matrix(rObj_RNA)

with h5sparse.File("count_matrix_RNA.h5", "w") as h5f:
     h5f.create_dataset('sparse/matrix', data=sp_mat_RNA)

# ATAC gene score matrix to h5
import h5py

h5 = h5py.File('gene_scores_ATACassays.h5', 'r')
