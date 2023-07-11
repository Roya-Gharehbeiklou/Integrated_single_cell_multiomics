# Project Name
Using Portal in Integrating single-cells multiomics data 
Portal codes are modified by @Roya Gharehbeikloo - r.gharehbeiklou@st.hanze.nl
Project Description

## Table of Contents
- [Installation](#installation)
- [Quick Start](#quick-start)
  - [Basic Usage](#basic-usage)
  - [Integrating Multiple Datasets](#integrating-multiple-datasets)
  - [Preparing the output file of portal](#preparing-the-output-file-of-portal)
- [Running on the Gearshift Cluster](#running-on-the-gearshift-cluster)

## Installation


### Setting up the Environment

1. Create a new conda environment named 'protal2' with Python 3.9 from the conda-forge channel:

    ```shell
    conda create -c conda-forge python=3.9 -n protal2
    ```

2. Activate the 'protal2' environment:

    ```shell
    conda activate portal2
    ```

3. Install 'llvmlite' and 'numba' using pip:

    ```shell
    pip install llvmlite
    pip install numba
    ```

4. Install R version 4.2.0 from the conda-forge channel:

    ```shell
    conda install -c conda-forge r-base=4.2.0
    ```

5. Install PyTorch and torchvision for CPU usage from the pytorch channel:

    ```shell
    conda install pytorch-cpu torchvision-cpu -c pytorch
    ```

6. Update conda to the latest version:

    ```shell
    conda update -n base -c defaults conda
    ```

7. Install 'portal-sc' package:

    ```shell
    conda install portal-sc
    ```

8. Set the R_HOME environment variable to the R installation in the conda environment:

    ```shell
    export R_HOME=~/.conda/envs/portal/bin/R
    ```

9. Create a symbolic link for the conda environment directory:

    ```shell
    ln -s /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Roya/.conda .conda
    ```

### Running the Portal Script

1. Execute the 'portal_pmbc.py' script:

    ```shell
    python portal_pmbc.py
    ```

### Copying Data

1. Copy the data to the specified destination:

    ```shell
    scp /data scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Roya/Portal_input/
    ```

### Retrieving Results

1. Retrieve the output.h5 file from the remote server:

    ```shell
    scp airlock+gearshift:/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Roya/Portal_input/output.h5
    ```

2. Retrieve the 're' directory from the remote server:

    ```shell
    scp airlock+gearshift:/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/students_hanze_2023/Users/Roya/Portal_input/re
    ```
## Quick Start
### Basic Usage
Starting with raw count matrices formatted as AnnData objects, Portal uses a standard pipline adopted by Seurat and Scanpy to preprocess data, followed by PCA for dimensionality reduction. After preprocessing, Portal can be trained via `model.train()`.
```python
import portal
import scanpy as sc

# read AnnData
adata_1 = sc.read_h5ad("adata_1.h5ad")
adata_2 = sc.read_h5ad("adata_2.h5ad")


model = portal.model.Model(training_steps=3000, npcs=20, n_latent=30, lambdacos=50.0)
#  At npcs=20, n_latent=30, lambdacos=50.0 and hvg_num=12000, we saw the best overlap
model.preprocess(adata_rna, adata_atac, hvg_num=12000) # perform preprocess and PCA

model.train() # train the model
model.eval() # get integrated latent representation of cells
```
The evaluating procedure `model.eval()` saves the integrated latent representation of cells in `model.latent`, which can be used for downstream integrative analysis.

#### Parameters in `portal.model.Model()`:
* `lambdacos`: Coefficient of the regularizer for preserving cosine similarity across domains. *Default*: `20.0`.
* `training_steps`: Number of steps for training. *Default*: `2000`. Use `training_steps=1000` for datasets with sample size < 20,000.
* `npcs`: Dimensionality of the embeddings in each domain (number of PCs). *Default*: `30`.
* `n_latent`: Dimensionality of the shared latent space. *Default*: `20`.
* `batch_size`: Batch size for training. *Default*: `500`.
* `seed`: Random seed. *Default*: `1234`.

The default setting of the parameter `lambdacos` works in general. We also enable tuning of this parameter to achieve a better performance, see [**Tuning `lambdacos` (optional)**](#tuning-lambdacos-optional). For the integration task where the cosine similarity is not a reliable cross-domain correspondance (such as cross-species integration), we recommend to use a lower value such as `lambdacos=10.0`.


### Integrating Multiple Datasets
Portal supports the integration of multiple datasets. If you have a list of AnnData objects, you can integrate them using the following commands:

```python
# loading RNA data
h5 = h5py.File('RNA_count.h5', "r")
h5_rna_barcodes = h5['10Xpbmc/barcodes']
h5_rna_features = h5['10Xpbmc/gene']
h5_rna_gene_names = h5['10Xpbmc/gene_names']

# loading ATAC data
data_atac = h5py.File('ATAC_scores.h5', "r")
h5_atac_barcodes = data_atac['matrix/barcodes']
h5_atac_features = data_atac['matrix/features']
h5_atac_cell_type = data_atac['matrix/celltypes']

# preprocess steps
rna_barcodes = np.array(h5_rna_barcodes).astype(str).copy()
rna_features = np.array(h5_rna_features).astype(str).copy()
rna_label = pd.read_csv('ancellTypes.csv', index_col = 0)

adata_rna = sc.AnnData(X=sp_rna_data.transpose())

adata_rna.obs.index = rna_barcodes
adata_rna.obs["cell_type"] = rna_label["x"].values.astype(str)
adata_rna.obs["data_type"] = "rna"
adata_rna.var.index = rna_features
print(adata_rna)



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

```

### Preparing the output file of portal

```python
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
```

## Running on the Gearshift Cluster
To run the Portal project on the Gearshift cluster, follow these steps:

- Set up a directory: Choose a location within your group's folder or designated storage area to store your project files. For example, create a directory in `/groups/your-group-name/`.

- Transfer files: Transfer all necessary project files, including the Portal code, data files, and any configuration files, to the cluster using tools like rsync or scp.

- Environment setup: Set up the required environment for running Portal on the cluster. Create a Conda environment with the necessary dependencies and activate it. Ensure that all required software and libraries are installed on the cluster.

- open a screen session on the cluster
    ```sh
    screen -S rserver
    ```
    Execute the 'portal_pmbc.py' script:

    ```shell
    python portal_pmbc.py
    ```

- Retrieve results: Once the job is completed, retrieve the output and results generated by Portal from the cluster. This may include integrated datasets, evaluation metrics, visualizations, or any other relevant files.



