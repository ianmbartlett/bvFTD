import scanpy as sc
import pandas as pd

# Initial Data Loading and Counts of the AnnData object
path = "/Users/ianbartlett/Desktop/school/uncw/mat495/00_pyth_data_loading/python_data/RNA_seq_fixed.h5ad"
adata = sc.read_h5ad(path)

print(adata)
print(adata.obs_names[:5])
print(adata.obs_names.is_unique)

print(adata.shape)

print(adata.obs["disease"].value_counts())

# Filtering out the diseases of no interest to leave bvFTD patients
adata_subset = adata[adata.obs["disease"].isin(
            ["normal",
             "Pick disease"]
)].to_memory()

# Confirming the filtering worked
print(adata_subset.obs["disease"].value_counts())

adata_ml = adata_subset.copy()

adata_ml.raw = None
adata_ml.obsm.clear()

print(adata_ml)
print("Raw exists:", adata_ml.raw is not None)
print("Obsm:", list(adata_ml.obsm.keys()))

# Filtering the top 2000 genes for a seurat-type object
sc.pp.highly_variable_genes(
    adata_ml,
    n_top_genes=2000,
    flavor="seurat_v3"
)

#Writing the smaller h5ad dataset

adata_small = adata_ml[:, adata_ml.var["highly_variable"]].copy()

adata_small.write_h5ad(
    "/Users/ianbartlett/Desktop/school/uncw/mat495/00_pyth_data_loading/python_data/normal_Pick_2000genes.h5ad"
)

