import scanpy as sc
import numpy as np
from scipy import sparse

path = "/Users/ianbartlett/Desktop/school/uncw/mat495/00_pyth_data_loading/python_data/normal_Pick_2000genes.h5ad"

top_2000_genes = sc.read_h5ad(path)

print(top_2000_genes)

print("Nonzero entries:", top_2000_genes.X.nnz)

print("Minimum nonzero:", top_2000_genes.X.data.min())
print("Maximum nonzero:", top_2000_genes.X.data.max())
print("First 20 nonzero values:", top_2000_genes.X.data[:20])

# Data already normalized and log-transformed --> Moving onto Scaling

# Data analysis for pseudobulk dataset creation.

print(top_2000_genes.obs.groupby(["donor_id", "disease"]).size())

celltype_counts = (
    top_2000_genes.obs
    .groupby(["donor_id", "cell_type"])
    .size()
    .unstack(fill_value=0)
)


