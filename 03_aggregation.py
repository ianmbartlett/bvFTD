import scanpy as sc
import pandas as pd

path = "/Users/ianbartlett/Desktop/school/uncw/mat495/00_pyth_data_loading/python_data/normal_Pick_2000genes.h5ad"
top_2000_genes = sc.read_h5ad(path)

# Create pseudobulk expression by donor and cell type
pb = sc.get.aggregate(
    top_2000_genes,
    by=["donor_id", "cell_type"],
    func=["sum", "mean"]
)

print(pb.layers.keys())

# Create donor-level metadata
donor_metadata = (
    top_2000_genes.obs[
        [
            "donor_id",
            "disease",
            "sex",
            "development_stage",
            "tissue_type",
            "self_reported_ethnicity"
        ]
    ]
    .drop_duplicates("donor_id")
)

pb.obs = pb.obs.reset_index(drop=True)

pb.obs = pb.obs.merge(
    donor_metadata,
    on="donor_id",
    how="left"
)

# Metadata columns to retain
metadata_cols = [
    "donor_id",
    "disease",
    "cell_type",
    "sex",
    "development_stage",
    "tissue_type",
    "self_reported_ethnicity",
    "n_obs_aggregated"
]

metadata = pb.obs[metadata_cols].reset_index(drop=True)

# Expression values from the summed pseudobulk layer
expression_df = pb.to_df(layer="sum").reset_index(drop=True)

# Combine metadata and gene expression
final_df = pd.concat(
    [metadata, expression_df],
    axis=1
)

print(final_df.head())
print(final_df.shape)

# Save CSV
final_df.to_csv(
    "/Users/ianbartlett/Desktop/school/uncw/mat495/r enviro/pseudobulk_expression_final.csv",
    index=False
)