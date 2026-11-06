#Pseudobulk Assessment
library(dplyr)
donor_cell.data <- read.csv("pseudobulk_expression_final.csv")
head(donor_cell.data)

# Removing synonymous columns


# Remove multiple columns
donor_cell.data_v2 <- donor_cell.data %>% select(-c(tissue_type, self_reported_ethnicity))
donor_cell.data_factor <- donor_cell.data_v2
donor_cell.data_factor$disease <- 
  factor(donor_cell.data_v2$disease,
         levels = c("normal","Pick disease"))

str(donor_cell.data_factor$disease)

#Creation of the mean expression matrix

# Metadata columns
metadata_cols <- c(
  "donor_id",
  "disease",
  "cell_type",
  "sex",
  "development_stage",
  "self_reported_ethnicity",
  "n_obs_aggregated"
)

# Identify all 2,000 gene columns
gene_cols <- setdiff(names(donor_cell.data_factor), metadata_cols)

# Divide every gene's summed expression by its corresponding cell count
mean_df <- donor_cell.data_factor %>%
  mutate(
    across(
      all_of(gene_cols),
      ~ .x / n_obs_aggregated
    )
  )
803750
head(mean_df)
dim(mean_df)

write.csv(donor_cell.data_factor, "sum_ct.csv")
write.csv(mean_df, "mean_ct.csv")

