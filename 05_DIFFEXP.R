# Differential Analysis
### Goal: Store genes and predictors in components for models
library(readr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(edgeR)


mean_df <- readr::read_csv("mean_ct.csv")
sum_df <- readr::read_csv("sum_ct.csv")
sum_df$disease <- as.factor(sum_df$disease)
sum_df <- sum_df2


#TOP 5 Cell types with the most observations

top5_celltypes <- sum_df2 %>%
  group_by(cell_type) %>%
  summarise(
    total_cells = sum(n_obs_aggregated, na.rm = TRUE)
  ) %>%
  arrange(desc(total_cells)) %>%
  slice_head(n = 5)

top5_celltypes

## NUMBER 1: DIFEXP W/ OLIGODENDROCYTES

oligo_data <- sum_df2 %>% subset(cell_type == "oligodendrocyte")
oligo_meta <- oligo_data[, c(
  "donor_id",
  "disease",
  "sex",
  "development_stage",
  "n_obs_aggregated"
)]

oligo_counts <- oligo_data %>% select(-c(
  "donor_id",
  "cell_type",
  "disease",
  "sex",
  "development_stage",
  "n_obs_aggregated"
))

oligo_t_counts <- t(as.matrix(oligo_counts)) %>% round()

group <- factor(oligo_data$disease)

table(group)

# Create edgeR object
y <- DGEList(
  counts = oligo_t_counts,
  group = group
)

# Remove genes with insufficient expression
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize library sizes
y <- calcNormFactors(y)

# Set up the statistical model
design <- model.matrix(~ group)

# Estimate variability between donors
y <- estimateDisp(y, design)

# Fit the quasi-likelihood model
fit <- glmQLFit(y, design, robust = TRUE)

# Test Pick disease versus healthy
test <- glmQLFTest(fit, coef = 2)

results <- topTags(
  test,
  n = Inf,
  sort.by = "PValue"
)$table

head(results, 10)

#Genes that meet the FDR threshold
significant_genes <- results[
  results$FDR < 0.05,
]

dim(significant_genes)

# Rank significant genes by absolute log fold change
largest_changes <- significant_genes[
  order(abs(significant_genes$logFC), decreasing = TRUE),
]

head(largest_changes, 10)

## OLIGODENDROCYTE ANALYSIS = 113 SIGNIFICANT GENES

## NUMBER 2: DIFEXP W/ GLUTAMATERGIC NEURONS

gn_data <- sum_df2 %>% subset(cell_type == "glutamatergic neuron")
gn_meta <- gn_data[, c(
  "donor_id",
  "disease",
  "sex",
  "development_stage",
  "n_obs_aggregated"
)]

gn_counts <- gn_data %>% select(-c(
  "donor_id",
  "cell_type",
  "disease",
  "sex",
  "development_stage",
  "n_obs_aggregated"
))

gn_t_counts <- t(as.matrix(gn_counts)) %>% round()

group <- factor(gn_data$disease)

table(group)

# Create edgeR object
y <- DGEList(
  counts = gn_t_counts,
  group = group
)

# Remove genes with insufficient expression
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize library sizes
y <- calcNormFactors(y)

# Set up the statistical model
design <- model.matrix(~ group)

# Estimate variability between donors
y <- estimateDisp(y, design)

# Fit the quasi-likelihood model
fit <- glmQLFit(y, design, robust = TRUE)

# Test Pick disease versus healthy
test <- glmQLFTest(fit, coef = 2)

results <- topTags(
  test,
  n = Inf,
  sort.by = "PValue"
)$table

head(results, 10)

#Genes that meet the FDR threshold
significant_genes <- results[
  results$FDR < 0.05,
]

dim(significant_genes)

# Rank significant genes by absolute log fold change
largest_changes <- significant_genes[
  order(abs(significant_genes$logFC), decreasing = TRUE),
]

head(largest_changes, 10)

## GLUTMATERGIC NEURON ANALYSIS = 69 SIGNIFICANT GENES

## DIFFEXP WITH ASTROCYTES

astro_data <- sum_df2 %>% subset(cell_type == "astrocyte")
astro_meta <- astro_data[, c(
  "donor_id",
  "disease",
  "sex",
  "development_stage",
  "n_obs_aggregated"
)]

astro_counts <- astro_data %>% select(-c(
  "donor_id",
  "cell_type",
  "disease",
  "sex",
  "development_stage",
  "n_obs_aggregated"
))

counts <- t(as.matrix(astro_counts)) %>% round()

group <- factor(astro_data$disease)

table(group)

# Create edgeR object
y <- DGEList(
  counts = counts,
  group = group
)

# Remove genes with insufficient expression
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize library sizes
y <- calcNormFactors(y)

# Set up the statistical model
design <- model.matrix(~ group)

# Estimate variability between donors
y <- estimateDisp(y, design)

# Fit the quasi-likelihood model
fit <- glmQLFit(y, design, robust = TRUE)

# Test Pick disease versus healthy
test <- glmQLFTest(fit, coef = 2)

results <- topTags(
  test,
  n = Inf,
  sort.by = "PValue"
)$table

head(results, 10)

#Genes that meet the FDR threshold
significant_genes <- results[
  results$FDR < 0.05,
]

dim(significant_genes)

# Rank significant genes by absolute log fold change
largest_changes <- significant_genes[
  order(abs(significant_genes$logFC), decreasing = TRUE),
]

head(largest_changes, 10)

## ASTROCYTE ANALYSIS = 23 SIGNIFICANT GENES

### PROJECT WILL PROCEED WITH OLIGODENDROCYTES

### OLIGODENDROCYTE FILTERING AND TOP 10 SIGNIFICANT GENE FILTERING
view(oligo_data)
view(oligo_counts)
view(oligo_meta)

oligo_data_filt <- oligo_data %>% select(
  "ENSG00000121966",
  "ENSG00000149257",
  "ENSG00000189056",
  "ENSG00000157404",
  "ENSG00000162645",
  "ENSG00000206172",
  "ENSG00000244734",
  "ENSG00000188536",
  "ENSG00000106211",
  "ENSG00000111907"
)

colnames(oligo_data_filt) <- c(
  "CXCR4",
  "SERPINH1",
  "RELN",
  "KIT",
  "GBP2",
  "HBA1",
  "HBB",
  "HBA2",
  "HSPB1",
  "TPD52L1")
oligo_df_final <- data.frame(oligo_meta,oligo_data_filt)

head(oligo_df_final)

write.csv(oligo_df_final, "oligo_df_final.csv")




