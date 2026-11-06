# OUR FINAL MODEL FOR COMPARISON: LINEAR SUPPORT VECTOR MACHINE 
library(e1071)
library(readr)
library(dplyr)

oligo_df <- read_csv("oligo_df_final.csv")

# Identify gene-expression columns
gene_cols <- names(oligo_df[,7:16])
oligo_df$disease <- ifelse(oligo_df$disease == "Pick disease",1,0)


# Keep the outcome numeric
oligo_df$disease <- as.integer(oligo_df$disease)

# Check the dataset
table(oligo_df$disease)
length(gene_cols)

library(e1071)

set.seed(123)

n <- nrow(oligo_df)

svm_results <- data.frame(
  donor_id = oligo_df$donor_id,
  actual_outcome = oligo_df$disease,
  predicted_outcome = NA_integer_,
  predicted_probability = NA_real_,
  decision_score = NA_real_
)

for (i in seq_len(n)) {
  
  # Split into training and test donors
  train_data <- oligo_df[-i, , drop = FALSE]
  test_data  <- oligo_df[i, , drop = FALSE]
  
  # Remove genes with no variation in this training fold
  gene_sd <- sapply(
    train_data[, gene_cols, drop = FALSE],
    sd,
    na.rm = TRUE
  )
  
  usable_genes <- gene_cols[
    is.finite(gene_sd) & gene_sd > 0
  ]
  
  # Standardize using training-fold means and SDs
  train_x <- scale(
    as.matrix(train_data[, usable_genes, drop = FALSE])
  )
  
  train_means <- attr(train_x, "scaled:center")
  train_sds   <- attr(train_x, "scaled:scale")
  
  test_x <- scale(
    as.matrix(test_data[, usable_genes, drop = FALSE]),
    center = train_means,
    scale = train_sds
  )
  
  # Define outcome classes
  train_y <- factor(
    train_data$disease,
    levels = c(0, 1),
    labels = c("Normal", "Pick")
  )
  
  # Fit linear SVM
  svm_fit <- svm(
    x = train_x,
    y = train_y,
    kernel = "linear",
    cost = 1,
    scale = FALSE,
    decision.values = TRUE
  )
  
  # Predict the held-out donor
  svm_pred <- predict(
    svm_fit,
    newdata = test_x,
    decision.values = TRUE
  )
  
  # Save predicted class
  svm_results$predicted_outcome[i] <-
    ifelse(as.character(svm_pred) == "Pick", 1L, 0L)
  
  # Save decision score
  decision <- attr(svm_pred, "decision.values")
  
  svm_results$decision_score[i] <-
    as.numeric(decision[1, 1])
}

# View final results
print(svm_results)

#Add the match colummn
svm_results$match <- ifelse(svm_results$actual_outcome == svm_results$predicted_outcome, 1,0)

library(pROC)

roc_svm <- roc(
  response = svm_results$actual_outcome,
  predictor = svm_results$decision_score,
  levels = c(0, 1),
  direction = "auto"
)

auc_svm <- (as.numeric(auc(roc_svm)) %>% round(3))*100

#Sensitivity
sensi_match <- ifelse(svm_results$actual_outcome == 1 & svm_results$predicted_outcome == 1, 1,0)
sensi <- sum(sensi_match)/sum(ifelse((svm_results$actual_outcome == 1),1,0)) * 100 %>% round(3)

#Specificity
speci <- speci_match <- ifelse((svm_results$predicted_outcome == 0 & svm_results$actual_outcome == 0) 
                               == TRUE, 1,0)
speci <- sum(speci_match)/sum(ifelse((svm_results$actual_outcome == 0),1,0)) * 100 %>% round(3)

bdfa <- (.3*(round(mean(svm_results$match) *100 ,3)))+(.1*sensi)+(.15*speci)+(auc_svm*.45)

mc <- read_csv("step2_modelcomp.csv")
mc <- mc %>% add_row(
  model = "SVM",
  accuracy = round(mean(svm_results$match) *100 ,3),
  sensi = sensi,
  speci = speci,
  AUC = auc_svm,
  BDFA = bdfa
)
print(mc)



