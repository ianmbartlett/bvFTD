#Read in oligo data frame
library(readr)
library(glmnet)
library(dplyr)
library(pander)
oligo_df <- readr::read_csv("oligo_df_final.csv")
gene_cols <- names(oligo_df[,7:16])
oligo_df$disease <- ifelse(oligo_df$disease == "Pick disease",1,0)


loocv_lasso_bic <- function(
    donor_data,
    gene_cols = NULL,
    id_col = "donor_id",
    outcome_col = "disease",
    threshold = 0.75
) {
  
  # If gene_cols is not supplied, use numeric columns
  # other than the donor ID and disease outcome.
  if (is.null(gene_cols)) {
    gene_cols <- names(oligo_df)[
      vapply(oligo_df, is.numeric, logical(1)) &
        !names(oligo_df) %in% c(id_col,outcome_col)
    ]
  }
  
  # Extract gene-expression predictors and outcome
  x <- as.matrix(donor_data[, gene_cols, drop = FALSE])
  y <- donor_data[[outcome_col]]
  
  # Basic checks
  if (!all(y %in% c(0, 1)) || length(unique(y)) != 2) {
    stop("The disease outcome must contain both 0 and 1.")
  }
  
  if (anyNA(x) || anyNA(y)) {
    stop("Handle missing values before running LOOCV.")
  }
  
  if (!all(is.finite(x))) {
    stop("Gene-expression predictors must be finite.")
  }
  
  if (nrow(x) != nrow(donor_data)) {
    stop("Predictor and outcome row counts do not match.")
  }
  
  # Store one prediction per donor
  n <- nrow(donor_data)
  
  results <- data.frame(
    LOFLD = as.character(donor_data[[id_col]]),
    predicted_probability = rep(NA_real_, n),
    predicted_outcome = rep(NA_integer_, n),
    actual_outcome = as.integer(y),
    match = rep(NA_integer_, n),
    best_lambda = rep(NA_real_, n),
    BIC = rep(NA_real_, n)
  )
  
  # Store selected genes and model objects by fold
  selected_genes <- vector("list", n)
  fitted_models <- vector("list", n)
  
  # LOOCV: leave one donor out at a time
  for (i in seq_len(n)) {
    
    # Training and test data
    train_x <- x[-i, , drop = FALSE]
    test_x  <- x[i, , drop = FALSE]
    
    train_y <- y[-i]
    test_y  <- y[i]
    
    # Both classes must remain in the training fold
    if (length(unique(train_y)) < 2) {
      stop(
        "A training fold contains only one outcome class. ",
        "Check the number of disease and normal donors."
      )
    }
    
    # Fit LASSO logistic regression
    lasso_fit <- glmnet(
      x = train_x,
      y = train_y,
      family = "binomial",
      alpha = 1,
      standardize = TRUE
    )
    
    # Candidate lambda values generated using training data
    lambda_values <- lasso_fit$lambda
    
    # Calculate BIC for every candidate lambda
    bic_values <- sapply(lambda_values, function(lam) {
      
      # Training probabilities at this lambda
      p <- as.numeric(predict(
        lasso_fit,
        newx = train_x,
        s = lam,
        type = "response"
      ))
      
      # Protect against log(0)
      p <- pmin(pmax(p, 1e-15), 1 - 1e-15)
      
      # Training log-likelihood
      log_likelihood <- sum(
        train_y * log(p) +
          (1 - train_y) * log(1 - p)
      )
      
      # Count nonzero coefficients, including intercept
      beta <- as.matrix(coef(lasso_fit, s = lam))
      k <- sum(beta[, 1] != 0)
      
      # Bayesian Information Criterion
      -2 * log_likelihood + k * log(length(train_y))
    })
    
    # Lambda with minimum BIC in this fold
    best_index <- which.min(bic_values)
    best_lambda <- lambda_values[best_index]
    best_bic <- bic_values[best_index]
    
    # Predict the held-out donor
    pred_prob <- as.numeric(predict(
      lasso_fit,
      newx = test_x,
      s = best_lambda,
      type = "response"
    ))
    
    # Apply the predetermined 0.75 threshold
    pred_outcome <- as.integer(pred_prob >= threshold)
    
    # Record prediction results
    results$predicted_probability[i] <- pred_prob
    results$predicted_outcome[i] <- pred_outcome
    results$match[i] <- as.integer(pred_outcome == test_y)
    results$best_lambda[i] <- best_lambda
    results$BIC[i] <- best_bic
    
    # Extract selected genes and their coefficients
    beta <- as.matrix(coef(lasso_fit, s = best_lambda))
    
    selected <- beta[-1, 1] != 0
    
    selected_genes[[i]] <- data.frame(
      LOFLD = as.character(donor_data[[id_col]][i]),
      gene = rownames(beta)[-1][selected],
      coefficient = beta[-1, 1][selected]
    )
    
    # Save model for later inspection
    fitted_models[[i]] <- lasso_fit
  }
  
  # Combine selected genes from all folds
  selected_genes <- do.call(rbind, selected_genes)
  rownames(selected_genes) <- NULL
  
  # Return all outputs
  return(list(
    results = results,
    selected_genes = selected_genes,
    models = fitted_models
  ))
}

lasso_loocv <- loocv_lasso_bic(donor_data=oligo_df,
                               gene_cols = gene_cols)

View(lasso_loocv$results)

# The Penalized logistic regression automatically selects genes and standardizes (log1p) them.
# This table shows us the most frequent genes that were assumed significant by the models
# upon each re-simulation. RELN was included in all 19 iterations of the models.
# This does not indicated RELN causes Pick disease, but rather that it would be a strong candidate
# for further analysis in determining if it does.

table(lasso_loocv$selected_genes$gene)

#Adding this to our model comparison df
#AUC 
roc_model <- roc(
  response = lasso_loocv$results$actual_outcome,
  predictor = lasso_loocv$results$predicted_probability,
  levels = c(0, 1),
  direction = "<"
)
auc <- (auc(roc_model) %>% round(3))*100

accuracy <- round(mean(lasso_loocv$results$match) *100 ,3)
sensi <- (6/9)*100
speci <- (8/10)*100

bdfa <- (.3*accuracy)+(.1*sensi)+(.15*speci)+(auc*.45)

model_comparison <- read_csv("final.model_comparison.csv") 
model_comparison <- model_comparison %>%
  add_row(
    model = "LASSO_LR",
    accuracy = accuracy,
    sensi = (6/9) * 100,
    speci = (8/10) * 100,
    AUC = auc,
    BDFA = bdfa
  )



write_csv(model_comparison,"step2_modelcomp.csv")




