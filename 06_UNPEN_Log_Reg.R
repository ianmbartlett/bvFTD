#### STEP 6 - LOGISTIC WITH DIFFEXP
library(readr)
library(dplyr)
library(pROC)

oligo_df <- read_csv("oligo_df_final.csv")
oligo_df$disease <- ifelse(oligo_df$disease == "Pick disease",1,0)


# Function for LOOCV logistic regression

# Function for LOOCV logistic regression
loocv_kit <- function(df) {
  
  # Number of donors
  n <- nrow(df)
  
  # Store KIT coefficients from each fold
  coefficients <- numeric(n)
  intercepts <- numeric(n)
  
  # Empty data frame to store predictions
  results <- data.frame(
    LOFLD = character(n),
    predicted_probability = numeric(n),
    predicted_outcome = integer(n),
    actual_outcome = integer(n),
    match = integer(n)
  )
  
  # Loop through each donor
  for (i in 1:n) {
    
    # Leave one donor out
    train_data <- df[-i, ]
    test_data  <- df[i, , drop = FALSE]
    
    # Fit simple logistic regression
    model <- glm(
      disease ~ KIT,
      data = train_data,
      family = binomial(link = "logit")
    )
    
    # Save KIT coefficient
    coefficients[i] <- coef(model)["KIT"]
    
    # Save model intercept
    intercepts[i] <- coef(model)["(Intercept)"]
    
    # Predict probability for the left-out donor
    pred_prob <- predict(
      model,
      newdata = test_data,
      type = "response"
    )
    
    # Classify using 0.75 threshold
    pred_outcome <- ifelse(pred_prob >= 0.75, 1, 0)
    
    # Actual outcome
    actual <- test_data$disease
    
    # Check prediction
    correct <- as.integer(pred_outcome == actual)
    
    # Store results
    results$LOFLD[i] <- as.character(test_data$donor_id)
    results$predicted_probability[i] <- as.numeric(pred_prob)
    results$predicted_outcome[i] <- pred_outcome
    results$actual_outcome[i] <- actual
    results$match[i] <- correct
  }
  
  # Calculate mean KIT coefficient
  mean_KIT_coefficient <- mean(coefficients, na.rm = TRUE)
  mean_intercept <- mean(intercepts, na.rm = TRUE)
  
  # Return predictions, coefficients, and mean coefficient
  return(list(
    predictions = results,
    coefficients = coefficients,
    mean_KIT_coefficient = mean_KIT_coefficient,
    mean_intercept = mean_intercept
  ))
}

# Run LOOCV
loocv_results <- loocv_kit(oligo_df)

# View predictions
print(loocv_results$predictions)

# View all 19 coefficients
print(loocv_results$coefficients)

# View mean coefficient
print(loocv_results$mean_KIT_coefficient)


# MODEL ACCURACY FOR COMPARISON
model <- c("UP_LogReg")
accuracy <- mean(loocv_results$predictions$match) * 100
sensi_match <- ifelse((loocv_results$predictions$predicted_outcome == 1 & loocv_results$predictions$actual_outcome == 1) 
                == TRUE, 1,0)
sensi <- sum(sensi_match)/sum(ifelse((loocv_results$predictions$actual_outcome == 1),1,0)) * 100 %>% round(3)
speci_match <- ifelse((loocv_results$predictions$predicted_outcome == 0 & loocv_results$predictions$actual_outcome == 0) 
                      == TRUE, 1,0)
speci <- sum(speci_match)/sum(ifelse((loocv_results$predictions$actual_outcome == 0),1,0)) * 100 %>% round(3)

print(loocv_results$mean_intercept)

#AUC Calculation
roc_model <- roc(
  response = loocv_results$predictions$actual_outcome,
  predictor = loocv_results$predictions$predicted_probability,
  levels = c(0, 1),
  direction = "<"
)

auc <- (auc(roc_model) %>% round(3))*100

bdfa <- (.3*accuracy)+(.1*sensi)+(.15*speci)+(auc*.45)


model_comparison <- data.frame(model,accuracy,sensi,speci,AUC=auc,BDFA=bdfa)

write.csv(model_comparison, "final.model_comparison.csv", row.names = FALSE)


