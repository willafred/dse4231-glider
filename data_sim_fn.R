
library(tidyverse)

# Overall function to generate dataframe

df_generator = function(n, p, sigma, treatment_eqn, outcome_eqn) {
  # Generate covariates
  Xorig = matrix(rnorm(n * p, mean = 1, sd = sigma), n, p)
  X_df = as.data.frame(Xorig)
  colnames(X_df) = paste0("X", 1:ncol(X_df))
  
  # Compute treatment probabilities
  P = treatment_eqn(X_df)
  
  # Generate treatment assignment
  treatment = rbinom(n, 1, P)
  while(length(unique(treatment)) < 2 | sum(treatment == 0) < 5 | sum(treatment == 1) < 5) {
    treatment = rbinom(n, 1, P)  # Regenerate treatment
  }
  
  # Compute outcomes
  outcome = outcome_eqn(X_df, treatment, sigma)
  
  # Create final dataframe
  df = X_df %>%
    mutate(treatment = treatment, outcome = outcome)
  
  final_df = df %>%
    dplyr::select(outcome, treatment, starts_with("X"))  # Ensure correct column order
  return(final_df)
}

overall_scenario_generator = function(num_params, num_params_scene10){
  # --- ## SCENARIO 1 ## ---
  treatment_eqn = function(X) {
    linear_comb = 0.4*X$X1 + 0.3*X$X2 + 0.2*X$X5 +0.1*X$X4
    P = exp(linear_comb) / (1 + exp(linear_comb))  # Logistic function
    return(P)
  }
  
  outcome_eqn = function(X, A, sigma) {
    g_V = 0
    Y = A + g_V + rnorm(nrow(X), mean = 0, sd = sigma)  # Adding normal noise
    return(Y)
  }
  df_scenario_1 = df_generator(n = 500, p = num_params, sigma = 1, treatment_eqn, outcome_eqn)
  
  # --- ## SCENARIO 2 ## ---
  treatment_eqn = function(X) {
    linear_comb = 0.5*X$X1 + 0.5*X$X2 + 0.5*X$X3 +0.1*X$X4
    P = exp(linear_comb) / (1 + exp(linear_comb))  # Logistic function
    return(P)
  }
  
  outcome_eqn = function(X, A, sigma) {
    g_V = 0.5*X$X1 + X$X3 + 0.5*X$X4
    Y = A + g_V + rnorm(nrow(X), mean = 0, sd = sigma)  # Adding normal noise
    return(Y)
  }
  
  df_scenario_2 = df_generator(n = 500, p = num_params, sigma = 1, treatment_eqn, outcome_eqn)
  
  # --- ## SCENARIO 3 ## ---
  
  treatment_eqn = function(X) {
    linear_comb = 0.1*X$X1 + 0.1*X$X2 + 1*X$X3 + 1*X$X4 + 1*X$X5
    P = exp(linear_comb) / (1 + exp(linear_comb))  # Logistic function
    return(P)
  }
  
  outcome_eqn = function(X, A, sigma) {
    g_V = 2*X$X1 + 2*X$X3
    Y = A + g_V + rnorm(nrow(X), mean = 0, sd = sigma)  # Adding normal noise
    return(Y)
  }
  
  df_scenario_3 = df_generator(n = 500, p = num_params, sigma = 1, treatment_eqn, outcome_eqn)
  
  # --- ## SCENARIO 4 ## ---
  treatment_eqn = function(X) {
    linear_comb = 0.5*X$X1 + 0.4*X$X2 + 0.3*X$X3 + 0.2*X$X4 + 0.1*X$X5
    P = exp(linear_comb) / (1 + exp(linear_comb))  # Logistic function
    return(P)
  }
  
  outcome_eqn = function(X, A, sigma) {
    g_V = 0.5*X$X1 + 1*X$X2 + 1.5*X$X3 + 2*X$X4 + 2.5*X$X5
    Y = A + g_V + rnorm(nrow(X), mean = 0, sd = sigma)  # Adding normal noise
    return(Y)
  }
  
  df_scenario_4 = df_generator(n = 500, p = num_params, sigma = 1, treatment_eqn, outcome_eqn)
  
  # --- ## SCENARIO 5 ## ---
  treatment_eqn = function(X) {
    linear_comb = 0.5*X$X1 + 0.5*X$X2 + 0.1*X$X3
    P = exp(linear_comb) / (1 + exp(linear_comb))  # Logistic function
    return(P)
  }
  
  outcome_eqn = function(X, A, sigma) {
    # Manually compute interaction terms for the first five variables
    interaction_terms = NULL
    for (i in 1:5) {
      for (j in i:5) {
        interaction_terms = cbind(interaction_terms, X[, i] * X[, j])  # Create pairwise interaction terms
      }
    }
    interaction_sum = rowSums(interaction_terms)  # Summing interaction terms row-wise
    g_V = 1*X$X3 + 1*X$X4 + interaction_sum
    Y = A + g_V + rnorm(nrow(X), mean = 0, sd = sigma)  # Adding normal noise
    return(Y)
  }
  
  df_scenario_5 = df_generator(n = 500, p = num_params, sigma = 1, treatment_eqn, outcome_eqn)
  
  # --- ## SCENARIO 6 ## ---
  treatment_eqn = function(X) {
    linear_comb = 1*X$X1 + 1*X$X2 + 1*X$X5
    P = exp(linear_comb) / (1 + exp(linear_comb))  # Logistic function
    return(P)
  }
  
  outcome_eqn = function(X, A, sigma) {
    interaction_terms = NULL
    for (i in 1:5) {
      for (j in i:5) {
        interaction_terms = cbind(interaction_terms, X[, i] * X[, j])  # Create pairwise interaction terms
      }
    }
    g_V = rowSums(interaction_terms)  # Summing interaction terms row-wise
    Y = A + g_V + rnorm(nrow(X), mean = 0, sd = sigma)  # Adding normal noise
    return(Y)
  }
  
  df_scenario_6 = df_generator(n = 500, p = num_params, sigma = 1, treatment_eqn, outcome_eqn)
  
  # --- ## SCENARIO 7 ## ---
  
  treatment_eqn = function(X) {
    linear_comb = 0.2*X$X1 + 0.2*X$X2 + 0.2*X$X5
    P = exp(linear_comb) / (1 + exp(linear_comb))  # Logistic function
    return(P)
  }
  
  outcome_eqn = function(X, A, sigma) {
    g_V = 0.25*X$X3 + (X$X1 + X$X2)^2 - ((X$X1)^2 - X$X3)^2 + ((X$X4)^2 - 0.5*X$X5)*(X$X3 - 0.5*X$X4)
    Y = A + g_V + rnorm(nrow(X), mean = 0, sd = sigma)  # Adding normal noise
    return(Y)
  }
  
  df_scenario_7 = df_generator(n = 500, p = num_params, sigma = 1, treatment_eqn, outcome_eqn)
  
  # --- ## SCENARIO 8 ## ---
  
  treatment_eqn = function(X) {
    interaction_terms = NULL
    for (i in 1:5) {
      for (j in i:5) {
        interaction_terms = cbind(interaction_terms, X[, i] * X[, j])  # Create pairwise interaction terms
      }
    }
    interaction_sum = rowSums(interaction_terms)
    linear_comb = 1*X$X3 + 1*X$X4 + 1*X$X5 + interaction_sum
    P = exp(linear_comb) / (1 + exp(linear_comb))  # Logistic function
    return(P)
  }
  
  outcome_eqn = function(X, A, sigma) {
    g_V = 0.5*X$X1 + 0.5*X$X2 + 0.1*X$X3
    Y = A + g_V + rnorm(nrow(X), mean = 0, sd = sigma)  # Adding normal noise
    return(Y)
  }
  
  df_scenario_8 = df_generator(n = 500, p = num_params, sigma = 1, treatment_eqn, outcome_eqn)
  
  # --- ## SCENARIO 9 ## ---
  treatment_eqn = function(X) {
    linear_comb = (1*X$X1 + 1*X$X2 + 0.5*X$X3)^2
    P = exp(linear_comb) / (1 + exp(linear_comb))  # Logistic function
    return(P)
  }
  
  outcome_eqn = function(X, A, sigma) {
    g_V = 0.5*X$X1 + 0.5*X$X3 + 0.5*X$X4
    Y = A + g_V + rnorm(nrow(X), mean = 0, sd = sigma)  # Adding normal noise
    return(Y)
  }
  
  df_scenario_9 = df_generator(n = 500, p = num_params, sigma = 1, treatment_eqn, outcome_eqn)
  
  # --- ## SCENARIO 10 ## ---
  treatment_eqn = function(X) {
    linear_comb = 0.2*X$X1 - 2*X$X2 + X$X5 - X$X6 + X$X7 - X$X8
    P = exp(linear_comb) / (1 + exp(linear_comb))  # Logistic function
    return(P)
  }
  
  outcome_eqn = function(X, A, sigma) {
    g_V = 2*X$X1 + 0.2*X$X2 + 5*X$X3 + 5*X$X4
    Y = A + g_V + rnorm(nrow(X), mean = 0, sd = sigma)  # Adding normal noise
    return(Y)
  }
  df_scenario_10 = df_generator(n = 500, p = num_params_scene10, sigma = 2, treatment_eqn, outcome_eqn)
  
  # --- ## SCENARIO 11 ## ---
  treatment_eqn = function(X) {
    linear_comb = 2*X$X1 + 1*X$X2 + 5*X$X3
    P = exp(linear_comb) / (1 + exp(linear_comb))  # Logistic function
    return(P)
  }
  
  outcome_eqn = function(X, A, sigma) {
    g_V = sin(0.5*X$X1 + 0.5*X$X3 + 0.5*X$X4)
    Y = A + g_V + rnorm(nrow(X), mean = 0, sd = sigma)  # Adding normal noise
    return(Y)
  }
  
  df_scenario_11 = df_generator(n = 500, p = num_params, sigma = 1, treatment_eqn, outcome_eqn)
  
  # --- ## SCENARIO 12 ## ---
  
  treatment_eqn = function(X) {
    linear_comb = cos(2*X$X1 + 1*X$X2 + 5*X$X3)
    P = exp(linear_comb) / (1 + exp(linear_comb))  # Logistic function
    return(P)
  }
  
  outcome_eqn = function(X, A, sigma) {
    g_V = 1*X$X1 + 1*X$X3 + 1*X$X4 - 5*X$X5
    Y = A + g_V + rnorm(nrow(X), mean = 0, sd = sigma)  # Adding normal noise
    return(Y)
  }
  
  df_scenario_12 = df_generator(n = 500, p = num_params, sigma = 1, treatment_eqn, outcome_eqn)
  
  # --- ## SCENARIO 13 ## ---
  df_generator_categorical= function(n, p, sigma, treatment_eqn, outcome_eqn, cat_levels = 3) {
    # Generate continuous covariates
    Xorig = matrix(rnorm(n * p, mean = 1, sd = sigma), n, p)
    
    # Generate categorical covariates
    cat_cols = sample(1:cat_levels, n * 2, replace = TRUE)  # Assume 2 categorical variables
    cat_matrix = matrix(cat_cols, nrow = n, ncol = 2)  # 2 categorical variables
    
    # Combine continuous and categorical covariates
    X_full = cbind(Xorig, cat_matrix)
    
    X_df = as.data.frame(X_full)
    colnames(X_df) = c(paste0("X", 1:ncol(Xorig)), paste0("Cat", 1:ncol(cat_matrix)))
    
    # Compute treatment probabilities
    P = treatment_eqn(X_df)
    
    # Generate treatment assignment
    treatment = rbinom(n, 1, P)
    
    # Compute outcomes
    outcome = outcome_eqn(X_df, treatment, sigma)
    
    # Create final dataframe
    df = X_df %>%
      mutate(treatment = treatment, outcome = outcome)
    
    final_df = df %>%
      dplyr::select(outcome, treatment, starts_with("X"), starts_with("Cat"))  # Ensure correct column order
    return(final_df)
  }
  
  treatment_eqn = function(X) {
    # Include the categorical variables in the treatment equation (e.g., as dummy variables)
    linear_comb = 2*X$X1 + 1*X$X2 + 5*X$X3
    P = exp(linear_comb) / (1 + exp(linear_comb))  # Logistic function
    return(P)
  }
  
  outcome_eqn = function(X, A, sigma) {
    g_V = 1*X$X1 + 1*X$X3 + 1*X$X4 + (X$X5)^2
    Y = A + g_V + rnorm(nrow(X), mean = 0, sd = sigma)  # Adding normal noise
    return(Y)
  }
  df_scenario_13 = df_generator_categorical(n = 500, p = num_params, sigma = 1, treatment_eqn, outcome_eqn)
  
  combined_scenarios = list(
    scenario_1 = df_scenario_1,
    scenario_2 = df_scenario_2,
    scenario_3 = df_scenario_3,
    scenario_4 = df_scenario_4,
    scenario_5 = df_scenario_5,
    scenario_6 = df_scenario_6,
    scenario_7 = df_scenario_7,
    scenario_8 = df_scenario_8,
    scenario_9 = df_scenario_9,
    scenario_10 = df_scenario_10,
    scenario_11 = df_scenario_11,
    scenario_12 = df_scenario_12,
    scenario_13 = df_scenario_13
  )
  return(combined_scenarios)
}