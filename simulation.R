# Function generates a correlation matrix corresponding to AR(1) dependence structure.
ar1_matrix <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n-1))
  rho^exponent
}

# Function generates an equal correlation matrix.
equal_cor_matrix <- function(n, rho) {
  matrix <- rho*matrix(rep(1, n ^ 2), nrow = n, ncol = n) + (1-rho)*diag(n)
  if (!is.positive.definite(matrix)) {
    # Use solve() to take the inverse of the equal_cor_matrix with opposite rho.
    matrix <- solve(equal_cor_matrix(n, -rho))
    # Rescale the matrix.
    matrix <- matrix / matrix[1,1]
  }
  matrix
}

# Fits linear regression model without intercept.
# Inputs: vars = a vector of strings representing selected variables; 
# data = data frame representing an input data and response values.
# Output: a matrix of fitted coefficients.
fit_linear_model <- function (vars, data) 
{
  # Create a vector of response values.
  y <- data.matrix(data$Y)
  # Initialize a design matrix X to 0.
  X <- matrix(, nrow = 50, ncol = 0)
  # Append columns corresponding to selected variables to X.
  for (var in vars) {X <- cbind(X, data[var])}
  # Use solve() to calculate an inverse.
  X <- data.matrix(X)
  betas <- solve(t(X) %*% X) %*% t(X) %*% y
  fitted_y <- X %*% betas
}

# Generates a vector of standard errors..
# An element at index i is a standard error of the ith data sample.
# Parameters: vars = vector of variables; data is a dataframe of data samples.
get_standard_errors <- function(vars, data, variance) {
  # Generate a matrix containing only selected variables.
  X <- matrix(, nrow = 50, ncol = 0)
  # Append columns corresponding to selected variables to X.
  for (var in vars) {X <- cbind(X, data[var])}
  X <- data.matrix(X)
  std_error_matrix <- (variance) * (X %*% solve(t(X) %*% X) %*% t(X))
  std_error <- c()
  # Take the entries on the diagonal.
  for (i in 1:dim(std_error_matrix)[1]) {
    std_error <- c(std_error, sqrt(std_error_matrix[i,i]))
  }
  std_error
}

# Returns the residual mean square.
estimate_variance <- function(vars, data) {
  # Find the fitted response value.
  fitted_y <- fit_linear_model(vars, data)
  # Compute residual sum of squares.
  res_sum <- sum((data$Y - fitted_y)^2)
  # Estimate (sigma^2) with residual mean square.
  variance = res_sum / (50 - length(vars))
}

# Finds the 95% two-sided confidence intervals for the given data set with selected variables.
# Inputs: vars = a vector of strings representing selected variables; 
# data = data frame representing an input data and response values;
# variance = "true" variance or -1 if is_real_sigma == FALSE;
# is_real_sigma = a boolean of whether we want to use a real sigma for CIs.
# Output: a matrix of lower and upper bounds for each data point in the experiment.
get_cis <- function (vars, data, is_real_sigma, variance) 
{
  # Fit the linear regression model for the given subset of variables.
  fitted_y <- fit_linear_model(vars, data)
  if (!is_real_sigma) {
    # Estimate (sigma^2) using residual mean square.
    variance = estimate_variance(vars, data)
    # Find a t-score corresponding to 2.5% of the upper tail.
    # Degree of freedom for selected model  = 50 - size.
    score <- qt(p = 0.025, df = 50 - length(vars), lower.tail = FALSE)
  } else {
    # Find a Z-score corresponding to 2.5% of the upper tail.
    score <- qnorm(p = 0.025, mean = 0, sd = 1, lower.tail=FALSE)
  }
  # Find the standard errors.
  se <- get_standard_errors(vars, data, variance)
  # Generate a matrix of CIs.
  cis <- matrix(, nrow = dim(data)[1], ncol = 0) # Create an empty matrix.
  # Put the column of lower bounds in the matrix.
  cis <- cbind(cis, fitted_y - score * se)
  
  # Put the column of upper bounds in the matrix.
  cis <- cbind(cis, fitted_y + score * se)
  # Name columns.
  colnames(cis) <- c("lwr", "upr")
  cis
}
# Determines whether the mean response value lies in the given interval.
in_interval <- function(interval, mean_y, i) {
  return((mean_y[i]>=interval[i,"lwr"]) & (mean_y[i]<=interval[i,"upr"]))
}

# Returns TRUE if the selected model overfits, FALSE otherwise.
strict_overfit <- function(vars, oracle_vars) {
  all(oracle_vars %in% vars) & (length(vars) > length(oracle_vars))
}

# Returns TRUE if the selected model contains all true variables, FALSE otherwise.
contains_true_vars <- function(vars, oracle_vars) {
  all(oracle_vars %in% vars)
}

# Function simulates N datasets.
# Required: to install packages "MASS", "leaps", "MuMIn", "ggplot2", "ggpubr", "rgl", "matlib", "lqmm".
# N = number of simulations
# B = a vector of coefficients in front of X1, X2, X3 (B = c(1, 2, 3) in our case)

# statistic = one of the following strings:
# "aic" = Akaike Information Criteria 
# "aicc" = Akaike Information Criteria Corrected
# "bic" = Schwarz Bayesian Information Criteria 
# "adjr2" = adjusted R^2
# "cp" = Mallow's Cp

# scenario is an integer corresponding to the following:
# 1: AR(1) with rho = 0.5
# 2: AR(1) with rho = -0.5
# 3: IID error
# 4: equal correlation with rho = 0.5 
# 5: equal correlation with rho = -0.5 
simulate <- function(N, B, statistic, scenario){
  library(MASS)
  library(leaps)
  library(MuMIn)
  library(ggplot2)
  library(ggpubr)
  options(rgl.useNULL = TRUE)
  library(rgl)
  library(matlib)
  library(lqmm) # must cite this package: call citation("lqmm")
  
  # Correlation matrix
  cor_matrix <- switch(scenario, ar1_matrix(10, 0.5), ar1_matrix(10, -0.5), 
                       diag(10), equal_cor_matrix(10, 0.5), equal_cor_matrix(10, -0.5))
  if (!is.positive.definite(cor_matrix)) {
    return("ERROR: a given correlation matrix is not positive definite!")
  }
  # A list of all possible variables.
  all_vars <- c(paste("X", 1:10, sep=""))
  oracle_vars <- c()
  for (i in 1:10) {
    if (B[i] != 0) {
      oracle_vars <- c(oracle_vars, paste("X", i, sep=""))
    }
  }
  # Initialize the arrays of sigmas of the oracle model and the estimated model. 
  real_sigma <- c()
  est_sigma <- c()
  with_true_vars = 0
  overfitted = 0
  oracle_ci_coverage = 0
  est_ci_coverage = 0
  oracle_ci_true_coverage = 0
  est_ci_true_coverage = 0
  # Repeat the experiment N times
  for(j in 1:N) {
    # Generate data set.
    df <- data.frame(mvrnorm(n = 50, mu = rep(0, 10), Sigma = cor_matrix, empirical = TRUE))
    # Mean-centering input data.
    data_matrix <- sapply(df, function(x) scale(x, scale=FALSE))
    df <- as.data.frame(data_matrix)
    # Generate sample response values.
    y <- data_matrix %*% B + rnorm(50, 0, 1)
    # Mean-centering response values to get rid of intercept.
    y <- scale(y, scale = FALSE)
    # Generate fitted response values and mean-center them.
    mean_y <- data_matrix %*% B
    mean_y <- scale(mean_y, scale = FALSE)
    df$Y <- y
    
    # Choose the best subset of variables.
    # Since we centered the data, we can now disregard an intercept.
    models <- regsubsets(Y~., data = df, nvmax = 10, intercept=FALSE)
    reg_summary = summary(models)
    
    # Choose the most optimal size of the subset of variables.
    if (statistic == 'bic') {
      size = which.min(reg_summary$bic)
    } else if (statistic == 'cp') {
      size = which.min(reg_summary$cp)
    } else if (statistic == 'adjr2') {
      size = which.max(reg_summary$adjr2)
    } else if ((statistic == 'aic') | (statistic == 'aicc')) {
      aics <- c()
      for (i in 1:10) {
        vars <- all_vars[reg_summary$which[i,]]
        # Fit the linear regression with the variables in subset.
        fmla <- as.formula(paste("Y~", paste(vars, collapse="+")))
        # Remove the intercept: add "-1" to formula.
        fmla <- update(fmla, ~ . -1)
        est_fit <- lm(fmla, data = df)
        if(statistic == 'aic') {
          aics <- c(aics, AIC(est_fit))
        } else {
          aics <- c(aics, AICc(est_fit))
        }
      }
      size = which.min(aics)
    } else {
      print("Invalid statistic!")
      break;
    }
    
    # Create an array of selected variables for the most optimal size.
    vars <- all_vars[reg_summary$which[size,]]
    
    # Get 95% CIs for oracle model using estimated variance.
    est_cis_for_oracle <- get_cis(oracle_vars, df, is_real_sigma = FALSE, variance = -1)
    # Get 95% CIs for oracle model using true variance.
    real_cis_for_oracle <- get_cis(oracle_vars, df, is_real_sigma = TRUE, variance = 1) 
    
    # Get 95% CIs for selected model using estimated variance.
    est_cis_for_selected <- get_cis(vars, df, is_real_sigma = FALSE, variance = -1) 
    # Get 95% CIs for selected model using true variance.
    real_cis_for_selected <- get_cis(vars, df, is_real_sigma = TRUE, variance = 1)
    
    # Check a coverage of confidence intervals.
    for (i in 1:50) {
      # For oracle model with estimated variance.
      if (in_interval(est_cis_for_oracle, mean_y, i)) {
        oracle_ci_coverage = oracle_ci_coverage + 1
      }
      # For oracle model with real (true) variance.
      if (in_interval(real_cis_for_oracle, mean_y, i)) {
        oracle_ci_true_coverage = oracle_ci_true_coverage + 1
      }
      # For selected model with estimated variance.
      if (in_interval(est_cis_for_selected, mean_y, i)) {
        est_ci_coverage = est_ci_coverage + 1
      }
      # For selected model with real (true) variance.
      if (in_interval(real_cis_for_selected, mean_y, i)) {
        est_ci_true_coverage = est_ci_true_coverage + 1
      }
    }
    
    # Count the number of cases with non-strict overfitting:
    if(contains_true_vars(vars, oracle_vars)) {
      with_true_vars <- with_true_vars + 1
    }
    
    # Count the number of cases with strict overfitting:
    if(strict_overfit(vars, oracle_vars)) {
      overfitted <- overfitted + 1
      
      # Add the sigma of the selected model to the array.
      variance_selected = estimate_variance(vars, df)
      est_sigma <- c(est_sigma, sqrt(variance_selected))
      # Add the sigma of the oracle model to the array.
      variance_oracle = estimate_variance(oracle_vars, df)
      real_sigma <- c(real_sigma, sqrt(variance_oracle))
    }
  }
  
  # Plot the graph.
  plot(real_sigma, est_sigma, main="Scatterplot of oracle vs estimated model sigmas", xlab="Oracle model sigma", ylab="Estimated model sigma")
  
  sigma_ratios <- real_sigma/est_sigma
  
  hist(sigma_ratios, 
       main="Histogram of ratio of real sigma to estimated sigma", 
       xlab="Ratio[ratio > 1]", 
       ylab="Density",
       breaks = 15)
  print(paste("Number of simulations with non-strict overfitting =", with_true_vars))
  print(paste("Number of simulations with strict overfitting =", overfitted))
  print(paste("Mean of sigma ratios =", mean(sigma_ratios)))
  print(sprintf("Coverage of oracle model's CIs with estimated sigma is %f%%", (2*oracle_ci_coverage)/N))
  print(sprintf("Coverage of selected model's CIs with estimated sigma is %f%%", (2*est_ci_coverage)/N))
  
  print(sprintf("Coverage of oracle model's CIs with true sigma is %f%%", (2*oracle_ci_true_coverage)/N))
  print(sprintf("Coverage of selected model's CIs with true sigma is %f%%", (2*est_ci_true_coverage)/N))
  
}

# Call the function.
N <- 1000
# Vector of length 10.
B <- c(1, 2, 3, rep(0, 7))
B <- data.matrix(B)
# Possible statistics: "aic", "aicc", "bic", "adjr2", "cp".
statistic <- "aic"
# Possible scenarios: 1, 2, 3, 4, 5 (read lines 113-118 for explanation).
scenario <- 5
set.seed(7)
simulate(N, B, statistic, scenario)

