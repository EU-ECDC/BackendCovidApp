hhh4_cov <- function(sts_object, contact, equations, covariates, weight_fun,
                     lag_dist, parameter_lag, empty = FALSE, start = NULL){
  # Extract formula for each component
  formula_ar <- equations[["ar"]]
  formula_end <- equations[["end"]]
  formula_ne <- equations[["ne"]]
  ## Check column names match
  if (length(covariates) > 0){
    col_cov <- colnames(covariates[[1]])
    test_column <- lapply(covariates, function(X) 
      if (!all(colnames(X) == col_cov)) stop("column names do not match"))
    if (!all(colnames(sts_object@observed) == col_cov))
      stop("column names do not match")
  }
  regs <- unique(sub("[.].*", "", colnames(sts_object@observed)))
  nb_reg <- length(regs)
  nb_date <- nrow(sts_object@observed)
  # When the auto-regressive and neighbourhood components are combined, we want 
  # the power law to act on o+1, where o is the neighbour order
  if (is.null(formula_ar)){
    # If the distance matrix is actual distance, then do not add +1 (i.e. the 
    # distance matrix is not the degree of connectivity)
    if (all(neighbourhood(sts_object) == round(neighbourhood(sts_object))) &
       max(neighbourhood(sts_object) < 20))
      neighbourhood(sts_object) <- neighbourhood(sts_object) + 1
    formula_ar <- formula(~ -1)
  } else
    if (any(neighbourhood(sts_object) != round(neighbourhood(sts_object))) ||
       max(neighbourhood(sts_object) > 20))
      neighbourhood(sts_object)[neighbourhood(sts_object) == 0] <- -1
  

  # Set maximum lag for distributed lags
  max_lag <- 20L

  if (empty){
    sts_object@neighbourhood[sts_object@neighbourhood > 1] <- 0
    weight_fun <- sts_object@neighbourhood
  }
  
  # Fit a gravity model (power-law decay in neighbour order and population-
  # dependence in coefficient) for the neighbour weights with the given contact 
  # matrix (without normalisation) 
  control_pop <- list(
    ar = list(f = formula_ar),
    end = list(f = formula_end),
    ne = list(f = formula_ne, 
              weights = weight_fun, scale = expandC(contact,nb_reg)),
    data = covariates, 
    funct_lag = lag_dist, par_lag = parameter_lag, max_lag = max_lag, 
    family = "NegBin1", subset = (max_lag+1):nb_date, verbose = TRUE,
    start = list(fixed = start)
  )
  
  # Fit the model
  m1 <- hhh4_lag(sts_object,control = control_pop)
  return(m1)
}

