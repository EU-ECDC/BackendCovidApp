# Function to predict n days ahead with rolling window and generate output for evaluating predictive power
nStepAhead <- function(model, t_start, country, t_end = NULL, type = "first", nparam = 1,
                       n = 28, nsim = 100, cumu_inc = TRUE, delay = 30, t_start_cov = NULL,
                       keep.estimates = FALSE, actual = FALSE, transmissibility = 1, importation = 1,
                       default = TRUE, type_cov = NULL, delay_cov = 30, total = FALSE){
    if (cumu_inc & missing(delay)){
        warning("cumu_inc is TRUE but delay has not been set. Using default value for delay of 30.")
    }
    
    # Extract observed case counts from model object
    obs <- model$stsObj@observed
    # Get number of regions/region-age groups
    nb_reg <- ncol(obs)
    # Get last time point
    t_max <- nrow(obs)
    
    if (t_start > t_max){
        stop("Start time t_start for predictions must be within observed period.")
    }
    
    if (is.null(t_end)){
        t_end <- t_max - n
    }
    
    if (t_end < t_start){
        stop("End time t_end for predictions must be after start time t_start.")
    }
    
    if (t_end > t_max){
        stop("End time t_end for predictions must be within observed period.")
    }
    
    # If end of prediction period t_end is beyond end of observed data add NA 
    # rows to sts object so that future counts can be simulated
    if (t_end + n > t_max){
        sts_object <- model$stsObj
        model$stsObj <- sts(observed = rbind(obs, 
                                             matrix(NA, nrow = t_end + n - t_max, ncol = nb_reg)),
                            start = sts_object@start,
                            freq = sts_object@freq,
                            population = rbind(sts_object@populationFrac, 
                                               matrix(sts_object@populationFrac[t_max,], 
                                                      nrow = t_end + n - t_max, ncol = nb_reg, byrow = T)),
                            neighbourhood = sts_object@neighbourhood,
                            epoch = c(sts_object@epoch,sts_object@epoch[t_max] + (1:(t_end + n - t_max))))
    }
    max_lag <- model$max_lag
    
    # Predict covariates for prediction period and overwrite covariate data in model object
    if (is.null(t_start_cov)) t_start_cov <- t_start
    covariates_pred <- predict_covariates(model, t_start_cov, t_end + n, country = country, actual = actual, 
                                          default = default, type = type_cov, delay = delay_cov)
    nms <- names(model$control$data)
    for (i in seq_along(model$control$data)){
        nm <- nms[i]
        if (is.null(dim(model$control$data[[i]]))){
            if (t_end + n <= t_max){ 
                # if final predicted time point is within observed data, then 
                # overwrite covariate data for prediction period
                model$control$data[[i]][(t_start_cov+1):(t_end+n)] <- covariates_pred[[nm]]
            } else { 
                # otherwise bind predicted covariate data to actual covariate data up to t_start_cov
                covariate <- model$control$data[[i]]
                model$control$data[[i]] <- c(covariate[1:t_start_cov],covariates_pred[[nm]])
            }
        } else {
            if (t_end + n <= t_max){
                model$control$data[[i]][(t_start_cov+1):(t_end+n),] <- covariates_pred[[nm]]    
            } else {
                covariate <- model$control$data[[i]]
                model$control$data[[i]] <- rbind(covariate[1:t_start_cov,],covariates_pred[[nm]])
            }
        }
    }
    
    # Extract model terms from model object
    term <- model$terms
    if (is.null(term)){
        term <- model$terms <- terms(model)
    }
    
    # Get first fitted time point
    t_min <- model$control$subset[1]
    
    # Make vector of predicted time points
    tps <- (t_start:t_end) + n
    ntps <- length(tps)
    
    # Initialise matrices for storing predictions and dispersion parameter values, 
    # and variables for evaluating predictive performance if keep.estimates = TRUE
    pred <- matrix(NA, nrow = ntps, ncol = nb_reg, 
                   dimnames = list(tps, colnames(obs)))
    psi <- matrix(NA, nrow = ntps, ncol = 1)
    if (keep.estimates) {
        coefficients <- matrix(NA, nrow = ntps, ncol = length(term$initialTheta), 
                               dimnames = list(tps, names(term$initialTheta)))
        Sigma.orig <- matrix(NA, nrow = ntps, ncol = term$nSigma, 
                             dimnames = list(tps, names(model$Sigma.orig)))
        logliks <- matrix(NA, nrow = ntps, ncol = 2L, 
                          dimnames = list(tps, c("loglikelihood", "margll")))
    }
    
    # Refit model up to starting day
    # N.B. Currently not estimating lag weights so set refit_par_lag to false
    mod <- update(model, refit_par_lag = F, subset = t_min:t_start)
    if(!mod$convergence) stop("Model did not converge")
    for(i in 1:ntps){
        # Time point up to which model is refitted
        t_fit <- t_start + i - 1
        # Time point to be predicted
        t_pred <- t_start + i - 1 + n
        # Print current prediction window
        window <- c(t_fit, t_pred)
        print(window)
        
        # If making rolling n-step-ahead predictions refit model up to the current day
        if (type=="rolling" & i>1){ 
            mod <- update(model, refit_par_lag = F, subset = t_min:t_fit)
        }
        if (nparam > 1){
            # If nparam > 1, sample nparam parameter sets from multivariate normal 
            # distribution using parameter covariance matrix
            coef_ref <- mod$coefficients
            cov_mat <- mod$cov
            # Compute the different parameter sets
            matrix_params <- coef_ref + (cov_mat %>% chol %>% t) %*%
                matrix(rnorm(n = ncol(cov_mat) * nparam), nrow = ncol(cov_mat), ncol = nparam)
            matrix_params[mod$se > 10, ] <- coef_ref[mod$se > 10]
            list_sim <- list()
            for(j in seq_len(nparam)){
                model_j <- mod
                # Set the coefficients used for model_j
                model_j$coefficients <- matrix_params[,j]
                # Set the offsets for the epidemic and endemic components
                model_j$control$ne$offset <- transmissibility
                model_j$control$end$offset <- importation
                
                # Simulate n days ahead of current fit
                sim_j <- simulate(model_j, nsim = nsim,
                                  subset = t_fit + (1:n),
                                  y.start = obs[t_fit - ((max_lag - 1):0),],
                                  simplify = T,
                                  cumu_inc = cumu_inc,
                                  delay = delay,
                                  total = total)
                list_sim[[j]] <- sim_j
            }
            # Convert list of hhh4sims objects to hhh4sims object
            sim <- list_to_sims(list_sim = list_sim)
        } else{
            # Set the offsets for the epidemic and endemic components
            mod$control$ne$offset <- transmissibility
            mod$control$end$offset <- importation
            
            # Simulate n days ahead of current fit
            sim <- simulate(mod, nsim = nsim,
                            subset = t_fit + (1:n),
                            y.start = obs[t_fit - ((max_lag - 1):0),],
                            simplify = T,
                            cumu_inc = cumu_inc,
                            delay = delay,
                            total = total)
        }    
        # Calculate predicted mean from simulations
        mu <- apply(sim, c(1,2), mean)[n,]
        # Extract dispersion parameter from model fit
        disp <- coef(mod)["overdisp"]
        
        # Store the predicted means and dispersion
        pred[i,] <- mu
        psi[i,] <- disp
        
        if (keep.estimates){
            coefficients[i, ] <- mod$coefficients
            Sigma.orig[i, ] <- unname(mod$Sigma.orig)
            logliks[i, ] <- c(mod$loglikelihood,mod$margll)
        }
    }
    
    res <- c(list(sim = sim, pred = pred, observed = model$stsObj@observed[tps,,drop = F], psi = psi),
             if (keep.estimates) list(coefficients = coefficients, 
                                      Sigma.orig = Sigma.orig,
                                      logliks = logliks))
    # Make class of output "oneStepAhead" so oneStepAhead methods can be used on it
    class(res) <- "oneStepAhead"
    
    return(res)
    
}

# Function to convert list of hhh4sims into hhh4sims
list_to_sims <- function(list_sim){
    sims_out <- do.call("abind", list_sim)
    attr(sims_out, "initial") <- attr(list_sim[[1]], "initial")
    attr(sims_out, "stsObserved") <- attr(list_sim[[1]], "stsObserved")
    attr(sims_out, "class") <- attr(list_sim[[1]], "class")
    return(sims_out)
}
