## Custom version of simulate.hhh4lag function from hhh4addon package to enable 
## simulation of distributed-lag age-stratified models with separate 
## autoregressive and neighbourhood components and to account for cumulative 
## incidence being covariate in autoregressive and neighbourhood predictors
my.simulate.hhh4lag <- function(object, # result from a call to hhh4lag
                                nsim = 1, # number of replicates to simulate
                                seed = NULL, 
                                y.start = NULL, # initial counts for epidemic components
                                subset = 1:nrow(object$stsObj), 
                                coefs = coef(object), # coefficients used for simulation
                                components = c("ar", "ne", "end"), # which comp to include
                                simplify = nsim > 1, # counts array only (no full sts)
                                cumu_inc = FALSE, # calculate cumulative incidence as simulation progresses
                                delay = 30, # delay between date cumulative incidence is calculated up to and current date
                                total = FALSE,...) 
{
    ## Determine seed (this part is copied from stats:::simulate.lm with
    ## Copyright (C) 1995-2012 The R Core Team)
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1) # initialize the RNG if necessary
    if (is.null(seed)){
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)        
    } else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    ## END seed
    
    cl <- match.call()
    theta <- if (missing(coefs)) 
        coefs
    else surveillance:::checkCoefs(object, coefs)
    control <- object$control
    
    # Get lags
    lag.ar <- object$control$ar$lag
    lag.ne <- object$control$ne$lag
    maxlag <- max(c(control$max_lag, control$ar$lag, control$ne$lag), 
                  na.rm = TRUE)
    minlag <- min(c(control$min_lag, control$ar$lag, control$ne$lag), 
                  na.rm = TRUE)
    
    # Get initial counts
    nUnits <- object$nUnit
    if (is.null(y.start)) {
        y.means <- ceiling(colMeans(observed(object$stsObj)[subset, 
                                                            , drop = FALSE]))
        y.start <- matrix(y.means, maxlag, nUnits, byrow = TRUE)
    } else {
        if (is.vector(y.start)) 
            y.start <- t(y.start)
        if (ncol(y.start) != nUnits) 
            stop(sQuote("y.start"), " must have nUnits=", nUnits, 
                 " columns")
        if (nrow(y.start) != maxlag) 
            stop("need 'y.start' values for lag=", maxlag, " initial time points (provide a matrix with ", 
                 maxlag, " rows).")
    }
    
    # Extract overdispersion parameters (simHHH4 assumes psi->0 means Poisson)
    model <- hhh4addon:::terms.hhh4lag(object)
    psi <- surveillance:::splitParams(theta, model)$overdisp
    if (length(psi) > 1) 
        psi <- psi[model$indexPsi]
    neweights <- surveillance:::getNEweights(object, coefW(theta))
    if (cumu_inc){
        # If cumulative incidence is a covariate in the model, calculate 
        # cumulative incidence and update covariates as simulation progresses
        if (total){
            pop <- population(object$stsObj)
        } else {
            pop <- object$control$data$pop * object$control$data$pop_age * 1e5    
        }
        # Get initial value cumulative incidence in last month (inc_new)
        init_inc_new <- object$control$data$inc_new[subset[1], ]
        # Get initial value of cumulative incidence between 6 months and 1 month ago (inc_old)
        init_inc_old <- object$control$data$inc_old[subset[1], ]
        # Initialise values of inc_new and inc_old covariates in object$control
        object$control$data$inc_new[] <- 0
        object$control$data$inc_old[] <- 0
        # Get values of predictors (lambda_ait, phi_ait, nu_ait) with offsets
        exppreds <- hhh4addon:::get_exppreds_with_offsets_lag(object, subset = subset,
                                                              theta = theta)
        object$control$data$inc_new[] <-
            matrix(init_inc_new, byrow = T, nrow = nrow(object$control$data$inc_new), 
                   ncol = length(init_inc_new))
        object$control$data$inc_old[] <-
            matrix(init_inc_old, byrow = T, nrow = nrow(object$control$data$inc_old), 
                   ncol = length(init_inc_old))
        # Call simulation function
        simcall <- quote(my.simHHH4lag_cumu_inc(psi = psi, neW = neweights, start = y.start, lag.ar = lag.ar, 
                                                lag.ne = lag.ne, funct_lag = control$funct_lag, par_lag = control$par_lag, 
                                                min_lag = control$min_lag, max_lag = control$max_lag, scale = control$ne$scale, 
                                                object = object, subset = subset, theta = theta, components = components, pop = pop,
                                                delay = delay, exppreds = exppreds))
    } else {
        # Get values of predictors with offsets
        exppreds <- hhh4addon:::get_exppreds_with_offsets_lag(object, 
                                                              subset = subset, 
                                                              theta = theta)
        ## Set predictor to zero if not included ('components' argument)
        stopifnot(length(components) > 0, components %in% c("ar", "ne", "end"))
        # Extract autoregressive, neighbourhood, and endemic predictors
        ar <- getComp("ar", exppreds, components)
        ne <- getComp("ne", exppreds, components)
        end <- getComp("end", exppreds, components)
        # Call simulation function
        simcall <- quote(my.simHHH4lag(ar = ar, ne = ne, end = end, 
                                       psi = psi, neW = neweights, start = y.start, lag.ar = lag.ar, 
                                       lag.ne = lag.ne, funct_lag = control$funct_lag, par_lag = control$par_lag, 
                                       min_lag = control$min_lag, max_lag = control$max_lag, scale = control$ne$scale))        
    }
    if (!simplify) {
        # If only 1 replicate, add observed counts to output object
        res0 <- object$stsObj[subset, ]
        setObserved <- function(observed) {
            res0@observed[] <- observed
            res0
        }
        simcall <- call("setObserved", simcall)
    }
    res <- if (nsim == 1) 
        eval(simcall)
    else replicate(nsim, eval(simcall), simplify = if (simplify) 
        "array"
        else FALSE)
    if (simplify) {
        dimnames(res)[1:2] <- list(subset, colnames(model$response))
        attr(res, "initial") <- y.start
        attr(res, "stsObserved") <- object$stsObj[subset, ]
        class(res) <- "hhh4sims"
    }
    attr(res, "call") <- cl
    attr(res, "seed") <- RNGstate
    res
}

### Internal auxiliary function, which performs the actual simulation
my.simHHH4lag <- function(ar, ne, end, psi, neW, start, lag.ar, lag.ne, funct_lag, 
                          par_lag, min_lag, max_lag, scale) 
{
    nTime <- nrow(end)
    nUnits <- ncol(end)
    
    ## simulate from Poisson or NegBin model
    rdistr <- if (length(psi) == 0 || 
                  isTRUE(all.equal(psi, 0, check.attributes = FALSE))) {
        rpois
    } else {
        psi.inv <- 1/psi   # since R uses different parametrization
        ## draw 'n' samples from NegBin with mean vector 'mean' (length=nUnits)
        ## and overdispersion psi such that Variance = mean + psi*mean^2
        ## where 'size'=1/psi and length(psi) == 1 or length(mean)
        function(n, mean) rnbinom(n, mu = mean, size = psi.inv)
    }
    
    ## if only endemic component -> simulate independently
    if (all(ar + ne == 0)) {
        return(matrix(rdistr(length(end), end), nTime, nUnits))
    }
    
    ## weighted sum of counts of other (neighbouring) regions
    ## params: y - vector with (lagged) counts of regions #BJ: distributed lags act here.
    ##         W - nUnits x nUnits adjacency/weight matrix (0=no neighbour)
    wSumNE <- if (is.null(neW) || all(neW == 0)) { # includes the case nUnits==1
        function(y, W) numeric(nUnits)
    } else {
        function(y, W) .colSums(W * y, nUnits, nUnits)
    }
    
    # Initialize matrices for means mu_ait and simulated data y_ait
    mu <- y <- matrix(0, nTime, nUnits)
    y <- rbind(start, y)
    nStart <- nrow(y) - nrow(mu)        # usually just 1 for lag=1
    
    ## simulate
    timeDependentWeights <- length(dim(neW)) == 3
    if (!timeDependentWeights) 
        neWt <- neW
    for (t in seq_len(nTime)) {
        if (timeDependentWeights) 
            neWt <- neW[, , t]
        # Calculate lagged local cases
        Ylagged <- hhh4addon:::weightedSumAR(observed = y[nStart +
                                                              t - (max_lag:0), , drop = FALSE], lag = lag.ar,
                                             funct_lag = funct_lag, par_lag = par_lag, min_lag = min_lag,
                                             max_lag = max_lag, sum_up = TRUE)[max_lag + 1, , drop = FALSE]

        if (!is.null(scale)){
            # If the model is age stratified, compute the within-region matrix 
            # of connectivity between age groups
            C <- scale
            # Extract the number of regions 
            regs <- (sub("[.].*", "", colnames(y)))
            nb_reg <- length(unique(regs))
            rownames(C) <- colnames(C) <- colnames(y)
            for(i in seq_len(nb_reg)){
                # Elements linking different regions are set to 0
                reg_i <- regs[i]
                C[regs == reg_i, regs != reg_i] <- 0 
            }
            Ylagged.ar <- Ylagged %*% C
        } else {
            Ylagged.ar <- Ylagged
        }
        # Calculate lagged neighbourhood cases
        if (!is.null(neW)) {
            Ylagged.ne <- hhh4addon:::weightedSumNE(y[nStart + 
                                                          t - (max_lag:0), , drop = FALSE], weights = neWt, 
                                                    lag = lag.ne, funct_lag = funct_lag, par_lag = par_lag, 
                                                    min_lag = min_lag, max_lag = max_lag, sum_up = TRUE)[max_lag + 
                                                                                                             1, ]
        } else {
            Ylagged.ne <- 0
        }
        # Calculate mean of distribution
        mu[t, ] <- ar[t, ] * Ylagged.ar + ne[t, ] * Ylagged.ne + end[t, ]
        # Sample from Poisson/NegBin with that mean
        y[nStart + t, ] <- rdistr(nUnits, mu[t, ])
    }
    
    ## return simulated data without initial counts
    y[-seq_len(nStart), , drop = FALSE]
}

## Internal auxiliary function, which performs the actual simulation, while 
## calculating and updating the cumulative incidence covariates as the simulation 
## progresses
my.simHHH4lag_cumu_inc <- function(psi, neW, start, lag.ar, lag.ne, funct_lag, 
                                   par_lag, min_lag, max_lag, scale, object, 
                                   subset, theta, components, pop, delay, exppreds)
{
    nTime <- length(subset)
    nUnits <- ncol(start)
    
    ## simulate from Poisson or NegBin model
    rdistr <- if (length(psi) == 0 ||
                  isTRUE(all.equal(psi, 0, check.attributes = FALSE))) {
        rpois
    } else {
        psi.inv <- 1/psi   # since R uses different parametrization
        ## draw 'n' samples from NegBin with mean vector 'mean' (length=nUnits)
        ## and overdispersion psi such that Variance = mean + psi*mean^2
        ## where 'size'=1/psi and length(psi) == 1 or length(mean)
        function(n, mean) rnbinom(n, mu = mean, size = psi.inv)
    }

    ## weighted sum of counts of other (neighbouring) regions
    ## params: y - vector with (lagged) counts of regions #BJ: distributed lags act here.
    ##         W - nUnits x nUnits adjacency/weight matrix (0=no neighbour)
    wSumNE <- if (is.null(neW) || all(neW == 0)) { # includes the case nUnits==1
        function(y, W) numeric(nUnits)
    } else {
        function(y, W) .colSums(W * y, nUnits, nUnits)
    }
    ytot <- object$stsObj@observed[seq_len(min(subset) - 1),]
    row_prev <- nrow(ytot)
    
    ## Initialize matrices for means mu_ait and simulated data y_ait
    y <- matrix(0, nTime, nUnits)
    y <- rbind(start, y)
    nStart <- nrow(y) - nTime        # usually just 1 for lag=1
    
    ## Simulate
    timeDependentWeights <- length(dim(neW)) == 3
    
    # Extract autoregressive, neighbourhood, and endemic predictors
    ar_new <- getComp("ar", exppreds, components)
    ne_new <- getComp("ne", exppreds, components)
    end_new <- getComp("end", exppreds, components)
    ar_t <- ar_new[1,]
    ne_t <- ne_new[1,]
    end_t <- end_new[1,]
    
    # Create temporary variables for cumulative incidence covariates
    cumu_new <- object$control$data$inc_new
    cumu_old <- object$control$data$inc_old
    
    # Calculate lag weights
    lag_weights <- funct_lag(par_lag = par_lag, min_lag = min_lag, max_lag = max_lag)
    mat_lag <- t(matrix(lag_weights))
    if (!is.null(scale)){
        # If the model is age stratified, compute the within-region matrix 
        # of connectivity between age groups
        C <- scale
        # Extract the number of regions 
        regs <- (sub("[.].*", "", colnames(y)))
        nb_reg <- length(unique(regs))
        rownames(C) <- colnames(C) <- colnames(y)
        for(i in seq_len(nb_reg)){
            # Elements linking different regions are set to 0
            reg_i <- regs[i]
            C[regs == reg_i, regs != reg_i] <- 0 
        }
    }
    
    if (is.na(object$coefficients["ar.log(1 - inc_new)"])) 
        object$coefficients["ar.log(1 - inc_new)"] <- 0
    if (is.na(object$coefficients["ar.log(1 - inc_old)"])) 
        object$coefficients["ar.log(1 - inc_old)"] <- 0
    if (is.na(object$coefficients["ne.log(1 - inc_new)"])) 
        object$coefficients["ne.log(1 - inc_new)"] <- 0
    if (is.na(object$coefficients["ne.log(1 - inc_old)"])) 
        object$coefficients["ne.log(1 - inc_old)"] <- 0

    if (!timeDependentWeights)
        neWt <- neW
    for (t in seq_len(nTime)) {
        prop_suscep_new <- 1 - object$control$data$inc_new[nStart + t, ]
        prop_suscep_old <- 1 - object$control$data$inc_old[nStart + t, ]
        ar_t <- ar_new[t,] * 
            prop_suscep_new ^ object$coefficients["ar.log(1 - inc_new)"] *
            prop_suscep_old ^ object$coefficients["ar.log(1 - inc_old)"]
        ne_t <- ne_new[t,] * 
            prop_suscep_new ^ object$coefficients["ne.log(1 - inc_new)"] *
            prop_suscep_old ^ object$coefficients["ne.log(1 - inc_old)"]
        
        if (timeDependentWeights)
            neWt <- neW[, , t]
        Ylagged <- mat_lag %*% y[nStart + t - (0:max_lag), , drop = FALSE][-1,]
        if (!is.null(scale)){
            Ylagged.ar <- Ylagged %*% C
        } else {
            Ylagged.ar <- Ylagged
        }
        if (!is.null(neW)) {
            Ylagged.ne <- Ylagged %*% neW
        } else {
            Ylagged.ne <- 0
        }
        # Calculate mean of distribution
        mu <- ar_t * Ylagged.ar + ne_t * Ylagged.ne + end_new[t,]
        ## Sample from Poisson/NegBin with that mean
        y[nStart + t, ] <- rdistr(nUnits, mu)
        # Update cumulative incidence covariates
        if (nStart + t - delay > 0){
            cumu_new[nStart + t, ] <- cumu_new[nStart + t - 1, ] + 
                y[nStart + t - delay, ]/pop[nStart + t, ] - 
                ytot[row_prev + t - 30, ]/pop[nStart + t, ]
            if (row_prev + t - 180 > 0){
                cumu_old[nStart + t, ] <- cumu_old[nStart + t - 1, ] + 
                    ytot[row_prev + t - 30, ]/pop[nStart + t, ] - 
                    ytot[row_prev + t - 180, ]/pop[nStart + t, ]
            } else{
                cumu_old[nStart + t, ] <- cumu_old[nStart + t - 1, ] + 
                    ytot[row_prev + t - 30, ]/pop[nStart + t, ]
            }
        }
        if (t < nTime){
            object$control$data$inc_new[nStart + t + 1, ] <- pmin(cumu_new[nStart + t, ],1) # pmin as cumulative incidence shouldn't go above 1,
            object$control$data$inc_old[nStart + t + 1, ] <- pmin(cumu_old[nStart + t, ],1) # pmin as cumulative incidence shouldn't go above 1,
            # may have an issue with log(0) in model if inc=1? subtract small number from 1 in pmin if so?
        }
        ytot <- rbind(ytot, y[nStart + t,])
    }
    
    ## return simulated data without initial counts
    y[-seq_len(nStart), , drop = FALSE]
}

# Function for extracting autoregressive, neighbourhood and endemic 
# predictors from exppreds object
getComp <- function(comp, exppreds, components) {
    exppred <- exppreds[[comp]]
    if (comp %in% components) 
        exppred
    else `[<-`(exppred, value = 0)
}