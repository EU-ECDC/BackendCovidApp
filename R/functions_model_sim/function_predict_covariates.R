### function_predict_covariates.R contains functions for predicting covariates for case forecasts:
## - predict_covariates: Wrapper function for handling different types of covariate (continous, binary, etc.) and different options for predicting covariates
## - predict_cov: Predict covariate between two time points with a specified prediction type (current/average/actual)

predict_covariates <- function(model, t_start, t_end, country, actual = FALSE, default = TRUE, type = NULL, delay = 30){
    # Extract covariate data from model object
    covariates <- model$control$data
    # Drop time variable for the time being
    covariates$t <- NULL
    
    # Get last observed day
    t_max <- nrow(model$stsObj@observed)
    
    # Get names of covariates
    nms <- names(covariates)
    
    # Predict values of covariates
    covariates_pred <- covariates
    if (actual){
        # Use actual values if actual = TRUE and prediction window is within observed period
        stopifnot("t_end must be within observed time period to use type = 'actual'." = (t_end <= t_max))
        if (!missing(default) & default){
            warning("Ignoring 'default' as actual = TRUE.")    
        }
        for (i in seq_along(covariates_pred)){
            covariates_pred[[i]] <- predict_cov(covariates[[i]], t_start, t_end, type = "actual")     
        }
    } else {
        # Otherwise use:
        if (default){
            # default types if default = TRUE
            if (!is.null(type)){
                warning("Ignoring 'type' as default = TRUE.")
            }
            for (i in seq_along(covariates_pred)){
                if (nms[i] %in% c("pop","pop_age","cov2","cov3","cov2_1","cov_tot",
                                  "test","inc_WT","inc_delta","inc_tot", "test_prop", "test_age")){
                    # populations, vaccine coverages, previous incidence, test proportions: current values (type = "current")
                    covariates_pred[[i]] <- predict_cov(covariates[[i]], t_start, t_end, type = "current")
                } else if (nms[i] %in% c("europe","SI","MI")){
                    # number of cases in the rest of Europe over the past month, NPI stringency index, mobility index: mean of values in past 'delay' days (type = "average")    
                    covariates_pred[[i]] <- predict_cov(covariates[[i]], t_start, t_end, type = "average", delay = delay)
                } else if (nms[i] %in% c("sat", "sun", "monday", "tues", "wed", "thu", "fri")) {
                    # day-of-week effect
                    n_weeks <- (t_end - t_start)%/%7 + 1
                    covariates_pred[[i]] <- covariates_pred[[i]][((t_start%%7) : (t_start%%7 + 28)) + 1,]
                    covariates_pred[[i]] <- covariates_pred[[i]][seq_len(t_end - t_start),]
                    days <- as.character(as.Date(epoch(model$stsObj), origin = "1970-01-01")[(t_start + 1):t_end])
                    if(nms[i] == "sun" & country != "IT"){
                      covariates_pred[[i]][is.element(days, bank_holiday_country(country)), ] <- 1
                    } else
                      covariates_pred[[i]][is.element(days, bank_holiday_country(country)), ] <- 0
                } else { 
                    # if covariate isn't in lists above, use current value for predicted values
                    covariates_pred[[i]] <- predict_cov(covariates[[i]], t_start, t_end, type = "current")
                }
            }
        } else {
            # given type if 'type' is specified
            stopifnot("type must be specified if default = TRUE" = !is.null(type))
            for (i in seq_along(covariates_pred)){
                if (names(covariates_pred)[i] %in% names(type)){
                    covariates_pred[[i]] <- predict_cov(covariates[[i]], t_start, t_end, type = type[names(type) == names(covariates_pred)[i]], delay = delay)    
                } else { # if type of projection isn't specified, use current value
                   covariates_pred[[i]] <- predict_cov(covariates[[i]], t_start, t_end, type = "current")
                }
            } 
        }        
    }
    
    covariates_pred <- c(list(t = ((t_start+1):t_end) - 1), covariates_pred)
        
    return(covariates_pred)
}

predict_cov <- function(x, t_start, t_end, type = "current", delay = 30){
    ntps <- t_end - t_start
    if (type == "current"){
        x_pred <- matrix(x[t_start,], nrow = ntps, ncol = ncol(x), byrow = T)    
    } else if (type == "average"){
        x_pred <- matrix(colMeans(x[t_start-((delay-1):0),]), nrow = ntps, ncol = ncol(x), byrow = T)
    } else if (type == "actual"){
        x_pred <- x[(t_start+1):t_end,]
    } else {
        stop("type must be one of 'current', 'average' or 'actual'.")
    }
}
