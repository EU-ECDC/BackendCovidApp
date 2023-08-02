## Function generating a data frame describing the forecasts in list_pred. 
## The data table contains the prediction quantiles in each region, age group, date,
## for the different scenarios contained in list_pred
format_pred <- function(list_pred, groups, deaths, groups_death){
  ## Initialise output data table
  long_pred_tot <- data.frame(matrix(ncol = 9, nrow = 0))
  colnames(long_pred_tot) <- c("quant", "reg", "age", "t", "t_pred", "age_strat", 
                               "level", "type", "n_cases")
  ## For each element of list_pred (i.e. scenario of forecasts), extract the sim object
  ## compute summary statistics, and add to long_pred_tot
  for(k in seq_along(list_pred)){
    pred_k <- list_pred[[k]]$sim
    array_pred_k <- array(pred_k, dim = dim(pred_k), dimnames = dimnames(pred_k))
    array_pred_k_death <- array(deaths[[k]], dim = dim(deaths[[k]]), dimnames = dimnames(deaths[[k]]))
    ## In each group (nuts level and age-stratification), aggregate the prediction using 
    ## groups[[i]]
    for (i in seq_along(groups)){
      names_agg <- aggregate(t(array_pred_k[, , 1]), list(groups[[i]]), sum)[,1]
      ## Create empty 3-d array for case forecasts
      pred_agg <- array(dim = c(nrow(array_pred_k), length(names_agg), dim(array_pred_k)[3]),
                        dimnames = list(dimnames(array_pred_k)[[1]], names_agg, 
                                        dimnames(array_pred_k)[[3]]))
      ## In each z matrix of pred_agg, aggregate the predictions using groups[[i]]
      for(j in seq_len(dim(array_pred_k)[3])) 
        pred_agg[,,j] <- t(aggregate(t(array_pred_k[, , j]), list(groups[[i]]), sum)[,-1])
      
      ## Create empty 3-d array for death forecasts
      pred_agg_death <- array(dim = c(nrow(array_pred_k_death), length(names_agg), dim(array_pred_k_death)[3]),
                              dimnames = list(dimnames(array_pred_k_death)[[1]], names_agg, 
                                              dimnames(array_pred_k_death)[[3]]))
      
      ## In each z matrix of pred_agg, aggregate the predictions using groups[[i]]
      for(j in seq_len(dim(array_pred_k)[3])){ 
        if(names(groups)[i] %in% c("nuts3", "nuts3_age", "nuts3_tot")){
          pred_agg_death[,,j] <- -1
        } else{
          pred_agg_death[,,j] <- t(aggregate(t(array_pred_k_death[, , j]), 
                                             list(groups_death[[names(groups)[i]]]), sum)[,-1])
        }
      }
      ## Compute the 2.5, 25, 50, 75, and 97.5 case forecast centiles at each date / region 
      quant_pred <- apply(pred_agg, c(1,2), function(X) quantile(X, c(.025, .25, .5, .75, .975)))
      
      quant_pred_death <- quant_pred * NA
      ## Compute the 2.5, 25, 50, 75, and 97.5 death forecast centiles at each date / region 
      quant_pred_death[, rownames(pred_agg_death),] <- apply(pred_agg_death, c(1,2), function(X) 
        quantile(X, c(.025, .25, .5, .75, .975)))
      quant_pred_death[quant_pred_death < 0] <- NA
      
      ## Extract region names and age groups
      reg <- sub("[.].*", "", dimnames(quant_pred)[[3]])
      age <- sub(".*[.]", "", dimnames(quant_pred)[[3]])
      
      ## Create long_pred, data_frame with same columns as long_pred_tot
      # quant: centile 
      # reg: region
      # age; age_group
      # t: date
      # conditions: scenario (taken from the kth name of list_pred)
      # level: nuts / age stratification level
      # type: number of daily cases (will be used in the time-series plot)
      # n_cases: number of cases
      # n_deaths: number of deaths
      long_pred <- cbind.data.frame(
        quant = rep(dimnames(quant_pred)[[1]], prod(dim(quant_pred)[c(2,3)])),
        reg = rep(reg, each = prod(dim(quant_pred)[c(1,2)])),
        age = rep(age, each = prod(dim(quant_pred)[c(1,2)])),
        t = rep(rep(dimnames(quant_pred)[[2]], each = nrow(quant_pred)), dim(quant_pred)[3]),
        conditions = names(list_pred)[k],
        level = names(groups)[i], type = "daily_cases",
        n_cases = c(quant_pred), n_deaths = c(quant_pred_death))
      
      # If non-age-stratified nuts2 and nuts3 level, compute map data
      if (names(groups)[i] == "nuts2_tot"|| names(groups)[i] == "nuts3_tot"){
        ## Compute number of cases in the next two weeks
        pred_week <- colSums(pred_agg[1:14,,])
        pred_death_week <- colSums(pred_agg_death[1:2,,])
        quant_week <- apply(pred_week, 1, function(X) quantile(X, seq(0,1,.1)))
        quant_week_death <- apply(pred_death_week, 1, function(X) quantile(X, seq(0,1,.1)))
        quant_week_death[quant_week_death < 0] <- NA
        ## Create week_pred, data_frame with same columns as long_pred_tot, and merge it to 
        ## long_pred
        week_pred <- cbind.data.frame(
          quant = rep(dimnames(quant_week)[[1]], dim(quant_week)[2]),
          reg = rep(reg, each = prod(dim(quant_week)[1])),
          age = rep(age, each = prod(dim(quant_week)[1])),
          t = NA, conditions = names(list_pred)[k],
          level = names(groups)[i], type = "14_day_cases",
          n_cases = c(quant_week), n_deaths = c(quant_week_death))
        long_pred <- rbind.data.frame(long_pred, week_pred)
      }
      ## Add long_pred to long_pred_tot
      long_pred_tot <- rbind.data.frame(long_pred_tot, long_pred)
    }
  }
  
  return(long_pred_tot)
}
