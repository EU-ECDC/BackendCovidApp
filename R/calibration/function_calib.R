compute_calib <- function(country, range_dates, download, data_file, total, dates, nsim = 1000, empty = F){
  #### Import and process data
  ## Import all data, make sts object, and generate covariates
  list_import <- import_all_files(country = country, total = total,
                                  range = range_dates, delay_cumu = 0,
                                  download_case = download, file_case = data_file, last_date = TRUE)
  ## Extract elements of list_import
  sts_obj_nei <- list_import$sts_obj_nei
  covariates <- list_import$covariates
  C <- list_import$C
  ## Create covariate containing age groups
  GROUPS <- names(covariates)[substr(names(covariates), 1, 3) == "age"]
  
  #### Run model 
  ### Parameter of the serial interval model: 80% of the serial interval comes from direct tranmission
  par_dist <- -log(1/.8 -1)
  ### Define the equation of the model:
  if (total){
    ## Equation if the model is not age-stratified
    ar_terms <- NULL
    end_terms <- c("1", 
                   "rural", "int_rur", "int_urb"
                   , "log(pop)", "log(europe/1e8)"
    )
    ne_terms <- c("1", "log(pop)"
                  , "log(test_prop)", "tues", "wed", "thu", "fri", "sat", "sun"
                  , "rural", "int_rur", "int_urb", "log(1 - cov_tot)"
                  , "log(1 - inc_old)", "log(1 - inc_new)"
                  , "delta", "omicron"
    )
    if(empty) end_terms <- ne_terms <- c("1")
    ## Remove covariates that are always equal to 0
    for (i in seq_along(covariates)){
      if (sum(covariates[[i]]) == 0){
        nm <- names(covariates)[i]
        ar_terms <- ar_terms[!grepl(nm,ar_terms)]
        end_terms <- end_terms[!grepl(nm,end_terms)]
        ne_terms <- ne_terms[!grepl(nm,ne_terms)]
      }
    }
    ## If all values of cov_tot are 0 as subnational vaccination data is not 
    ## available, add log(1 - cov_tot_nat_tot) (overall national coverage) covariate
    if (sum(covariates$cov_tot) == 0){
      ne_terms <- c(ne_terms, "log(1 - cov_nat_tot)")
    }
  } else {
    ## Equation if the model is age-stratified
    ar_terms <- NULL
    end_terms <- c("1", "log(europe/1e8)"
                   , "log(pop * pop_age)", "rural", "int_rur", "int_urb"
    )
    ne_terms <- c("1", GROUPS[-1] 
                  , "log(test_age)", "log(test_prop)", "log(pop_age)", "log(pop)"
                  , "tues", "wed", "thu", "fri", "sat", "sun"
                  , "rural", "int_rur", "int_urb"
                  , "log(1 - cov_tot)"
                  , "log(1 - inc_old)"
                  , "log((1 - inc_new))"
                  , "delta", "omicron"
    )
    if(empty) end_terms <- ne_terms <- c("1")
    
    ## If all values of test_age are equal to 1 as age-stratified national testing
    ## data is not available, then drop log(test_age) covariate
    if (all(covariates$test_age == 1)){
      ne_terms <- ne_terms[ne_terms != "log(test_age)"]
    }
    ## If all values of cov_tot are 0 as subnational vaccination data is not 
    ## available, add log(1 - cov_nat_tot) (national age-stratified coverage) covariate
    if (sum(covariates$cov_tot) == 0){
      ne_terms <- c(ne_terms, "log(1 - cov_nat_tot)")
    }
    
    ## Remove covariates that are always equal to 0
    for (i in seq_along(covariates)){
      if (sum(covariates[[i]]) == 0){
        nm <- names(covariates)[i]
        ar_terms <- ar_terms[!grepl(nm,ar_terms)]
        end_terms <- end_terms[!grepl(nm,end_terms)]
        ne_terms <- ne_terms[!grepl(nm,ne_terms)]
      }
    }
  }
  
  ## Add seasonality to ne components
  if(empty){
    equations <- 
      list(ar = ar_terms,
           ne = reformulate(ne_terms, intercept = TRUE),
           end = reformulate(end_terms, intercept = TRUE))
  } else { 
    equations <- 
      list(ar = ar_terms,
           ne = addSeason2formula(reformulate(ne_terms, intercept = TRUE), period = 365),
           end = addSeason2formula(reformulate(end_terms, intercept = TRUE))
    )
  }
  
  ## Run the model
  model_fit <-
    hhh4_cov(sts_object = sts_obj_nei, contact = C/max(C), equations = equations,
             covariates = covariates, lag_dist = lag_daily, parameter_lag = par_dist,
             weight_fun = W_powerlaw(maxlag = 5, log = T, normalize = T), empty = empty)
  
  #### Generate and analyse predictions
  ### Initialise list of predictions
  pred_list <- list()
  q_pred_list <- list()
  q_pred_reg_list <- list()
  
  ### 28-day forecasts since March 2022
  ## Initialise vector of prediction dates, prediction and observation matrices
  t_pred <- which(is.element(as.Date(epoch(sts_obj_nei), origin = "1970-01-01"), dates))
  pred_week <- low_pred_week <- up_pred_week <- matrix(, nrow = length(t_pred), ncol = 4)
  obs_week <- numeric()
  
  obs <- sts_obj_nei@observed
  obs_week_reg <- matrix(nrow = length(t_pred), ncol = ncol(obs))
  colnames(obs_week_reg) <- colnames(obs)
  rownames(obs) <- as.character(as.Date(epoch(sts_obj_nei), origin = "1970-01-01"))
  
  for(i in seq_along(t_pred)){
    ## Prediction date
    t_i <- t_pred[i]
    model_pred <- model_fit
    ## Set random number seed for simulations
    set.seed(1)
    ## Set extreme values of standard deviations to 0
    if(any(model_pred$se > 10)) {
      model_pred$coefficients[-length(model_pred$coefficients)] <- 0
      model_pred$coefficients[length(model_pred$coefficients)] <- 2
    }
    ## Generate predictions
    pred_list[[i]] <- nStepAhead(model = model_pred, n = 28, nsim = nsim/10, t_start = t_i,
                                 delay = 0, t_end = t_i, total = total, nparam = 10, country = country)
    ## Aggregate: compute country-wide weekly predictions
    pred_i <- apply(pred_list[[i]]$sim, 3, rowSums)
    rownames(pred_i) <- ceiling((as.numeric(rownames(pred_i)) - as.numeric(rownames(pred_i))[1] + 1)/7)
    tmp <- rowsum(pred_i, group = rownames(pred_i))
    dt_tmp <- melt(as.data.table(tmp, keep.rownames = "horizon"), id.vars = "horizon")
    q_pred_i <- dt_tmp[, as.list(
      quantile(value, 
               probs = c(0.01, 0.025, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.45, 0.5, 
                         0.55, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.975, 0.99))), 
      by = .(horizon)]
    pred_week[i,] <- t(q_pred_i[, `50%`])
    low_pred_week[i,] <- t(q_pred_i[, `2.5%`])
    up_pred_week[i,] <- t(q_pred_i[, `97.5%`])
    # Calculate weekly number of cases in data
    if(t_i + 7 <= nrow(obs)) obs_week[i] <- sum(obs[t_i + seq_len(7), ])
    if(t_i + 7 <= nrow(obs)) obs_week_reg[i,] <- colSums(obs[t_i + seq_len(7), ])
    q_pred_list[[i]] <- q_pred_i
    q_pred_reg_list_j <- list()
    for(j in seq_len(ncol(pred_list[[i]]$sim))){
      pred_j <- apply(pred_list[[i]]$sim[,j,], 3, rowSums)
      rownames(pred_j) <- ceiling((as.numeric(rownames(pred_j)) - as.numeric(rownames(pred_j))[1] + 1)/7)
      tmp <- rowsum(pred_j, group = rownames(pred_j))
      dt_tmp <- melt(as.data.table(tmp, keep.rownames = "horizon"), id.vars = "horizon")
      dt_tmp$reg <- colnames(pred_list[[i]]$sim)[j]
      # q_pred_reg_list_j[[j]] <- q_pred_j
      q_pred_reg_list_j[[j]] <- dt_tmp
    }
    dt_i <- rbindlist(q_pred_reg_list_j)
    if(!all(sub("[.].*", "", dt_i$reg) == dt_i$reg)){
      dt_tot <- dt_i[, lapply(.SD, sum), by = .(horizon, variable, reg = paste0(sub("[.].*", "", reg), ".tot")), 
                     .SDcols = "value"]
      dt_i <- as.data.table(rbind.data.frame(dt_i, dt_tot))
      dt_age <- dt_i[, lapply(.SD, sum), by = .(horizon, variable, reg = paste0("tot.", sub(".*[.]", "", reg))), 
                     .SDcols = "value"]
      dt_i <- as.data.table(rbind.data.frame(dt_i, dt_age))
      dt_i[, age := sub(".*[.]", "", reg)]
      dt_i[age == "0-9" | age == "10-19", age := "0-20"]
      dt_i[age == "20-29" | age == "30-39" | age == "40-49" | age == "50-59", 
           age := "20-60"]
      dt_i[age == "60-69" | age == "70-79", age := "60-80"]
      dt_i[, reg := sub("[.].*", "", reg)]
      dt_i <- dt_i[, lapply(.SD, sum), by = .(horizon, variable, reg, age), .SDcols = "value"]
      dt_i[, reg := paste(reg, age, sep = ".")]
      dt_i[, age := NULL]
    }
    q_pred_reg_i <- dt_i[, as.list(
      quantile(value,
               probs = c(0.01, 0.025, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.45, 0.5,
                         0.55, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.975, 0.99))),
      by = .(horizon, reg)]
    
    q_pred_reg_list[[i]] <- q_pred_reg_i
  }
  names(obs_week) <- rownames(obs_week_reg) <- 
    as.character(as.Date(epoch(sts_obj_nei), origin = "1970-01-01"))[t_pred][seq_along(obs_week)]
  rownames(pred_week) <- as.character(as.Date(epoch(sts_obj_nei), origin = "1970-01-01"))[t_pred]
  q_pred <- rbindlist(q_pred_list, idcol = "i")
  q_pred[, date := dates[i]]
  q_pred[, i := NULL]
  q_pred <- melt(q_pred, id.vars = c("date", "horizon"), variable.name = "quantile", value.name = "prediction")
  q_pred[, quantile := as.numeric(sub("%", "", quantile))/100]
  
  q_pred_reg <- rbindlist(q_pred_reg_list, idcol = "i")
  q_pred_reg[, date := dates[i]]
  q_pred_reg[, i := NULL]
  ## Generate calib, containing the number of weekly cases (predicted and in data)
  calib <- list(obs = obs_week, pred_week = pred_week, low = low_pred_week, up = up_pred_week, q_pred = q_pred,
                q_pred_reg = q_pred_reg, obs_reg = obs_week_reg)
  names(pred_list) <- t_pred
  pred <- pred_list[[1]]
  
  ## Compute prediction scores
  # Produces a list of lists (scores_list), each element of which is a list of 
  # the scores for all region(-age) groups for the 1-, 2-, 3-, and 4-week-ahead
  # forecasts at each of the prediction dates (dates)
  scores_list <- list()
  for (i in seq_along(pred_list)){
    sim_i <- pred_list[[i]]$sim
    class(sim_i) <- "array"
    t_i <- as.Date(epoch(sts_obj_nei), origin = "1970-01-01")[as.numeric(names(pred_list))[i]]
    t_pred_i <- as.numeric(names(pred_list[i]))
    if(t_pred_i + 1 < nrow(sts_obj_nei@observed)){
      # Get daily number of cases for the next 28 days following the prediction date
      obs_i <- sts_obj_nei@observed[seq(t_pred_i + 1, min(t_pred_i + 28, nrow(sts_obj_nei@observed))),]
      sim_i_dt <- as.data.table(sim_i)
      names(sim_i_dt) <- c("target_end_date", "location", "sample", "prediction")
      sim_i_dt[, forecast_date := as.Date(t_i)]
      sim_i_dt[, target_end_date := as.Date(epoch(sts_obj_nei), origin = "1970-01-01")[as.integer(target_end_date)]]
      sim_i_dt[, horizon := as.numeric(target_end_date - as.Date(t_i))/7] # N.B. horizon is in weeks here for consistency with death forecasts
      tmp <- as.data.table(obs_i)
      tmp[, target_end_date := as.Date(epoch(sts_obj_nei), origin = "1970-01-01")[seq(t_pred_i + 1, min(t_pred_i + 28, nrow(sts_obj_nei@observed)))]]
      obs_i_dt <- melt(tmp, id.vars = "target_end_date", variable.name = "location", value.name = "true_value")
      nsa_i <- merge(obs_i_dt, sim_i_dt, by = c("target_end_date","location"))
      
      # scores_list[[i]] <- score(nsa_i[target_end_date %in% seq.Date(t_i + 7, min(t_i + 28, max(target_end_date)), by = 7)], metrics = c("crps","dss","se_mean"))
      scores_list[[i]] <- score(nsa_i, metrics = c("crps","dss","se_mean"))
    }
  }
  # Combine the scores for all region(-age) groups, prediction dates and forecast horizons into one data table
  scores <- rbindlist(scores_list)
  scores[, country := country]
  scores[, `:=`(reg_nb = sub("\\..*", "", location),
                age = ifelse(country %in% c("CZ","FR"), sub(".*\\.", "", location), "Total"))]
  
  # Round up the forecast horizon, so that scores are averaged over daily forecasts in the same forecast week
  scores[, horizon := ceiling(horizon)]
  # Calculate the mean scores over the different prediction dates for each region(-age) group and forecast horizon
  mean_scores <- scores[, lapply(.SD, mean), .SDcols = c("crps", "dss", "se_mean"), 
                        by = .(country, location, reg_nb, age, horizon)]
  
  ## Compute pit histogram
  pxm1 <- px <- 
    array(NA, dim = c(dim(pred_list[[1]]$sim)[c(1, 2)], length(pred_list)))
  
  for(i in seq_along(pred_list)){
    # Extract predictions at i
    pred_i <- pred_list[[i]]$sim
    # Compute prediction date 
    t_pred_i <- as.numeric(names(pred_list)[i])
    if(t_pred_i + 1 < nrow(sts_obj_nei@observed)){
      # Get daily number of cases for the next 28 days following the prediction date
      obs_i <- sts_obj_nei@observed[seq(t_pred_i + 1, min(t_pred_i + 28, nrow(sts_obj_nei@observed))),]
      pred_i <- pred_i[seq_len(nrow(obs_i)), , ]
      
      # Create array with daily number of cases in the data per age / region
      x_i <- array(obs_i, dim = dim(pred_i))
      xm1_i <- array(obs_i - 1, dim = dim(pred_i))
      
      # Compute proportion of simulations where the number of predicted case is below x_i
      px[seq_len(nrow(obs_i)), , i] <- t(apply(pred_i <= x_i, 1, function(X) return(rowSums(X)/ncol(X))))
      # Compute proportion of simulations where the number of predicted case is below x_i - 1
      pxm1[seq_len(nrow(obs_i)), , i] <- t(apply(pred_i <= xm1_i, 1, function(X) return(rowSums(X)/ncol(X))))
      # If x is 0, set pxm1 to 0
      pxm1[seq_len(nrow(obs_i)), , i][obs_i == 0] <- 0
    }
  }
  colnames(px) <- colnames(pxm1) <- colnames(sts_obj_nei)
  pit <- list(px = px, pxm1 = pxm1)
  
  ## Compute comparison between median estimate and case data
  list_compar_cases <- list()
  for(j in 1:4){
    obs_vec <- numeric()
    med_pred_vec <- numeric()
    max_axes <- quantile(obs[-c(1:t_pred[1]),], .99)
    for(i in seq_along(t_pred)){
      ## Prediction date
      t_i <- t_pred[i]
      pred_i <- pred_list[[i]]$sim
      med_pred <- pred_i %>% apply(1, function(X) apply(X, 1, median)) %>% t
      
      if(max(t_i + 1:7 + (j -1) * 7) <= nrow(obs)){
        obs_vec <- c(obs_vec, obs[t_i + 1:7 + (j -1) * 7,])
        med_pred_vec <- c(med_pred_vec, med_pred[1:7 + (j -1) * 7,])
      }
    }
    list_compar_cases[[j]] <- rbind(obs_vec, med_pred_vec)
  }
  
  #### Generate predictions of numbers of deaths
  ## Compute cfr
  dt_cfr <- calculate_cfr(country = country, download = T, data_file = data_file, 
                          death_data_file = "", total = total, level = "nuts2")
  min_date_omicron <- min(as.Date(epoch(sts_obj_nei), origin = "1970-01-01")[covariates$omicron[,1] == 1])
  ## Prediction of number of deaths
  death_pred <- forecast_death(predictions = list(pred_list), 
                               dt_cfr = dt_cfr[date > min_date_omicron & !is.na(CFR),], 
                               cases = sts_obj_nei@observed, country = country,
                               dates = as.Date(epoch(sts_obj_nei), origin = "1970-01-01"),
                               total = total)[[1]]
  ## Import data on number of deaths
  data_death_mat <- data_death(country = country, download = T, death_data_file = "", 
                               total = total, pred_date = max(as.Date(names(obs_week)) + 7))
  data_death_mat <- data_death_mat[, colnames(death_pred[[1]])]
  # Compute national number of deaths
  obs_death <- rowSums(data_death_mat)[names(obs_week)]
  # Compute median predicted number of weekly deaths, and 95% prediction interval
  death_pred_nat <- lapply(death_pred, function(x) apply(x, c(1,3), sum))
  death_pred_reg_list <- list()
  for(j in seq_len(ncol(death_pred[[1]]))){
    death_pred_j <- lapply(death_pred, function(x) x[,j,])
    tmp_reg <- lapply(death_pred_j, function(x) {
      melt(as.data.table(x, keep.rownames = "t"), id.vars = "t")
      })
    death_pred_reg_list[[j]] <- rbindlist(tmp_reg, idcol = "i")
    names(death_pred_reg_list)[j] <- colnames(death_pred[[1]])[j]
  }
  tmp <- lapply(death_pred_nat, function(x) melt(as.data.table(x, keep.rownames = "t"), id.vars = "t"))
  dt_tmp <- rbindlist(tmp, idcol = "i")
  dt_tmp[, `:=`(date = as.Date(epoch(sts_obj_nei), origin = "1970-01-01")[as.integer(i)],
                horizon = (as.integer(t) - as.integer(i))/7,
                sample = as.integer(sub("V", "", variable)))]
  dt_tmp[, (c("variable", "i", "t")) := NULL]
  q_death_pred <- dt_tmp[, as.list(
    quantile(value, probs = c(0.01, 0.025, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.45, 0.5, 
                              0.55, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.975, 0.99))), 
    by = .(date, horizon)
  ]
  q_death_pred <- melt(q_death_pred, id.vars = c("date", "horizon"), variable.name = "quantile", value.name = "prediction")
  mat_death_pred <- dcast(q_death_pred[quantile == "50%",], date ~ horizon, value.var = "prediction")
  mat_death_pred <- as.matrix(mat_death_pred[,!"date"], rownames.value = as.character(mat_death_pred$date))
  low_death_pred <- dcast(q_death_pred[quantile == "2.5%"], date ~ horizon, value.var = "prediction")
  low_death_pred <- as.matrix(low_death_pred[,!"date"], rownames.value = as.character(low_death_pred$date))
  up_death_pred <- dcast(q_death_pred[quantile == "97.5%"], date ~ horizon, value.var = "prediction")
  up_death_pred <- as.matrix(up_death_pred[,!"date"], rownames.value = as.character(up_death_pred$date))
  q_death_pred[, quantile := as.numeric(sub("%", "", quantile))/100]
  dt_death_reg <- rbindlist(death_pred_reg_list, idcol = "reg")
  if(!all(sub("[.].*", "", dt_death_reg$reg) == dt_death_reg$reg)){
    dt_death_tot <- dt_death_reg[, lapply(.SD, sum), by = .(i, t, variable, reg = paste0(sub("[.].*", "", reg), ".tot")), 
                               .SDcols = "value"]
    dt_death_reg <- as.data.table(rbind.data.frame(dt_death_reg, dt_death_tot))
    dt_death_age <- dt_death_reg[, lapply(.SD, sum), by = .(i, t, variable, reg = paste0("tot.", sub(".*[.]", "", reg))), 
                               .SDcols = "value"]
    dt_death_reg <- as.data.table(rbind.data.frame(dt_death_reg, dt_death_age))
  }
  dt_death_reg[, `:=`(date = as.Date(epoch(sts_obj_nei), origin = "1970-01-01")[as.integer(i)],
                horizon = (as.integer(t) - as.integer(i))/7,
                sample = as.integer(sub("V", "", variable)))]
  dt_death_reg[, (c("variable", "i", "t")) := NULL]
  
  if(!all(sub("[.].*", "", dt_death_reg$reg) == dt_death_reg$reg)){
    dt_death_reg[, age := sub(".*[.]", "", reg)]
    dt_death_reg[age == "0-9" | age == "10-19", age := "0-20"]
    dt_death_reg[age == "20-29" | age == "30-39" | age == "40-49" | age == "50-59", 
                 age := "20-60"]
    dt_death_reg[age == "60-69" | age == "70-79", age := "60-80"]
    dt_death_reg[, reg := sub("[.].*", "", reg)]
    dt_death_reg <- dt_death_reg[, lapply(.SD, sum), by = .(horizon, sample, reg, age, date), 
                                 .SDcols = "value"]
    dt_death_reg[, reg := paste(reg, age, sep = ".")]
    dt_death_reg[, age := NULL]
  }
  
  
  q_death_reg <- dt_death_reg[, as.list(
    quantile(value, probs = c(0.01, 0.025, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.45, 0.5, 
                              0.55, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.975, 0.99))), 
    by = .(date, horizon, reg)
  ]
  
  # Put death predictions and data in the list calib_death 
  calib_death <- list(obs = obs_death, pred_week = mat_death_pred, low = low_death_pred, 
                      up = up_death_pred, q_pred = q_death_pred, q_pred_reg = q_death_reg,
                      obs_reg = data_death_mat[is.element(rownames(data_death_mat), 
                                                          names(obs_week)),])
  
  ## Compute death prediction scores
  # Produces a list of data tables (scores_death_list), each of which contains 
  # the scores for all region(-age) groups for the 1-, 2-, 3-, and 4-week-ahead
  # death forecasts at each of the prediction dates (dates)
  scores_death_list <- list()
  for (i in seq_along(death_pred)){
    sim_i <- death_pred[[i]]
    t_i <- as.Date(epoch(sts_obj_nei), origin = "1970-01-01")[as.numeric(names(death_pred))[i]]
    t_pred_i <- which(rownames(data_death_mat) == as.character(t_i))
    if(length(t_pred_i) == 0) t_pred_i <- -1
    if(t_pred_i != -1 & (t_pred_i + 1) <= nrow(data_death_mat)){
      # Extract observed number of death in the next 4 weeks
      obs_i <- data_death_mat[seq(t_pred_i + 1, min(t_pred_i + 4, nrow(data_death_mat))), , drop = F]
      if(!is.matrix(obs_i)){
        obs_i <- t(as.matrix(obs_i))
        sim_i <- array(sim_i[1, , ], dim = c(1, dim(sim_i)[c(2,3)]), 
                       dimnames = list(rownames(sim_i)[1], dimnames(sim_i)[[2]], 
                                       dimnames(sim_i)[[3]]))
      } else {
        sim_i <- sim_i[seq_len(nrow(obs_i)), , , drop = F]
      }
      
      ## Aggregate predictions per week / nuts2 regions (since most death dataset are at nuts2 level)
      if(ncol(sim_i) != ncol(data_death_mat)){
        ## Extract nuts2 regions from colnames(pred_list[[1]][[1]])
        if(total){ 
          colnames(sim_i) <- paste0(get_nuts2_reg(x = colnames(pred_list[[1]][[1]]), country))
        } else
          colnames(sim_i) <- paste0(get_nuts2_reg(x = colnames(pred_list[[1]][[1]]), country), 
                                    sub(".*[.]", ".", colnames(pred_list[[1]][[1]])))
        # Initialise aggregated predictions
        sim_agg <- array(, dim = c(nrow(sim_i), length(unique(colnames(sim_i))), 
                                   dim(sim_i)[3]))
        colnames(sim_agg) <- unique(colnames(sim_i))
        
        # Compute aggregated predictions
        if(nrow(sim_i) == 1){
          for(j in seq_len(dim(sim_i)[3])) 
            sim_agg[, , j] <- sapply(by(sim_i[, , j], names(sim_i[, , j]), sum),identity)[colnames(sim_agg)]
        } else{
          for(j in seq_len(dim(sim_i)[3])) 
            sim_agg[, , j] <- sapply(by(t(sim_i[, , j]),colnames(sim_i[, , j]),colSums),identity)[, colnames(sim_agg)]
        }
      } else {
        sim_agg <- sim_i
      }
      sim_agg_dt <- as.data.table(sim_agg)
      names(sim_agg_dt) <- c("target_end_date", "location", "sample", "prediction")
      sim_agg_dt[, forecast_date := as.Date(t_i)]
      sim_agg_dt[, target_end_date := as.Date(epoch(sts_obj_nei), origin = "1970-01-01")[as.integer(target_end_date)]]
      sim_agg_dt[, horizon := as.numeric(target_end_date - as.Date(t_i))/7]
      obs_i_dt <- melt(as.data.table(obs_i, keep.rownames = "target_end_date"), id.vars = "target_end_date", variable.name = "location", value.name = "true_value")
      obs_i_dt[, target_end_date := as.Date(target_end_date)]
      nsa_i <- merge(obs_i_dt, sim_agg_dt, by = c("target_end_date","location"))
      
      scores_death_list[[i]] <- score(nsa_i, metrics = c("crps","dss","se_mean"))
    }
  }
  # Combine the scores for all region(-age) groups, prediction dates and forecast horizons into one data table
  scores_death <- rbindlist(scores_death_list)
  scores_death[, country := country]
  scores_death[, `:=`(reg_nb = sub("\\..*", "", location),
                      age = ifelse(country %in% c("CZ","FR"), sub(".*\\.", "", location), "Total"))]
  
  # Calculate the mean scores over the different prediction dates for each region(-age) group and forecast horizon
  mean_scores_death <- scores_death[, lapply(.SD, mean), .SDcols = c("crps", "dss", "se_mean"), 
                                    by = .(country, location, reg_nb, age, horizon)]
  
  ## Compute pit histogram death forecasts
  pxm1_death <- px_death <- 
    array(NA, dim = c(dim(death_pred[[1]])[1], ncol(data_death_mat), length(death_pred)))
  for(i in seq_along(death_pred)){
    ## Extract death predictions and prediction date
    pred_i <- death_pred[[i]]
    t_i <- as.Date(epoch(sts_obj_nei), origin = "1970-01-01")[as.numeric(names(death_pred))[i]]
    t_pred_i <- which(rownames(data_death_mat) == as.character(t_i))
    if(length(t_pred_i) == 0) t_pred_i <- -1
    if(t_pred_i != -1 & (t_pred_i + 1) <= nrow(data_death_mat)){
      # Extract observed number of death in the next 4 weeks
      obs_i <- data_death_mat[seq(t_pred_i + 1, min(t_pred_i + 4, nrow(data_death_mat))), ]
      if(!is.matrix(obs_i)){
        obs_i <- t(as.matrix(obs_i))
        pred_i <- array(pred_i[1, , ], dim = c(1, dim(pred_i)[c(2,3)]), 
                        dimnames = list(rownames(pred_i)[1], dimnames(pred_i)[[2]], 
                                        dimnames(pred_i)[[3]]))
      } else
        pred_i <- pred_i[seq_len(nrow(obs_i)), , ]
      
      ## Aggregate predictions per week / nuts2 regions (since most death dataset are at nuts2 level)
      if(ncol(pred_i) != ncol(data_death_mat)){
        ## Extract nuts2 regions from colnames(pred_list[[1]][[1]])
        if(total){ 
          colnames(pred_i) <- paste0(get_nuts2_reg(x = colnames(pred_list[[1]][[1]]), country))
        } else
          colnames(pred_i) <- paste0(get_nuts2_reg(x = colnames(pred_list[[1]][[1]]), country), 
                                     sub(".*[.]", ".", colnames(pred_list[[1]][[1]])))
        # Initialise aggregated predictions
        pred_agg <- array(, dim = c(nrow(pred_i), length(unique(colnames(pred_i))), 
                                    dim(pred_i)[3]))
        colnames(pred_agg) <- unique(colnames(pred_i))
        
        # Compute aggregated predictions
        if(nrow(pred_i) == 1){
          for(j in seq_len(dim(pred_i)[3])) 
            pred_agg[, , j] <- sapply(by(pred_i[, , j], names(pred_i[, , j]), sum),identity)[colnames(pred_agg)]
        } else{
          for(j in seq_len(dim(pred_i)[3])) 
            pred_agg[, , j] <- sapply(by(t(pred_i[, , j]),colnames(pred_i[, , j]),colSums),identity)[, colnames(pred_agg)]
        }
      } else pred_agg <- pred_i
      # Create array with weekly number of deaths in the data per age / region
      x_i <- array(obs_i, dim = dim(pred_agg))
      xm1_i <- array(obs_i - 1, dim = dim(pred_agg))
      
      # Compute proportion of simulations where the number of predicted death is below x_i
      px_death[seq_len(nrow(obs_i)), , i] <- t(apply(pred_agg <= x_i, 1, function(X) return(rowSums(X)/ncol(X))))
      # Compute proportion of simulations where the number of predicted death is below xm1_i
      pxm1_death[seq_len(nrow(obs_i)), , i] <- t(apply(pred_agg <= xm1_i, 1, function(X) return(rowSums(X)/ncol(X))))
      pxm1_death[seq_len(nrow(obs_i)),,i][obs_i == 0] <- 0
    }
  }
  colnames(px_death) <- colnames(pxm1_death) <- colnames(data_death_mat)
  pit_death <- list(px = px_death, pxm1 = pxm1_death)
  
  ## Compute comparison between median estimate and death data 
  list_compar_deaths <- list()
  for(j in 1:4){
    obs_vec <- numeric()
    med_pred_vec <- numeric()
    for(i in seq_along(t_pred)){
      ## Prediction date
      t_i <- as.Date(epoch(sts_obj_nei), origin = "1970-01-01")[as.numeric(names(death_pred))[i]]
      t_pred_i <- which(rownames(data_death_mat) == as.character(t_i))
      pred_i <- death_pred[[i]]
      med_pred <- pred_i %>% apply(1, function(X) apply(X, 1, median)) %>% t
      if(length(t_pred_i) == 0) t_pred_i <- Inf
      if((t_pred_i + j) <= nrow(data_death_mat)){
        obs_vec <- c(obs_vec, data_death_mat[t_pred_i + j,])
        med_pred_vec <- c(med_pred_vec, med_pred[j,])
      }
    }
    list_compar_deaths[[j]] <- rbind(obs_vec, med_pred_vec)
  }
  return(list(calib = calib, calib_death = calib_death, scores = mean_scores, 
              scores_death = mean_scores_death, pit = pit, pit_death = pit_death,
              compar_cases = list_compar_cases, compar_deaths = list_compar_deaths))
}

run_calib <- function(country, all_total, empty, run, pred_date, min_date = NA){
  range_dates <- c("2020-09-01", pred_date)
  
  download <- T
  data_file <- c("")
  
  # Output directory
  dir_out <- if(empty){
    if(all_total){
      "Output/total_empty/"
    } else {
      "Output/age_structured_empty/"
    }
  } else {
    if(all_total){
      "Output/total/"
    } else {
      "Output/age_structured/"
    }
  }
  dir.create(dir_out, recursive = T)
  
  #### Import ensemble forecasts ####
  
  ## Set forecast dates
  dates <- seq(as.Date("2022-04-25"), as.Date(pred_date), 7)
  dates <- dates[dates > (as.Date(pred_date) - 180)]
  # dates <- dates[1:5]
  dates_pred <- dates - 2
  dates <- dates[dates_pred <= max(range_dates)]
  dates_pred <- dates_pred[dates_pred <= max(range_dates)]
  
  if (dates_pred[length(dates_pred)] - dates_pred[1] < 28){
    stop("Final prediction date must be at least 4 weeks after first prediction date")
  }
  
  ## Set quantiles of interest
  quants <- c(0.01, 0.025, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.45, 0.5,
              0.55, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.975, 0.99)
  
  ## Create list for storing ensemble case and death forecasts for all prediction dates
  list_ensemble <- vector("list", length(dates))
  for(i in seq_along(dates)){
    ## Import ensemble predictions
    dt_ensemble_i <- as.data.table(read.csv2(
      paste0("https://raw.githubusercontent.com/covid19-forecast-hub-europe/covid19-forecast-hub-europe/main/data-processed/EuroCOVIDhub-ensemble/", 
             dates[i], "-EuroCOVIDhub-ensemble.csv"),
      sep = ","))
    list_ensemble[[i]] <- dt_ensemble_i[!is.na(quantile),]
  }
  
  # Bind ensemble forecasts for different forecast dates into one data table 
  dt_ensemble <- rbindlist(list_ensemble)
  
  # Remove type column to avoid conflict with type in compar_scores function
  dt_ensemble[, type := NULL]
  
  #### Run calibration ####
  countries <- c("IT", "CZ", "FR")
  if(run){
    if(country == "all"){
      calib <- list()
      ## Generate calibration for each country
      for(i in seq_len(3)){
        country_i <- countries[i]
        if (country_i == "IT" || all_total == T) total <- T else total <- F
        message(sprintf("Running calibration for %s", country_i))
        calib[[i]] <- compute_calib(country_i, range_dates, download, data_file, total, dates_pred, nsim = 100, empty)
      }
      saveRDS(calib, paste0("Output/all_calib_countries", if(all_total) "_total" else "", if(empty) "_empty" else "", ".rds"))
    } else{
      if (country == "IT" || all_total == T) total <- T else total <- F
      calib <- compute_calib(country, range_dates, download, data_file, total, dates_pred, nsim = 100, empty)
      saveRDS(calib, paste0("Output/calib_", country, if(all_total) "_total" else "", if(empty) "_empty" else "", ".rds"))
    }
  } else{
    if(country == "all"){
      calib <- readRDS(paste0("Output/all_calib_countries", if(all_total) "_total" else "", if(empty) "_empty" else "", ".rds"))
    } else {
      calib <- readRDS(paste0("Output/calib_", country, if(all_total) "_total" else "", 
                              if(empty) "_empty" else "", ".rds"))
    }
  }
  
  
  #### Generate Main Tables and Figures ####
  
  #### Generate prediction scores for weekly forecasts from our model and ensemble model
  scores_case <- compar_scores(calib, country, "case", quants, dt_ensemble, countries) 
  scores_case_quantile <- summarise_scores(scores_case$scores_quantile)
  scores_case_point <- scores_case$scores_point
  
  compar_scores_table(scores_case_quantile, country, "case", "quantile", dir_out)
  compar_scores_table(scores_case_point, country, "case", "point", dir_out)
  
  scores_death <- compar_scores(calib, country, "death", quants, dt_ensemble, countries)
  scores_death_quantile <- summarise_scores(scores_death$scores_quantile)
  scores_death_point <- scores_death$scores_point
  
  compar_scores_table(scores_death_quantile, country, "death", "quantile", dir_out)
  compar_scores_table(scores_death_point, country, "death", "point", dir_out)
  
  # Tables comparing scores for forecasts from our model and ensemble model for paper
  table_compar_scores(all_total, empty, "case")
  table_compar_scores(all_total, empty, "death")
  
  if(is.na(min_date)){
    if(country == "all"){
      min_date <- as.Date(max(calib[[1]]$calib$obs %>% names)) - 180  
    } else min_date <- as.Date(max(calib$calib$obs %>% names)) - 180
  }
  
  #### Generate table of prediction scores for all countries at different forecast horizons
  scores_table(calib, country, "case", dir_out)
  scores_table(calib, country, "death", dir_out)
  
  if(country == "all"){
    #### Plot prediction scores
    png(filename = paste0(dir_out, "scores_all_case.png"), width = 1200, height = 300)
    print(figure_scores(calib, "case", which = "crps"))
    dev.off()
    png(filename = paste0(dir_out, "scores_all_death.png"), width = 1200, height = 300)
    print(figure_scores(calib, "death", which = "crps"))
    dev.off()
    
    #### Plot comparison with ensemble forecasts
    figure_country(calib, "case", dt_ensemble, min_date = min_date, dir_out = dir_out)
    figure_country(calib, "death", dt_ensemble, min_date = min_date, dir_out = dir_out)
    
    #### Plot comparison of scores over time with ensemble forecasts
    png(paste0(dir_out, "compar_scores_all_case_quantile.png"), width = 800, height = 700)
    print(figure_compar_scores(scores_case_quantile, "interval_score", "WIS"))
    dev.off()
    png(paste0(dir_out, "compar_scores_all_case_point.png"), width = 800, height = 700)
    print(figure_compar_scores(scores_case_point, "se_point", "SE median"))
    dev.off()
    png(paste0(dir_out, "compar_scores_all_death_quantile.png"), width = 800, height = 700)
    print(figure_compar_scores(scores_death_quantile, "interval_score", "WIS"))
    dev.off()
    png(paste0(dir_out, "compar_scores_all_death_point.png"), width = 800, height = 700)
    print(figure_compar_scores(scores_death_point, "se_point", "SE median"))
    dev.off()
    
    #### Heat maps: Comparison data and median predictions
    png(filename = paste0(dir_out, "case_country_median.png"), width = 1200, height = 900)
    figure_median(calib, "case", min_date = min_date)
    dev.off()
    png(filename = paste0(dir_out, "death_country_median.png"), width = 1200, height = 900)
    figure_median(calib, "death", min_date = min_date)
    dev.off()
    png(filename = paste0(dir_out, "case_country_median_nb.png"), width = 1200, height = 900)
    figure_median(calib, "case", prop = F, min_date = min_date)
    dev.off()
    png(filename = paste0(dir_out, "death_country_median_nb.png"), width = 1200, height = 900)
    figure_median(calib, "death", prop = F, min_date = min_date)
    dev.off()
    
    #### PIT histogram
    figure_pit(calib, "case", min_date = min_date, dir_out = dir_out)
    figure_pit(calib, "death", min_date = min_date, dir_out = dir_out)
  }
  
  #### Generate calibration figures per country ####
  #### Figures of prediction scores
  scores_figure(calib, country, "case", which = "crps", total = all_total, countries = countries, dir_out = dir_out)
  scores_figure(calib, country, "death", which = "crps", total = all_total, countries = countries, dir_out = dir_out)
  
  #### Figures of comparison of prediction scores with ensemble forecasts
  compar_scores_figure(scores_case_quantile, country, "case", "interval_score", "WIS", countries, dir_out)
  compar_scores_figure(scores_case_point, country, "case", "se_point", "SE median", countries, dir_out)
  compar_scores_figure(scores_death_quantile, country, "death", "interval_score", "WIS", countries, dir_out)
  compar_scores_figure(scores_death_point, country, "death", "se_point", "SE median", countries, dir_out)
  
  #### Analysis PIT histograms - cases + deaths
  ## The plots generated will be saved in Output
  pit_figure(calib, country, "case", countries, dir_out)
  pit_figure(calib, country, "death", countries, dir_out)
  
  #### Comparison with ensemble forecasts - case + deaths
  ## The plots generated will be saved in Output
  compar_figure(calib, country, "case", dt_ensemble, countries, dir_out)
  compar_figure(calib, country, "death", dt_ensemble, countries, dir_out)
  
  #### Comparison between median forecasts and data
  compar_median(calib, country, "case", countries, dir_out)
  compar_median(calib, country, "death", countries, dir_out)
  
}
