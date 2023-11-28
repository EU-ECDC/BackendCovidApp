generate_list_output <- function(country, range_dates, download, data_file, total, nsim = 1000, 
                                 last_date = FALSE, fit = TRUE, simulations = TRUE){
  delay_cumu <- 0
  #### Import and process data
  ## Import all data, make sts object, and generate covariates
  list_import <- import_all_files(country = country, total = total, 
                                  range = range_dates, delay_cumu = delay_cumu,
                                  download_case = download, file_case = data_file,
                                  last_date = last_date)
  
  ## Extract elements of list_import
  sts_obj_nei <- list_import$sts_obj_nei
  map <- list_import$map
  age_groups <- list_import$age_groups
  covariates <- list_import$covariates
  C <- list_import$C
  # If last_date is TRUE, fit to the last available date
  if(last_date == T){
    last_report <- max(as.Date(epoch(list_import$sts_obj_nei), origin = "1970-01-01"))
    range_dates[2] <- as.character(last_report)
  }
  
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
    end_terms <- c("1", "log(europe/1e8)", "log(pop * pop_age)", 
                   "rural", "int_rur", "int_urb"
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
  equations <- 
    list(ar = ar_terms,
         ne = addSeason2formula(reformulate(ne_terms, intercept = TRUE), period = 365),
         end = addSeason2formula(reformulate(end_terms, intercept = TRUE))
    )
  ## Run the model
  model_fit <-
    hhh4_cov(sts_object = sts_obj_nei, contact = C/max(C), equations = equations,
             covariates = covariates, lag_dist = lag_daily, parameter_lag = par_dist,
             weight_fun = W_powerlaw(maxlag = 5, log = T, normalize = T))
  
  #### Generate and analyse predictions
  ### Initialise list of predictions
  pred_list <- list()
  pred_list[[1]] <- list()
  pred_list[[2]] <- list()
  
  if(fit){
    ### 28-day forecasts for the past three months
    t_pred <- nrow(sts_obj_nei@observed) - seq(0, 90, 7)
    for(i in seq_along(t_pred)){
      ## Prediction date
      t_i <- t_pred[i]
      model_pred <- model_fit
      ## Set extreme values of standard deviations to 0
      if(any(model_pred$se > 10)) {
        model_pred$coefficients[-length(model_pred$coefficients)] <- 0
        model_pred$coefficients[length(model_pred$coefficients)] <- 2
      }
      set.seed(1)
      ## Generate predictions
      pred_list[[1]][[i]] <- nStepAhead(model = model_pred, n = 28, nsim = nsim/10, t_start = t_i, 
                                        country = country, delay = delay_cumu, t_end = t_i, total = total, nparam = 10)
    }
    names(pred_list[[1]]) <- t_pred
    pred <- pred_list[[1]][[1]]
  }
  
  if(simulations){
    ### Generate scenarios with changes in transmission
    delay <- c(7)
    transmissibility <- c(1, 1.2, 1.4)
    NPIs <- c(.6, .8, 1)
    if (!total) target <- c("children", "work", "older", "all", "endemic") else 
      target <- c("all", "endemic")
    
    count <- 1
    
    ### Generate predictions for each value of changes in transmission, target groups, and delay
    for(j in seq_along(transmissibility)){
      t_i <- nrow(sts_obj_nei@observed)
      for(k in seq_along(NPIs)){
        for(l in seq_along(target)){
          ## Select age group targeted by drop in transmission
          if (target[l] == "children") cols <- c(grep("0-9", colnames(model_fit$stsObj@observed)),
                                                 grep("10-19", colnames(model_fit$stsObj@observed)))
          if (target[l] == "older") cols <- c(grep("60-69", colnames(model_fit$stsObj@observed)),
                                              grep("70-79", colnames(model_fit$stsObj@observed)),
                                              grep("80+", colnames(model_fit$stsObj@observed)))
          if (target[l] == "work") cols <- c(grep("20-29", colnames(model_fit$stsObj@observed)),
                                             grep("30-39", colnames(model_fit$stsObj@observed)),
                                             grep("40-49", colnames(model_fit$stsObj@observed)),
                                             grep("50-59", colnames(model_fit$stsObj@observed)))
          if (target[l] == "all" | target[l] == "endemic") cols <- seq_len(ncol(model_fit$stsObj@observed))
          ## If endemic, remove imports
          if (target[l] == "endemic") import <- 0 else import <- 1
          
          for(m in seq_along(delay)){
            ## Generate offset transmissibility matrix, corresponding to the change in transmission
            ## in this scenario
            trans_mat <- matrix(1, ncol = model_fit$nUnit, 
                                nrow = nrow(model_fit$control$data$pop) + 28)
            colnames(trans_mat) <- colnames(model_fit$control$data$pop)
            trans_mat[-(1:t_i),] <- transmissibility[j]
            trans_mat[-(1:(t_i + delay[m])), cols] <- transmissibility[j] * NPIs[k]
            
            ## Generate predictions
            set.seed(1)
            pred_list[[2]][[count]] <- nStepAhead(
              model = model_fit, n = 28, nsim = nsim/10, t_start = t_i, delay = delay_cumu, 
              t_end = t_i, total = total, transmissibility =  trans_mat, importation = import,
              nparam = 10, country = country)
            ## Set names of current scenario
            names(pred_list[[2]])[count] <- paste(transmissibility[j], NPIs[k], target[l], 
                                                  delay[m], sep = ";")
            count <- count + 1
          }
        }
      }
    }
  }
  
  #### Generate predictions number of deaths
  dt_cfr <- calculate_cfr(country = country, download = T, data_file = data_file, 
                          death_data_file = "", total = total, level = "nuts2")
  min_date_omicron <- min(as.Date(epoch(sts_obj_nei), origin = "1970-01-01")[covariates$omicron[,1] == 1])
  ## Prediction number of deaths
  death_pred <- forecast_death(predictions = pred_list, 
                               dt_cfr = dt_cfr[date > min_date_omicron & !is.na(CFR),], 
                               cases = sts_obj_nei@observed, country = country,
                               dates = as.Date(epoch(sts_obj_nei), origin = "1970-01-01"),
                               total = total)
  ## Import data on number of deaths
  data_death_mat <- data_death(country = country, download = T, death_data_file = "", 
                               total = total, pred_date = range_dates[2])
  
  if(all(substr(colnames(data_death_mat), 1, 2)[1] != country)) 
    colnames(data_death_mat) <- paste0(country, "_", colnames(data_death_mat))
  ## Extract values of predictors
  predictors <- meanHHH(coef(model_fit), terms(model_fit))
  
  ## Extract wide age groups
  if (!total){
    age_groups <- sub(".*[.]", "", colnames(sts_obj_nei))
    age_groups[age_groups == "0-9" | age_groups == "10-19"] <- "0-20"
    age_groups[is.element(age_groups, c("20-29", "30-39", "40-49", "50-59"))] <- "20-60"
    age_groups[is.element(age_groups, c("60-69", "70-79"))] <- "60-80"
    age_groups[is.element(age_groups, c("80-89", "90+"))] <- "80+"
    
    age_groups_death <- sub(".*[.]", "", colnames(data_death_mat))
    age_groups_death[age_groups_death == "0-9" | age_groups_death == "10-19"] <- "0-20"
    age_groups_death[is.element(age_groups_death, c("20-29", "30-39", "40-49", "50-59"))] <- "20-60"
    age_groups_death[is.element(age_groups_death, c("60-69", "70-79"))] <- "60-80"
    age_groups_death[is.element(age_groups_death, c("80-89", "90+"))] <- "80+"
  } else {
    age_groups_death <- rep(age_groups, ncol(data_death_mat))
    age_groups <- rep(age_groups, ncol(sts_obj_nei))
  }
  
  ## Extract NUTS-2 and NUTS-3 region names
  reg_groups_nuts2 <- get_nuts2_reg(colnames(sts_obj_nei), country)
  reg_groups_nuts3 <- get_nuts3_reg(colnames(sts_obj_nei), country)
  
  ## Create different grouping of columns (by country, nuts2, nuts3; age stratified or not)
  if (!total){
    groups <- list(country_age = paste0(country, ".", age_groups), 
                   country_tot = rep(paste0(country, ".tot"), length(age_groups)),
                   nuts2_age = paste0(reg_groups_nuts2, ".", age_groups), 
                   nuts2_tot = paste0(reg_groups_nuts2, ".tot"))
    
    groups_death <- list(country_age = paste0(country, ".", age_groups_death), 
                         country_tot = rep(paste0(country, ".tot"), length(age_groups_death)),
                         nuts2_age = paste0(sub("[.].*", "", colnames(data_death_mat)), ".", 
                                            age_groups_death), 
                         nuts2_tot = paste0(sub("[.].*", "", colnames(data_death_mat)), ".tot"))
    
    if (!is.null(reg_groups_nuts3)){
      groups <- c(groups, 
                  list(nuts3_age = paste0(reg_groups_nuts3, ".", age_groups), 
                       nuts3_tot = paste0(reg_groups_nuts3, ".tot")))
    }
  } else {
    groups <- list(country = rep(paste0(country,".tot"), ncol(sts_obj_nei)),
                   nuts2_tot = paste0(reg_groups_nuts2, ".tot"))
    groups_death <- list(country = rep(paste0(country,".tot"), ncol(data_death_mat)),
                         nuts2_tot = paste0(sub("[.].*", "", colnames(data_death_mat)), 
                                            ".", age_groups_death, ".tot"))
    if (!is.null(reg_groups_nuts3)){
      groups <- c(groups, list(nuts3_tot = paste0(reg_groups_nuts3, ".tot")))
    }
  }
  
  ## Generate summary statistics of pred_list objects
  if(fit == TRUE){
    pred_out1 <- format_pred(list_pred = pred_list[[1]], groups = groups, deaths = death_pred[[1]],
                             groups_death = groups_death)
  } else if(simulations == TRUE){
    pred_out1 <- readRDS(paste0("Output/output_model_", country, ifelse(total, "_total", ""), "_nosim.RDS"))$pred
  } else pred_out1 <- NULL

  if(simulations == TRUE){
    pred_out2 <- format_pred(list_pred = pred_list[[2]], groups = groups, deaths = death_pred[[2]],
                             groups_death = groups_death)
  } else pred_out2 <- NULL
  ## Build list output object
  list_output <- list(map = map, pred = pred_out1, obs = model_fit$stsObj,
                      scenario = pred_out2, pop = model_fit$control$data$pop, 
                      predictors = predictors, pop_age = model_fit$control$data$pop_age,
                      data_death = data_death_mat)
  return(list_output)
}
