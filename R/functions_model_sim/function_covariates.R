# - pop_covariate: Generate covariate matrix containing the number of inhabitant per region / age group
# - cumu_covariate: Generate covariate matrix with the cumulative number of cases over the next "delay" days
# - vacc_covariate: Generate covariate matrix containing the proportion of the population vaccinated in the past "waning" days
# - dow_covariate: Generate binary covariate matrix describing the day of the week effect
# - europe_covariate: Generate covariate matrix containing the number of reported cases in the rest of Europe over the past month
# - variant_covariate: Generate binary covariate matrix describing the active variant for each date
# - test_by_age_cov: Generate covariate matrix describing the proportion of the population that got tested in the past two weeks
# - urban_cov: Generate binary covariate matrices describing the urban-rural status of each region 
# - nb_cases: Function creating a 3d matrix containing the number of cases per day, region, and age
# - lag_daily: Lag distribution used in the Endemic-Epidemic model (describing the expected serial interval)
# - all_covariate: Main function generating all covariate matrices, returning a list of matrices.

## Create all covariates
all_covariate <- function(country, range, sts_object, dt_pop, dt_incidence, map, delay_inc, 
                          dt_test, dt_vacc, dt_europe, dt_urban, dt_variant, 
                          dt_vacc_nat = NULL, total = F){
  ## Create time-dependent covariates
  # cumulative incidence, for each variant, over the last 2 days, month, and 6 months
  incidence_new <- cumu_covariate(dt_incidence, c("WT+alpha", "delta", "omicron"), 
                                  sts_object, delay_inc, total)
  incidence_old <- cumu_covariate(dt_incidence, c("WT+alpha", "delta", "omicron"),
                                  sts_object, 30, total)
  incidence_rem <- cumu_covariate(dt_incidence, c("WT+alpha", "delta", "omicron"),
                                  sts_object, 3650, total)
  # New incidence = incidence 0 to 30 days ago
  incidence_new <- incidence_new - incidence_old
  # Old incidence = Incidence 1 to 6 months ago
  incidence_old <- incidence_old - incidence_rem
  
  ## day of the week effect
  dow_cov <- dow_covariate(sts_object, list(tues = "Tuesday", wed = "Wednesday", 
                                            thu = "Thursday", fri = "Friday",
                                            sat = "Saturday", sun = "Sunday"),
                           country)
  ### Vaccine covariates
  ## Regional
  # Second dose in the last 4 months
  dose2_cov <- vacc_covariate(dt_vacc, "dose2", sts_object, country, total = total, waning = 120)
  # Booster dose in the last 4 months
  dose3_cov <- vacc_covariate(dt_vacc, "dose3", sts_object, country, total = total, waning = 2000)
  
  ## Number of cases in Europe covariate
  europe_cov <- europe_covariate(dt_europe, sts_object)
  
  ## Urban / rural status
  list_urban_cov <- urban_cov(dt_urban, sts_object, map, country)
  
  # Variant covariate
  list_var_cov <- variant_covariate(dt_variant = dt_variant, sts_object = sts_object)
  
  # Ensure that the coverage covariate can't go above 1 (set max to 0.999 
  # to avoid errors from log(1 - cov_tot) = log(0))
  covariates <- c(
    list(pop = population(sts_object) / 1e5,
         europe = europe_cov,
         cov_tot = pmin(dose2_cov + dose3_cov, 0.999),
         inc_old = incidence_old, inc_new = incidence_new
    ),
    dow_cov,
    list_urban_cov, list_var_cov)
  
  if (!total){
    age_groups <- unique(dt_incidence$age)
    min_ages <- get_min_age(age_groups)
    nb_age <- length(min_ages)
    # Create covariate: local population by age group
    pop_age <- pop_covariate(dt_pop, dt_incidence, min_ages, tot = F)
    ## Generate population covariate matrix, and case data in hhh4 format
    case_by_age_arr <- nb_cases(dt_incidence)
    pop_by_age_cov <- create_sts(country, case_by_age_arr, pop_age/rowSums(pop_age), 
                                 map, by = "all", flatten = T, timeRange = range, 
                                 agegroups = rep(1,nb_age))@populationFrac
    ## Age group
    GROUPS <- unique(stratum(sts_object, 2))
    ## setup a model matrix with group indicators
    age <- sapply(GROUPS, function (g) {
      index <- which(stratum(sts_object, which = 2) == g)
      res <- col(sts_object)
      res[] <- res %in% index
      res
    }, simplify = FALSE, USE.NAMES = TRUE)
    names(age) <- GROUPS <- 
      c("age0.9", "age10.19", "age20.29", "age30.39", "age40.49", 
        "age50.59", "age60.69", "age70.79", "age80")
    
    ## TAG COUNTRY:
    ## AGE-STRATIFIED MODEL:
    ## Call test_by_age_cov. If testing data is linking to the case data (as in France), 
    ## use dt_incidence as argument, otherwise use dt_test.
    ## Proportion of tests by age group, in the past two weeks
    if (country == "FR"){
      test_age_cov <- test_by_age_cov(country, dt_incidence, sts_object, pop_age, 14)
    } else if (country %in% c("CZ")){
      test_age_cov <- test_by_age_cov(country, dt_test, sts_object, pop_age, 14)
    }
    ## Add pop and test to the covariates object
    covariates <- c(covariates, list(pop_age = pop_by_age_cov), age, test_age_cov)
  } else {
    pop_by_age_cov <- NULL
    
    pop <- matrix(dt_pop[,population], nrow = nrow(dt_pop), ncol = 1)
    ## TAG COUNTRY:
    ## NON-AGE-STRATIFIED MODEL:
    ## Call test_by_age_cov. If testing data is linking to the case data (as in France), 
    ## use dt_incidence as argument, otherwise use dt_test.
    if (country == "FR"){
      test_cov <- test_by_age_cov(country, dt_incidence, sts_object, pop, 14)  
    } else {
      test_cov <- test_by_age_cov(country, dt_test, sts_object, pop, 14)
    }
    covariates <- c(covariates, test_cov)
  }
  
  if (!is.null(dt_vacc_nat)){
    # National
    # Second dose in the last 4 months
    dose2_nat_cov <- vacc_covariate(dt_vacc_nat, "dose2", sts_object, country, national = T, waning = 120)
    # Booster dose in the last 4 months
    dose3_nat_cov <- vacc_covariate(dt_vacc_nat, "dose3", sts_object, country, national = T, waning = 2000)
    
    # Ensure that the _tot coverage covariate can't go above 1 (set max to 0.999 
    # to avoid errors from log(1 - cov_tot) = log(0))
    covariates <- c(
      covariates,list(cov_nat_tot = pmin(dose2_nat_cov + dose3_nat_cov, 0.999))
    )
  }

  return(covariates)
}

## Create pop covariate
pop_covariate <- function(pop_long, dt_incidence, min_ages, tot){
  age_groups <- unique(dt_incidence$age)
  nb_age <- length(age_groups)
  
  # Aggregate populations to age groups in case data
  agg_pop_long <- pop_long[,.(number,age,population)]
  
  agg_pop_long[,age_group := cut(get_min_age(age),c(min_ages,Inf),
                                 labels = age_groups,right = F)]
  
  agg_pop_long <- agg_pop_long[,.(population = sum(population)),by = .(number,age_group)]
  
  ## Create covariate pop
  agg_pop <- dcast(agg_pop_long, number ~ age_group,value.var = "population")
  agg_pop <- as.data.frame(agg_pop)
  rownames(agg_pop) <- agg_pop$number
  agg_pop <- as.matrix(agg_pop[,-1])
  
  if (tot == T) agg_pop <- agg_pop * 0 + rowSums(agg_pop)
  
  return(agg_pop)
}

## Cumulative incidence covariate
cumu_covariate <- function(dt_incidence, wave, sts_object, delay, total = F){
  # Extract columns of interest
  dt_function <- if (total) dt_incidence[, .(date, location_key, variant, cumu_incidence)]
    else dt_incidence[, .(date, location_key, age, variant, cumu_incidence)]
  # Extract dates from sts_object
  dates <- as.Date(epoch(sts_object), origin = "1970-01-01")
  # Last date where the variant "wave" was dominant
  max_date_variant <- max(dt_incidence[is.element(variant, wave), date])
  # Initialise matrix of cumulative incidence
  # Dates used in the loop (excluding cases that may participate in current transmission)
  dates_loop <- as.character(dates - delay)
  dates_loop[dates_loop > max_date_variant] <- max_date_variant
  cumu_incid <- matrix(0, nrow = length(unique(c(dates_loop, dt_function$date))), 
                       ncol = ncol(sts_object@observed))
  colnames(cumu_incid) <- colnames(sts_object@observed)
  rownames(cumu_incid) <- sort(unique(c(dates_loop, dt_function$date)))
  # Vector of labels, used to match dt_function to cumu_incid
  names_out <- if (total){
    apply(dt_function[date == date[1], .(location_key)], 1, 
          function(X) paste(X, collapse = "."))}
  else {
    apply(dt_function[date == date[1], .(location_key, age)], 1, 
          function(X) paste(X, collapse = "."))}
  # Select entries corresponding to the variant and dates of interest
  dt_function <- dt_function[is.element(variant, wave),]
  # For each date, extract the cumulative incidence up to date_i - delay
  for(i in seq_along(unique(dt_function$date))){
    date_i <- unique(dt_function$date)[i]
    out <- dt_function[date == date_i, cumu_incidence] 
    if (length(out) == 0) out <- rep(0, ncol(sts_object@observed))
    names(out) <- names_out
    cumu_incid[date_i, names(out)] <- out
  }
  if(length(wave) > 1){
    variant_change <- which(diff(cumu_incid[,1]) < 0)
    for(i in rev(variant_change)){
      cumu_incid[(i + 1):nrow(cumu_incid),] <- t(t(cumu_incid[(i + 1):nrow(cumu_incid),]) +
        cumu_incid[i,])
    }
  }
  cumu_incid <- cumu_incid[dates_loop,]
  # Return incidence per 100,000
  return(cumu_incid / 1e5)
}

## Coverage covariates
vacc_covariate <- function(dt_vaccine, which_dose, sts_object, country, waning = NULL, total = F, national = F){
  if (is.null(dt_vaccine)){
    out_vacc <- sts_object@observed * 0
  }  else {
    # Extract entries and columns of interest
    dt_function <- dt_vaccine[dose == which_dose,]
    ## Create ID column using date, location and age, to match dt_function and sts_object@observed
    if (total){
      if (national){
        dt_function[, ID := date_week]
      } else {
        dt_function[, ID := paste(date_week, region_nb, sep = ".")]
      }
    } else {
      if (national){
        dt_function[, ID := paste(date_week, age, sep = ".")]  
      } else {
        dt_function[, ID := paste(date_week, region_nb, age, sep = ".")]    
      }
    }
    setkey(dt_function, ID)
    # Initialise the return matrix
    out_vacc <- sts_object@observed * 0
    
    # Remove 2 weeks for the dose to be effective
    dates <- as.Date(epoch(sts_object), origin = "1970-01-01") - 14
    if (dt_function[,max(diff(as.integer(as.Date(date_week))))] == 1) {
      # If vaccination data is daily, just use it as it is
      week_dates <- dates    
    } else {
      # Value used to match the (daily) incidence to (weekly) vaccine uptake
      offset <- (as.Date(dt_function$date_week[1]) - 
                   lubridate::floor_date(as.Date(dt_function$date_week[1]), "weeks", 1)) %>% 
        as.numeric %%7
      week_dates <- lubridate::floor_date(dates, "weeks", 1 + offset)
    }
    # Extract region / age identified from sts_object
    nms <- colnames(sts_object@observed)
    dep_nb_sts <- if (national){
      get_age(nms, country)
    } else {
      if (total){
        get_reg_nb(nms, country)
      } else {
        get_reg_age_nb(nms, country)
      }
    }
    dates_vacc <- as.Date(unique(dt_function$date_week))
    nb_reg <- length(unique(get_reg_nb(nms, country)))
    
    if (is.null(waning)){
      # If waning is not considered, just return the column cumu_dose for each date
      for(i in seq_along(week_dates))
        if (total){
          if (national){
            out_vacc[i,] <- dt_function[as.character(week_dates[i]), cumu_dose]
          } else {
            out_vacc[i,] <- dt_function[paste0(week_dates[i], ".", dep_nb_sts), cumu_dose]
          }
        } else {
          if (national){
            out_vacc[i,] <- dt_function[paste0(week_dates[i], ".", dep_nb_sts), cumu_dose]
          } else {
            out_vacc[i,] <- dt_function[paste0(week_dates[i], ".", dep_nb_sts), cumu_dose]
          }
        }
    } else{
      # Otherwise, compute the cumulative proportion of doses over the waning period
      for(i in seq_along(week_dates)){
        dates_vacc_i <- 
          as.character(dates_vacc[dates_vacc <= week_dates[i] & 
                                    dates_vacc >= week_dates[i] - waning])
        if (any(dt_function[,date_week] %in% dates_vacc_i)){
          out_vacc[i,] <- if (total){
            if (national){
              dt_function[is.element(date_week, dates_vacc_i), lapply(.SD, sum), 
                          .SDcols = c("prop_dose")][
                            , prop_dose] 
            } else {
              dt_function[is.element(date_week, dates_vacc_i), lapply(.SD, sum), 
                          by = .(region_nb), .SDcols = c("prop_dose")][
                            order(region_nb), prop_dose]             
            }
          } else {
            if (national){
              rep(dt_function[is.element(date_week, dates_vacc_i), lapply(.SD, sum), 
                              by = .(age), .SDcols = c("prop_dose")][
                                order(age), prop_dose], each = nb_reg)
            } else {
              dt_function[is.element(date_week, dates_vacc_i), lapply(.SD, sum), 
                          by = .(region_nb, age), .SDcols = c("prop_dose")][
                            order(age, region_nb), prop_dose]            
            }
          }
        }
      }
    }
    
    # Overwrite any missing values with 0
    out_vacc[is.na(out_vacc)] <- 0    
  }

  return(out_vacc)
}

## day of the week 
dow_covariate <- function(sts_object, groups, country){
  # Compute the weekday of each row in sts_object@observed
  dow <- matrix(weekdays(as.Date(epoch(sts_object), origin = "1970-01-01")),
                ncol = ncol(sts_object@observed), nrow = nrow(sts_object@observed))
  ## TAG COUNTRY
  # Define the vector bank_holiday, which contains the dates of all bank holidays in the 
  # country in 2021 and 2022.
  bank_holiday <- bank_holiday_country(country)
  
  # Set bank holidays as Sundays
  if(country == "IT"){
    dow[is.element(as.Date(epoch(sts_object), origin = "1970-01-01"), 
                   as.Date(bank_holiday)),] <- "Monday"
  } else{
    dow[is.element(as.Date(epoch(sts_object), origin = "1970-01-01"), 
                   as.Date(bank_holiday)),] <- "Sunday"
  }
  colnames(dow) <- colnames(sts_object@observed)
  # Match dow to each element listed in groups
  dow_cov <- lapply(groups, function(X) {
    out <- matrix(0, nrow = nrow(dow), ncol = ncol(dow))
    colnames(out) <- colnames(dow)
    out[] <- is.element(dow, X)
    return(out)
  })
  names(dow_cov) <- names(groups)
  return(dow_cov)
}

bank_holiday_country <- function(country){
  if (country == "FR"){
    bank_holiday <- c("2021-01-01", "2021-04-04", "2021-04-05", "2021-05-01", "2021-05-08", "2021-05-13", 
                      "2021-05-24", "2021-07-14", "2021-08-15", "2021-11-01", "2021-11-11", "2021-12-25",
                      "2022-01-01", "2022-04-17", "2022-04-18", "2022-05-01", "2022-05-08", "2022-05-26",
                      "2022-06-06", "2022-07-14", "2022-08-15", "2022-11-01", "2022-11-11", "2022-12-25",
                      "2023-01-01", "2023-04-07", "2023-04-10", "2023-05-01", "2023-05-08", "2023-05-18",
                      "2023-05-28", "2023-05-29", "2023-07-14", "2023-08-15", "2023-11-01", "2023-11-11", 
                      "2023-12-25")
  } else if (country == "CZ"){
    bank_holiday <- c("2021-01-01", "2021-04-02", "2021-04-05", "2021-05-01", "2021-05-08", "2021-07-05",
                      "2021-07-06", "2021-09-28", "2021-10-28", "2021-11-17", "2021-12-25", "2022-01-01",
                      "2022-04-15", "2022-04-18", "2022-05-01", "2022-05-08", "2022-07-05", "2022-07-06",
                      "2022-09-28", "2022-10-28", "2022-11-17", "2022-12-25", "2023-01-01", "2023-04-07", 
                      "2023-04-10", "2023-05-01", "2023-05-08", "2023-07-05", "2023-07-06", "2023-08-21",
                      "2023-09-28", "2023-10-28", "2023-11-17", "2023-12-25")
  } else if (country == "IT"){
    bank_holiday <- c("2021-01-01", "2021-01-06", "2021-04-04", "2021-04-05", "2021-04-25", "2021-05-01", 
                      "2021-06-02", "2021-08-15", "2021-11-01", "2021-12-08", "2021-12-25", "2022-01-01", 
                      "2022-01-06", "2022-04-17", "2022-04-18", "2022-04-25", "2022-05-01", "2022-06-02",
                      "2023-08-15", "2023-11-01", "2023-12-08", "2023-12-25", "2023-01-01", "2023-01-06", 
                      "2023-04-09", "2023-04-10", "2023-04-25", "2023-05-01", "2023-06-02", "2023-08-15", 
                      "2023-11-01", "2023-12-08", "2023-12-25", "2023-12-26")
  }
  return(bank_holiday)
}

## Number of cases in Europe in the past month
europe_covariate <- function(dt_europe, sts_object, delay = 30){
  # Extract dates from sts_object
  dates <- as.Date(epoch(sts_object), origin = "1970-01-01")
  # Initialise matrix
  europe_cov <- matrix(, nrow = nrow(observed(sts_object)), ncol = ncol(observed(sts_object)))
  # For each date, compute the number of cases reported in Europe in the past 30 days
  for(i in seq_along(dates)){
    date_max <- dates[i] - 1
    date_min <- dates[i] - delay
    europe_cov[i,] <- sum(dt_europe[date > date_min & date <= date_max, cases])
  }
  return(europe_cov)
}

## Function to create a 3d matrix from dt_incidence. This matrix contains the 
## number of cases per day, region, and age
nb_cases <- function(dt_incidence){
  # Extract unique values of dates, regs and age_groups
  dates <- unique(dt_incidence[order(date),]$date)
  regs <- unique(dt_incidence[order(region),]$location_key) 
  age_groups <- unique(dt_incidence[order(age),]$age)
  # Compute the number of unique date, reg, and age_group
  nb_date <- length(dates)
  nb_reg <- length(regs)
  nb_age <- length(age_groups)
  # Make 3D array of case counts by date (rows) x region (columns) x age group (layers) 
  case_by_age_arr <- array(dt_incidence[order(region, date, age),nb_cases],
                              dim = c(nb_age, nb_date, nb_reg),
                              dimnames = list(age_groups, dates, regs))
  case_by_age_arr <- aperm(case_by_age_arr, c(2, 3, 1))
  return(case_by_age_arr)
}

## Lag distribution
lag_daily <- function(par_lag, min_lag, max_lag){
  prop1 <- 1/(1+exp(-par_lag))
  sd_w1 <- 1.5
  # Generate the distribution of the serial interval (no missing case)
  w1 <- dnorm(x = 0:max_lag, mean = 5, sd = sd_w1)
  w1 <- w1/sum(w1)
  # Normalise distribution (>= 1 day)
  # Generate the distribution of the serial interval (missing generation)
  w2 <- conv(w1, max_days = max_lag)
  w2[is.na(w2)] <- 0
  w1 <- w1[-1]
  # Compute the composite serial interval
  w_dens <- w1 * prop1 + w2 * (1-prop1)
  weights0 <- c(rep(0, min_lag - 1), w_dens)
  weights <- weights0/sum(weights0)
  return(weights)
}

## Test by age
test_by_age_cov <- function(country, dt_incidence, sts_object, dt_pop, delay){
  dt_function <- copy(dt_incidence)
  dt_function[, date := as.Date(date)]
  # Extract all dates
  dates <- as.Date(epoch(sts_object), origin = "1970-01-01")
  # Extract the number of regions in sts_object@observed
  nb_reg <- length(unique((sub("[.].*", "", colnames(sts_object@observed)))))
  # Initialise the output matrix
  prop_test_cov <- prop_age <- sts_object@observed * 0
  for(i in seq_along(dates)){
    # For each date, compute the number of test in the past "delay" days by age group
    # across the country
    date_i <- dates[i]
    if(any(colnames(dt_function) == "location_key")){
      nb_test_i <- dt_function[date <= date_i & date > date_i - delay, lapply(.SD, sum), 
                               by = .(age, region, location_key), .SDcols = "nb_tests"]
      vec_pop <- c(dt_pop)
      names(vec_pop) <- paste(rep(rownames(dt_pop), ncol(dt_pop)), 
                              rep(colnames(dt_pop), each = nrow(dt_pop)), sep = "_")
      nb_test_i[, prop_test := nb_tests / vec_pop[paste(region, age, sep = "_")]]
      vec_nb_test <- nb_test_i$prop_test
      if(all(nb_test_i$age == "Total")) names(vec_nb_test) <- nb_test_i$location_key else 
        names(vec_nb_test) <- paste(nb_test_i$location_key, nb_test_i$age, sep = ".")
      prop_age[i, names(vec_nb_test)] <- vec_nb_test
    } else {
      nb_test_i <- dt_function[date <= date_i & date > date_i - delay, lapply(.SD, sum), 
                               by = .(age), .SDcols = "nb_tests"]
      prop_test_i <- nb_test_i$nb_tests / colSums(dt_pop)
      
      prop_age[i, ] <- rep(prop_test_i, each = nb_reg)
    }
    
    # Divide it by the number of inhabitants per age group
    tot_test_i <- sum(nb_test_i$nb_tests)
    prop_test_cov[i, ] <- tot_test_i/sum(dt_pop)
  }
  prop_test_cov[is.na(prop_test_cov)] <- 0
  prop_age[is.na(prop_age)] <- 0
  return(list(test_prop = prop_test_cov, test_age = prop_age))
}

## Variant period
variant_covariate <- function(dt_variant, sts_object, prop_delta = 0.3, prop_omicron = 0.3){
  # Get dates from sts object
  dates <- as.Date(epoch(sts_object), origin = "1970-01-01")
  nb_date <- length(dates)
  # Compute beginning of Delta and Omicron period
  start_delta <- dt_variant[variant == "Delta" & percent_variant > prop_delta, date][1]
  start_omicron <- dt_variant[variant == "Omicron" & percent_variant > prop_omicron &
                                date > start_delta, date][1]
  # Binary vectors, equal to 1 when a given variant is dominant
  WT_alpha <- as.integer(dates < start_delta)
  delta <- as.integer(dates >= start_delta & dates < start_omicron)
  omicron <- as.integer(dates >= start_omicron)
  # Create covariate matrices
  WT_alpha_mat <- matrix(WT_alpha, nrow = nb_date, ncol = ncol(sts_object))
  delta_mat <- matrix(delta, nrow = nb_date, ncol = ncol(sts_object))
  omicron_mat <- matrix(omicron, nrow = nb_date, ncol = ncol(sts_object))
  
  return(list(WT_alpha = WT_alpha_mat, delta = delta_mat, omicron = omicron_mat))
}

# Create urban / rural covariate
urban_cov <- function(rural, sts_object, map, country){
  # Import rural / urban status (from https://ec.europa.eu/eurostat/web/rural-development/methodology)
  dep_vac <- map$reg_nb
  names(dep_vac) <- map$NUTS_ID
  rural[, dep := dep_vac[NUTS]]
  rural <- rural[!is.na(dep),]
  #1 Urban
  #21 intermediate - urban
  #31 Intermediate - rural
  #32 Rural
  
  rural_cov <- sts_object@observed * 0
  int_rur_cov <- sts_object@observed * 0
  int_urb_cov <- sts_object@observed * 0
  col_nms <- get_reg_nb(colnames(rural_cov), country)
  rural_cov[, is.element(col_nms, rural[code == 32, dep])] <- 1
  int_rur_cov[, is.element(col_nms, rural[code == 31, dep])] <- 1
  int_urb_cov[, is.element(col_nms, rural[code == 21, dep])] <- 1
  
  return(list(rural = rural_cov, int_rur = int_rur_cov, 
              int_urb = int_urb_cov))
}
