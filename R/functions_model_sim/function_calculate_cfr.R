# calculate_cfr:  Compute Case Fatality Ratio from case and death data. Can be 
# generated at a local level ("nuts2", by default), or at a national level ("national")
# data_death: Import death data, and output a matrix containing the number of 
# death reported per week and region
# regress_cfr: Implement linear regression models linking changes in CFR to 
# changes in the number of cases 1, 2, and 3 weeks ahead, using the data frame 
# generated in calculate_cfr.
# forecast_death: Generate 1, 2, 3, and 4 week-ahead forecasts, using the case data, and
# the parameters estimated in regress_cfr.

# CFR computed with a 3-week shift between case and death
# then summing the number of cases and shifted deaths in the last week
calculate_cfr <- function(country, download, data_file, death_data_file, total = F, level = "nuts-2"){
  
  # Import index
  index <- import_index()
  
  if (total){
    # Import  cases
    case <- import_case(country, download, data_file, index, total)
    setnames(case, "new_confirmed", "cases")
    # Import deaths
    death <- import_death(country, download, death_data_file, index, total)
    setnames(death, "new_deceased", "deaths")
    
    ## TAG COUNTRY
    ## if case data is at NUTS3 level and death data is at NUTS2, aggregate case
    ## data to match the spatial granularity of death data.
    # If country is Italy or France, aggregate case data at NUTS2 level
    if (country == "IT"){
      case[, location_key := get_nuts2_reg(location_key, country)]
      case <- case[, .(cases = sum(cases)), by = .(location_key, date)]
    } else if (country == "FR"){
      case[, location_key := substr(location_key, 1, 6)]
      case <- case[, .(cases = sum(cases)), by = .(location_key, date)]
    }
    
    # Merge cases and death data
    case_and_death <- merge(case, death, by = c("location_key", "date"))
    
    # Shift deaths back by 3 weeks to account for delay between case confirmation and death
    case_and_death[, deaths_shifted := shift(deaths, n = 21, type = "lead"), by = .(location_key)]
    cols <- c("cases", "deaths", "deaths_shifted")
    # Add NUTS2 region codes
    case_and_death[, nuts2 := get_nuts2_reg(location_key, country)]
    if (level == "nuts2"){
      # Aggregate cases and deaths by NUTS2 region
      case_and_death_agg <- case_and_death[, lapply(.SD, sum), .SDcols = cols, by = .(date, nuts2)]
      # Calculate CFR
      case_and_death_agg[, CFR := frollsum(deaths_shifted, n = 7)/frollsum(cases, n = 7), by = .(nuts2)]
    } else if (level == "national"){
      # Aggregate cases and deaths at national level
      case_and_death_agg <- case_and_death[, lapply(.SD, sum), .SDcols = cols, by = .(date)]
      # Calculate CFR
      case_and_death_agg[, CFR := frollsum(deaths_shifted, n = 7)/frollsum(cases, n = 7)]
    }
  } else {
    # Import cases
    case_by_age <- import_case(country, download, data_file, index, total)
    # Import deaths
    death_by_age <- import_death(country, download, death_data_file, index, total)
    
    # Melt to long format
    case_by_age_long <- melt(case_by_age, id.vars = c("location_key", "date"), measure.vars = patterns(cases = "new_confirmed", age_group = "age_bin"))
    ## TAG COUNTRY
    ## if case data is at NUTS3 level and death data is at NUTS2, aggregate case
    ## data to match the spatial granularity of death data.
    # If country is France, aggregate cases at NUTS2 level
    if (country == "FR"){
      case_by_age_long[, location_key := substr(location_key, 1, 6)]
      case_by_age_long <- case_by_age_long[, .(cases = sum(cases)), by = .(location_key, date, variable, age_group)]
    }
    death_by_age_long <- melt(death_by_age, measure.vars = patterns(deaths = "new_deceased", age_group = "age_bin"))
    
    # Merge case and death data
    case_and_death_long <- merge(case_by_age_long[, !"variable"], death_by_age_long[, !"variable"], by = c("location_key", "date", "age_group"))
    
    # Shift deaths back by 3 weeks to account for delay between case confirmation and death
    case_and_death_long[, deaths_shifted := shift(deaths, n = 21, type = "lead"), by = .(location_key, age_group)]
    
    # case_and_death_long[, year_mon := format(date, "%Y-%m")]
    
    cols <- c("cases", "deaths", "deaths_shifted")
    # Add NUTS2 region codes
    case_and_death_long[, nuts2 := get_nuts2_reg(location_key, country)]
    if (level == "nuts3"){
      case_and_death_agg <- copy(case_and_death_long)
      # Calculate CFR
      case_and_death_agg[, CFR := frollsum(deaths_shifted, n = 7)/frollsum(cases, n = 7), by = .(location_key, age_group)]
    } else if (level == "nuts2"){
      # Aggregate cases and deaths by NUTS2 region
      case_and_death_agg <- case_and_death_long[, lapply(.SD, sum), .SDcols = cols, by = .(date, nuts2, age_group)]
      case_and_death_agg[, CFR := frollsum(deaths_shifted, n = 7)/frollsum(cases, n = 7), by = .(nuts2, age_group)]
    } else if (level == "national"){
      # Aggregate cases and deaths at national level
      case_and_death_agg <- case_and_death_long[, lapply(.SD,  sum), .SDcols = cols, by = .(date, age_group)]
      case_and_death_agg[, CFR := frollsum(deaths_shifted, n = 7)/frollsum(cases, n = 7), by = .(age_group)]
    }
  }
  case_and_death_agg[!is.na(deaths_shifted) & is.nan(CFR), CFR := 0]
  
  return(case_and_death_agg)
}

# Forecast number of deaths, using the case data and CFR computed in calculate_cfr
forecast_death <- function(predictions, dt_cfr, cases, dates, country, total = F){
  # Apply the same function to each element of predictions
  predictions_week <- lapply(predictions, function(X){
    predictions_week <- list()
    for(i in seq_along(X)){
      # Access nb of cases before predictions
      if(grepl("^[0-9]{1,}$", names(X)[i])) pred_date_num <- as.numeric(names(X)[i]) else
        pred_date_num <- nrow(cases)
      date_i <- dates[pred_date_num]
      
      # Compute linear regression models linking changes in cfr to changes in the number
      # of cases 1, 2, and 3 weeks ahead
      cfr_i <- regress_cfr(dt_cfr, date_i)
      pred_i <- X[[i]]$sim
      
      # Aggregate case numbers by week (for death predictions 1, 2, and 3 weeks ahead)
      cases_week <- rbind(colSums(cases[pred_date_num - 14:20,]),
                          colSums(cases[pred_date_num - 7:13,]),
                          colSums(cases[pred_date_num - 0:6,]))
      
      # Aggregate predictions per week (for death predictions 4 weeks ahead)
      pred_week <- array(data = c(colSums(pred_i[1:7,,])
      ), dim = c(dim(pred_i)[c(2,3)], 1)) 
      
      ## Aggregate predictions by NUTS2 region
      # Get NUTS2 region for each column of cases_week
      col_nuts2 <- get_nuts2_reg(colnames(cases_week), country)
      
      # If the model is age-stratified, add the age group in each column 
      if(!total) col_nuts2 <- paste0(col_nuts2, ".", sub(".*[.]", "", colnames(cases_week)))
      
      # Aggregate case numbers (from data)
      cases_week <- vapply(unique(col_nuts2), function(x)
        rowSums(cases_week[,col_nuts2 == x, drop = FALSE]), numeric(nrow(cases_week)))
      
      # Aggregate predicted case numbers
      aggreg_pred <- vapply(unique(col_nuts2), function(x){
        mat <- pred_week[,,1]
        colSums(mat[col_nuts2 == x, , drop = FALSE])
      }, numeric(dim(pred_week)[2]))
      pred_week <- array(data = t(aggreg_pred), 
                         dim = c(sum(!duplicated(col_nuts2)), dim(pred_week)[c(2,3)]))
      
      ## Generate death forecasts
      # Create empty array (will be filled with death forecasts)
      pred_deaths <- array(, dim = c(4, nrow(pred_week), ncol(pred_week)))
      rownames(pred_deaths) <- pred_date_num + 7 * c(1, 2, 3, 4)
      colnames(pred_deaths) <- colnames(cases_week)
      # Draw forecasts from binomial distribution
      for(j in seq_len(ncol(cases_week))){
        # Select CFR in the NUTS2 region of interest
        if(any(colnames(cfr_i) == "nuts2")){
          nuts2_j <- sub("[.].*", "", colnames(cases_week)[j])
          if(total) age_j <- "Total" else age_j <- sub(".*[.]", "", colnames(cases_week)[j])
          # cfr_set contains the entries of cfr_i in the region / age-group of interest
          cfr_set <- cfr_i[nuts2 == nuts2_j & age_group == age_j, ]
          # ref_cfr corresponds to the last known value of the cfr in this region / age-group
          ref_cfr <- cfr_set[date == (max(as.Date(date)) - 21), CFR]
        } else cfr_set <- cfr_i$CFR

        ## Draw forecasts
        for(k in seq_len(nrow(cases_week))){
          # Draw values of cfr using the mean and standard deviation of delta_cfr
          cfr_vec <- cfr_set[date == (min(as.Date(date)) + 7 * (k)),
                             ref_cfr + rnorm(n = 10, mean = delta_cfr, sd = sd_cfr)]
          # If predicted values of the cfr are below 0, set them to 0
          cfr_vec[cfr_vec < 0] <- 0
          
          # Draw forecasts using a binomial distribution
          drawjk <- rbinom(n = ncol(pred_week), size = cases_week[k, j], prob = cfr_vec)
          pred_deaths[k, j, ] <- drawjk
        }
      }
      ## To generate 4-week ahead predictions, use cfr from 3-week ahead model
      for(j in seq_len(nrow(pred_week))){
        # Select CFR in the NUTS2 region
        if(any(colnames(cfr_i) == "nuts2")){ 
          nuts2_j <- sub("[.].*", "", colnames(cases_week)[j])
          if(total) age_j <- "Total" else age_j <- sub(".*[.]", "", colnames(cases_week)[j])
          cfr_set <- cfr_i[nuts2 == nuts2_j & age_group == age_j, ]
          ref_cfr <- cfr_set[date == (max(as.Date(date)) - 21), CFR]
        } else cfr_set <- cfr_i$CFR

        # Draw forecasts
        for(k in seq_len(dim(pred_week)[3])){
          n_cases_week <- nrow(cases_week)
          cfr_vec <- cfr_set[date == (min(as.Date(date)) + 7 * (n_cases_week)),
                             ref_cfr + rnorm(n = 10, mean = delta_cfr, sd = sd_cfr)]
          cfr_vec[cfr_vec < 0] <- 0
          drawjk <- rbinom(n = ncol(pred_week), size = pred_week[j, , k], prob = cfr_vec)
          pred_deaths[k + nrow(cases_week), j, ] <- drawjk
        }
      }
      predictions_week[[i]] <- pred_deaths
    }
    names(predictions_week) <- names(X)
    return(predictions_week)
  })
  return(predictions_week)
}

# Import death data
data_death <- function(country, download, death_data_file, total, pred_date){
  # Import index
  index <- import_index()
  
  # Import deaths
  dt_death <- import_death(country, download, death_data_file, index, total)
  
  # Remove entries after pred_date
  dt_death <- dt_death[date <= pred_date,]
  dt_death <- dt_death[order(location_key),]
  
  # If any(colnames(dt_death) == "new_deceased"), then the death data is age-stratified, and
  # dt_death has to be turned into a long data set
  if(!any(colnames(dt_death) == "new_deceased")){
    dt_death <- melt(dt_death, measure.vars = patterns(new_deceased = "new_deceased", age_group = "age_bin"))
    reg <- unique(paste(dt_death$location_key, dt_death$age_group, sep = "."))
  } else{
    reg <- unique(dt_death$location_key)
    dt_death[, age_group := "Total"]
  }
  
  day_start <- as.numeric(as.Date(pred_date)) %% 7
  
  # Add "week" column to dt_death
  dt_death[, week := as.Date(date) + (day_start - as.numeric(as.Date(date))) %% 7]
  dt_death[is.na(new_deceased), new_deceased := 0]
  
  # Group number of new deaths by week
  dt_death <- dt_death[, lapply(.SD, sum), by = .(week, location_key, age_group), .SDcols = "new_deceased"]
  if(!all(get_nuts2_reg(reg, country) == reg)){
    dt_death <- dt_death[, lapply(.SD, sum), 
                         by = .(week, location_key = get_nuts2_reg(location_key, country), age_group), 
                         .SDcols = "new_deceased"]
    reg <- unique(paste(dt_death$location_key, dt_death$age_group, sep = "."))
  }
  dates <- unique(dt_death$week)
  # Generate matrix of death count per region and week
  mat_death <- matrix(dt_death$new_deceased, nrow = length(dates), ncol = length(reg), byrow = F)

  rownames(mat_death) <- as.character(dates)
  if(total) colnames(mat_death) <- unique(dt_death$location_key) else colnames(mat_death) <- reg
  return(mat_death)
}

# Implement regression models
regress_cfr <- function(dt_cfr, t_pred){
  # Create the data table used for the regression analysis
  dt_regress <- copy(dt_cfr)
  if(!any(colnames(dt_regress) == "age_group")){
    dt_regress[, age_group := "Total"]
  }
  
  # Compute weekly number of cases
  dt_regress[, cases_week := frollsum(cases, n = 7), by = .(nuts2)]
  # Remove entries posterior to the prediction date
  dt_regress <- dt_regress[date <= as.Date(t_pred),]
  # Use one date per week
  dates_regress <- as.character(
    seq(from = max(as.Date(dt_regress$date)), to = min(as.Date(dt_regress$date)), -7))
  dt_regress <- dt_regress[is.element(date, as.character(dates_regress)),]
  
  # Sort by date + region/age
  dt_regress <- dt_regress[order(nuts2, age_group, date),]
  
  # Compute delta CFR, and delta cases
  dt_regress[, delta_cfr1 := c(NA, diff(CFR))]
  dt_regress[, delta_cases1 := c(NA, diff(cases_week, lag = 1))]
  dt_regress[, delta_cfr2 := c(NA, NA, diff(CFR, lag = 2))]
  dt_regress[, delta_cases2 := c(NA, NA, diff(cases_week, lag = 2))]
  dt_regress[, delta_cfr3 := c(NA, NA, NA, diff(CFR, lag = 3))]
  dt_regress[, delta_cases3 := c(NA, NA, NA, diff(cases_week, lag = 3))]
  
  # Set first entries to NA
  dt_regress[date == min(dates_regress), delta_cfr1 := NA]
  dt_regress[date == min(dates_regress), delta_cases1 := NA]
  dt_regress[is.element(date, sort(dates_regress)[c(1,2)]), delta_cases2 := NA]
  dt_regress[is.element(date, sort(dates_regress)[c(1,2)]), delta_cfr2 := NA]
  dt_regress[is.element(date, sort(dates_regress)[c(1, 2, 3)]), delta_cases3 := NA]
  dt_regress[is.element(date, sort(dates_regress)[c(1, 2, 3)]), delta_cfr3 := NA]
  dt_regress[is.element(date, sort(dates_regress, decreasing = T)[c(1, 2, 3)]), delta_cfr1 := NA]
  dt_regress[is.element(date, sort(dates_regress, decreasing = T)[c(1, 2, 3)]), delta_cfr2 := NA]
  dt_regress[is.element(date, sort(dates_regress, decreasing = T)[c(1, 2, 3)]), delta_cfr3 := NA]
  
  ## Remove entries with low case number
  dt_regress[cases_week < 100 & date < (max(as.Date(dates_regress)) - 14), delta_cases1 := NA]
  dt_regress[cases_week < 100 & date < (max(as.Date(dates_regress)) - 14), delta_cases2 := NA]
  dt_regress[cases_week < 100 & date < (max(as.Date(dates_regress)) - 14), delta_cases3 := NA]
  
  for(age in unique(dt_regress$age_group)){
    # Run regression models
    lm1 <- lm(delta_cfr1 ~ delta_cases1, data = dt_regress[age_group == age,])
    lm2 <- lm(delta_cfr2 ~ delta_cases2, data = dt_regress[age_group == age,])
    lm3 <- lm(delta_cfr3 ~ delta_cases3, data = dt_regress[age_group == age,])

    ## Compute new CFR values
    # Predict the mean
    dt_regress[date == (max(as.Date(dates_regress)) - 14) & age_group == age, 
               delta_cfr1 := predict(lm1, newdata = data.frame(delta_cases1), se.fit = T)$fit]
    dt_regress[date == (max(as.Date(dates_regress)) - 7) & age_group == age, 
               delta_cfr2 := predict(lm2, newdata = data.frame(delta_cases2), se.fit = T)$fit]
    dt_regress[date == (max(as.Date(dates_regress))) & age_group == age, 
               delta_cfr3 := predict(lm3, newdata = data.frame(delta_cases3), se.fit = T)$fit]
    
    # Compute the standard deviation (in the prediction interval)
    dt_regress[date == (max(as.Date(dates_regress)) - 14) & age_group == age, 
               sd_cfr1 := sqrt((predict(lm1, newdata = data.frame(delta_cases1), se.fit = T)$se.fit)^2 + 
                                 (predict(lm1, newdata = data.frame(delta_cases1), se.fit = T)$residual.scale)^2)]
    dt_regress[date == (max(as.Date(dates_regress)) - 7) & age_group == age, 
               sd_cfr2 := sqrt((predict(lm2, newdata = data.frame(delta_cases2), se.fit = T)$se.fit)^2 + 
                                 (predict(lm2, newdata = data.frame(delta_cases2), se.fit = T)$residual.scale)^2)]
    dt_regress[date == (max(as.Date(dates_regress))) & age_group == age, 
               sd_cfr3 := sqrt((predict(lm3, newdata = data.frame(delta_cases3), se.fit = T)$se.fit)^2 + 
                                 (predict(lm3, newdata = data.frame(delta_cases3), se.fit = T)$residual.scale)^2)]
  }
  
  # Remove entries prior to three weeks before the prediction date (i.e. keep only the latest data point)
  dt_regress <- dt_regress[date >= (max(as.Date(dates_regress)) - 21),]
  
  # Move the predicted mean and standard deviation to columns delta_cfr and sd_cfr 
  dt_regress[!is.na(sd_cfr1), delta_cfr := delta_cfr1]
  dt_regress[!is.na(sd_cfr2), delta_cfr := delta_cfr2]
  dt_regress[!is.na(sd_cfr3), delta_cfr := delta_cfr3]
  dt_regress[!is.na(sd_cfr1), delta_cases := delta_cases1]
  dt_regress[!is.na(sd_cfr2), delta_cases := delta_cases2]
  dt_regress[!is.na(sd_cfr3), delta_cases := delta_cases3]
  dt_regress[!is.na(sd_cfr1), sd_cfr := sd_cfr1]
  dt_regress[!is.na(sd_cfr2), sd_cfr := sd_cfr2]
  dt_regress[!is.na(sd_cfr3), sd_cfr := sd_cfr3]
  
  # Set CFR as NA in entries where the CFR is predicted
  dt_regress[date >= (max(as.Date(dates_regress)) - 14), CFR := NA]
  
  # Return columns of interest
  dt_regress <- dt_regress[, .(date, CFR, age_group, nuts2, cases_week, delta_cases, delta_cfr, sd_cfr)]
  return(dt_regress)
}
