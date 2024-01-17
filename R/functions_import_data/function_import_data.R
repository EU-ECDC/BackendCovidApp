### function_import_data.R contains functions to import the different databases:
## - import_case: Import case data, from Google Covid-19 Open database or national databases
## - import_death: Import death data from national databases
## - import_eu_case: Import and compute the number of daily cases reported in Europe from ecdc case data
## - import_vacc: Import vaccine data, using ecdc vaccine database or national databases
## - import_variant: Import the daily proportion of sequences from each variant 
## - import_index: Import keys, codes and names for each region from Google covid-19 open data
## - import_pop: Import the number of inhabitant per country or nuts-3 regions, can be age-stratified
## - import_map: Import the shapefiles containing the map describing EU countries (nuts1, nuts2, and nuts3 level)
## - import_contact: Import number of contacts between age groups.
## - import_urban_rural: Import the urban-rural status from Eurostat database
## - import_test: Import number of tests per day, from national databases or ECDC test database
## - import_all_files: Import all data files, generate sts object and covariates

## Function to import all files
import_all_files <- function(country, total, range, delay_cumu, last_date, download_case = T, 
                             file_case = "", min_percent_variant = 10){
  ## Import index of Google place names
  index <- import_index()
  # Get country name from index
  country_name <- index[location_key == country, country_name]
  
  ## Import case data
  case_by_age <- import_case(country = country, download = download_case, 
                             data_file = file_case, index = index, total = total)
  if(max(case_by_age$date) < range[2] & last_date == FALSE){ 
    stop(paste0("In this country, the last cases were reported on ", max(case_by_age$date),
                ". The prediction date (set with 'pred_date') cannot be later than this."))
  } else if(last_date == TRUE) range[2] <- max(case_by_age$date)
  
  # Import variant, and population datasets
  variant <- import_variant(country, min_percent = min_percent_variant)
  
  ## TAG COUNTRY
  ## If pop data is taken from google database, set google <- T
  ## If vaccine data is taken from ECDC vaccine database, set vacc_ecdc <- T
  if (country == "IT") vacc_ecdc <- google <- T else vacc_ecdc <- google <- F
  pop_regions <- import_pop(country, index, total, google, national = F)
  pop_national <- import_pop(country, index, total, google, national = T)
  # Import cases in Europe per day 
  case_europe <- import_eu_case(country_exclude = country)
  
  # Move case_by_age to long format, and compute cumulative incidence by
  # age, region and variant
  cumu_incidence_regions <- cumu_incidence(country = country, 
                                           dt_variant = variant, 
                                           dt_case = case_by_age, 
                                           dt_pop = pop_regions, total = total)
  
  # Get age groups from case data
  age_groups <- unique(cumu_incidence_regions$age)
  
  ## TAG COUNTRY
  # Import test data if it's not already present in the case data, define test <- NULL
  # If testing data comes from age-stratified national database, use import_test (total == F)
  # If testing data comes from the ECDC database, run import_tests with total = T
  if (country == "FR"){ # country has age-stratified subnational testing data (e.g. France)
    # 'test' is not used so set to NULL
    test <- NULL
  } else if (country == "CZ"){ # country has age-stratified national testing data
    test <- import_test(country, pop_national, age_groups)
  } else if (country == "IT"){ # country does not have age-stratified national testing data
    # use overall national testing data
    test <- import_test(country, pop_national, age_groups)
  }
  
  ## Import vaccination data
  vacc_age_regions <- import_vacc(country, vacc_ecdc = vacc_ecdc, total)
  ## Adjust vaccine uptake data to match incidence age groups
  if (!total){
    vacc_age <- adjust_age_group(vacc_age_regions, pop_regions, age_groups, country = country,
                                 total = F, national = F)
    min_ages <- get_min_age(age_groups)
    nb_age <- length(min_ages)
    agg_pop <- pop_covariate(pop_regions, cumu_incidence_regions, min_ages, tot = T)
  } else {
    agg_pop <- pop_regions[,population]
    if(country == "CZ") names(agg_pop) <- paste0("CZ_", pop_regions$number) else names(agg_pop) <- pop_regions$number
    if(country == "FR") agg_pop <- agg_pop[unique(cumu_incidence_regions$region)] else agg_pop <- agg_pop[unique(cumu_incidence_regions$location_key)]
    vacc_age <- adjust_age_group(vacc_age_regions, pop_regions, total = T, 
                                 national = F, country = country, vacc_ecdc = vacc_ecdc)
  }
  # Import urban/rural status (NUTS3 only)
  urban <- import_urban_rural(country)
  
  ## Import EU NUTS region map
  # url <- "https://gisco-services.ec.europa.eu/distribution/v2/nuts/shp/NUTS_RG_20M_2021_3035.shp.zip"
  # map <- import_map(country, url, index, cumu_incidence_regions)
  map <- import_map(country, url = FALSE, index, cumu_incidence_regions)
  
  ## Generate population covariate matrix, and case data in hhh4 format
  case_by_age_arr <- nb_cases(cumu_incidence_regions)
  
  ## Import contact matrix
  if (total){
    pop_age <- NULL
    C <- 1
  } else { 
    pop_age <- pop_covariate(pop_regions, cumu_incidence_regions, min_ages, tot = F)
    C <- import_contact(country, index, age_groups, pop_age, pop_national)
  }  
  
  ## Create sts object
  if (total) agegroups <- NULL else agegroups <- rep(1,nb_age)
  if (total) sts_by <- "districts" else sts_by <- "all"
  sts_obj_nei <- create_sts(country, case_by_age_arr, agg_pop, map, C = C, by = sts_by,  
                            flatten = T, timeRange = range, agegroups = agegroups)
  
  ## Generate all covariate matrices
  covariates <- all_covariate(country = country, range = range, 
                              sts_object = sts_obj_nei, dt_pop = pop_regions, 
                              dt_incidence = cumu_incidence_regions, map = map, 
                              delay_inc = delay_cumu, dt_test = test, 
                              dt_vacc = vacc_age, total = total,
                              dt_europe = case_europe, dt_urban = urban,
                              dt_variant = variant)
  
  return(list(sts_obj_nei = sts_obj_nei, C = C, map = map,
              covariates = covariates, age_groups = age_groups))
}

#### Import case data 
import_case <- function(country, download, data_file, index, total = F){
  if (total & !(country %in% c("FR", "CZ"))){
    # Download and import case file, if you want to import a local version, set download to FALSE
    ## TAG COUNTRY:
    ## NON-AGE-STRATIFIED MODEL:
    ## Import the case data from a country-specific URL (containing the number of daily cases)
    ## The output of this if section should be a data frame named "case1".
    ## Process the data so that case1 contains the following variables:
    # location_key: NUTS-3 region code
    # date: date in format YYYY-MM-DD
    # new_confirmed: daily number of new confirmed cases
    # new_tested: daily number of new individuals tested; This should be NA if the testing data is missing from the case data source
    # cumulative_confirmed: daily cumulative number of confirmed cases
    # cumulative_tested: daily cumulative number of individuals tested; this should be NA if the testing data is missing from the case data source
    # week_day: year and week number in the format YYYY-WW where WW is the week of the year as a decimal number (00-53), with the first Monday of the year as day 1 of week 1.
    if (country == "IT"){
      if (download){
        data_file <- "https://github.com/pcm-dpc/COVID-19/raw/master/dati-province/dpc-covid19-ita-province.csv"
      }
      case <- fread(data_file, fill = TRUE)
      case1 <- case[, .(date = as.Date(data), name_reg = denominazione_regione,
                        name_pro = denominazione_provincia, nuts2 = codice_nuts_2, 
                        nuts3 = codice_nuts_3, cumulative_confirmed = totale_casi)]
      case1 <- case1[name_reg != "",]
      # Add nuts2 entry to all rows
      for(i in unique(case1$name_reg)){
        nuts2_i <- unique(case1[name_reg == i & nuts2 != "", nuts2])
        case1[name_reg == i, nuts2 := nuts2_i]
      }
      # Add nuts3 entry to all rows
      for(i in unique(case1$name_pro)){
        nuts3_i <- unique(case1[name_pro == i & nuts3 != "", nuts3])
        if (length(nuts3_i) > 0) case1[name_pro == i, nuts3 := nuts3_i]
      }
      # Define column "location_key"
      case1[, location_key := nuts3]
      # Remove entries without nuts3 and select date after 1st September 2020
      case1 <- case1[nuts3 != "" & date >= "2020-09-01",]
      
      
      case_ref <- copy(case1)
      case_ref[, ID := paste0(nuts3, "_", date)]
      # select entries where cumulative confirmed == 0, and set cumulative_confirmed
      # to the value at previous date
      setkey(case_ref, "ID")
      id_zeroes <- paste0(case1[cumulative_confirmed == 0,]$nuts3, 
                          "_", case1[cumulative_confirmed == 0,]$date - 1)
      case1[cumulative_confirmed == 0, 
            cumulative_confirmed := case_ref[id_zeroes, cumulative_confirmed]]
      ## For each nuts3 entry, fix entries with drops in cumulative_confirmed
      for(i in unique(case1$nuts3)){
        # Compute number of new_confirmed
        case1[nuts3 == i, new_confirmed := c(0, diff(cumulative_confirmed))]
        # Select dates where new_confirmed is below 0
        dates_zeroes <- case1[nuts3 == i & new_confirmed < 0,]$date
        ## Set the number of cumulative cases at dates_zeroes and dates_zeroes - 1 to 
        ## the value at dates_zeroes - 2
        case1[nuts3 == i & is.element(date, dates_zeroes - 1), 
              cumulative_confirmed := case1[nuts3 == i & is.element(date, dates_zeroes - 2), 
                                            cumulative_confirmed]]
        case1[nuts3 == i & is.element(date, dates_zeroes), 
              cumulative_confirmed := case1[nuts3 == i & is.element(date, dates_zeroes - 2), 
                                            cumulative_confirmed]]
        ## Re-compute new_confirmed
        case1[nuts3 == i, new_confirmed := c(0, diff(cumulative_confirmed))]
        ## When new_confirmed < 0, set it to 0
        case1[nuts3 == i & new_confirmed < 0, new_confirmed := 0]
      }
      ## Remove extreme values (in ITH10 and ITG17), where the number of cases in
      ## a day exceeded 1% of the population
      case1[nuts3 == "ITH10" & new_confirmed > 8000, new_confirmed := 0]
      case1[nuts3 == "ITG17" & new_confirmed > 20000, new_confirmed := 0]
      ## Same in ITH52, remove recent extreme value (one-day-spike)
      case1[nuts3 == "ITH52" & new_confirmed > 5000, new_confirmed := 0]
      
      ## Select columns location_key, date, new_confirmed, and cumulative_confirmed
      case1 <- case1[, .(location_key = nuts3, date = as.character(date), 
                         new_confirmed, new_tested = NA, 
                         cumulative_confirmed, cumulative_tested = NA)]
      case1[, week_day := format(as.Date(date), "%Y-%W")]
      
    }
    ## Select required columns
    cols <- c("date","location_key",names(case1)[grep("_confirmed|_tested",names(case1))])
    case1 <- case1[,..cols]
    regs <- case1[,unique(location_key)]
    
    nb_regs <- length(regs)
    # Fill new cases for any missing date for each region with 0s
    dates <- as.character(seq.Date(case1[,as.Date(min(date))], 
                                   case1[,as.Date(max(date))], by=1))
    base <- CJ(location_key = regs,date = dates)
    case1 <- merge(base, case1, by = c("location_key","date"), all.x = T)
    case1[is.na(new_confirmed), new_confirmed := 0]
    
    ## Add column with week of confirmation
    case1[, week_day := format(as.Date(date), "%Y-%W")]
    return(case1)
  } else {
    ## TAG COUNTRY:
    ## AGE-STRATIFIED MODEL:
    ## Import the case data from a country-specific URL (containing the number 
    ## of daily cases, per 10 year age band).
    ## The output of this if section should be a data table named "case_by_age1".
    ## Process the data so that case_by_age1 contains the following variables:
    # date: date in format YYYY-MM-DD
    # location_key: NUTS-3 region code
    # 10 columns "new_confirmed_age_[0 to 9]: daily number of new confirmed cases in each age group
    # 10 columns "cumulative_confirmed_age_[0 to 9]": daily cumulative number of confirmed cases in each age group
    # 10 columns "new_tested_age_[0 to 9]": daily number of new individuals tested per age group; This should be NA if the testing data is missing from the case data source
    # 10 columns "cumulative_tested_age_[0 to 9]": daily cumulative number of individuals tested; this should be NA if the testing data is missing from the case data source
    # 10 columns "age_bin_[0 to 9]": 1 value per column for all rows, corresponding to the age group in this bin (age_bin_0 = "0-9" ; age_bin_1 = "10-19" etc..).   
    # week_day: year and week number in the format YYYY-WW where WW is the week of the year as a decimal number (00-53), with the first Monday of the year as day 1 of week 1.
    if (country == "FR"){
      # Download and import case file, if you want to import a local version, set download to FALSE
      if (download){
        data_file <- "https://www.data.gouv.fr/fr/datasets/r/674bddab-6d61-4e59-b0bd-0be535490db0"
      }
      case_by_age <- fread(data_file)
      
      ## Rename columns
      setnames(case_by_age, c("jour","cage10","P","T"), c("date","age_bin","new_confirmed_age","new_tested_age"))
      cols <- c("pop","new_confirmed_age","new_tested_age")
      case_by_age[, (cols) := lapply(.SD, function(x) as.integer(sub(",",".",x))), .SDcols = cols]
      setorder(case_by_age, dep, age_bin, date)
      case_by_age[, `:=`(cumulative_confirmed_age = cumsum(new_confirmed_age),
                         cumulative_tested_age = cumsum(new_tested_age)), by = .(dep, age_bin)]
      ## Change age_bin
      case_by_age[,age_bin := paste0(as.numeric(sub("\\[([0-9]+)-.*","\\1",age_bin)),"-",
                                     as.numeric(sub(".*-(.*))","\\1",age_bin)) - 1)]
      case_by_age[age_bin == "90-Inf", age_bin := "90-"]
      case_by_age[, grp := .GRP - 1, by = .(age_bin)]
      case_by_age <- dcast(case_by_age, dep + date ~ grp, value.var = c("new_confirmed_age","cumulative_confirmed_age","new_tested_age","cumulative_tested_age","age_bin"))
      case_by_age <- merge(index[country_code == "FR", .(location_key, subregion2_code)], case_by_age, by.x = "subregion2_code", by.y = "dep")
      case_by_age[, subregion2_code := NULL]
      case_by_age[, date := as.character(date)]
      ## Select subregion2 level
      case_by_age1 <- case_by_age[substr(location_key, 1, 2) == country & 
                                    nchar(location_key) > 6,]
      # Replace NA entry by 0
      case_by_age1[location_key == "FR_OCC_48" & 
                     date %in% c("2023-05-09", "2023-06-07", "2023-06-08", "2023-06-09", 
                                 "2023-06-10", "2023-06-11", "2023-06-19", "2023-06-20", 
                                 "2023-06-21", "2023-06-22", "2023-06-23", "2023-06-24"), 
                   new_tested_age_0 := 0]
      case_by_age1[location_key == "FR_NAQ_23" & 
                     date %in% c("2023-05-19", "2023-05-20", "2023-05-21", "2023-05-22"), 
                   new_tested_age_0 := 0]
      case_by_age1[location_key == "FR_NAQ_23" & date %in% c("2023-06-21"), 
                   new_tested_age_1 := 0]
      case_by_age1[location_key == "FR_OCC_48" & date %in% c("2023-06-19", "2023-06-20"), 
                   new_tested_age_2 := 0]
      case_by_age1[location_key == "FR_OCC_48" & date %in% c("2023-06-26", "2023-06-27"), 
                   new_tested_age_8 := 0]
      case_by_age1[location_key == "FR_OCC_48" & 
                     date %in% c("2023-06-02", "2023-06-03", "2023-06-04", "2023-06-05",
                                 "2023-06-06"),
                   new_tested_age_9 := 0]
      ## Remove Corsica
      case_by_age1 <- 
        case_by_age1[- grep("_COR_", case_by_age1$location_key),]      
    } else if (country == "CZ"){
      # Download and import case file, if you want to import a local version, set download to FALSE
      if (download){
        data_file <- "https://onemocneni-aktualne.mzcr.cz/api/v2/covid-19/osoby.csv"
      }
      case_by_age <- fread(data_file)
      colnames(case_by_age) <- c("id", "date", "age", "gender", "nuts", "lau", 
                                 "infection_abroad", "country_infection", "reported_by_khs")
      # Remove entries where age and nuts is not reported
      case_by_age <- case_by_age[nuts != "" & !is.na(age), .(date, age, nuts)]
      # Set age group
      case_by_age[, age_group := floor(age / 10)]
      case_by_age[, n := 1]
      # Compute the number of cases per day, nuts, and age group
      case_by_age <- case_by_age[, lapply(.SD, sum), by = .(date, nuts, age_group), .SDcols = "n"]
      case_by_age[, id := paste(date, nuts, age_group,sep = "_")]
      setkey(case_by_age, id)
      
      ## Create output matrix case_by_age1
      # Generate a vector containing all dates
      min_date <- min(case_by_age$date)
      max_date <- max(case_by_age$date)
      all_dates <- seq(min_date, max_date, "day")
      # Extract all nuts regions
      all_nuts <- unique(case_by_age$nuts)
      # Create case_by_age1, containing the number of new and cumulative cases, and 
      # new and cumulative tests (set to NA in Czechia)
      case_by_age1 <- as.data.frame(
        matrix(NA, ncol = 52, nrow = length(all_dates) * length(all_nuts)))
      colnames(case_by_age1) <- 
        c("date", "location_key", paste0("new_confirmed_age_", 0:9),
          paste0("cumulative_confirmed_age_", 0:9), paste0("new_tested_age_", 0:9), 
          paste0("cumulative_tested_age_", 0:9), paste0("age_bin_", 0:9))
      # Set the dates / nuts / age_bin columns
      case_by_age1$date <- as.character(rep(all_dates, each = length(all_nuts)))
      case_by_age1$location_key <- rep(all_nuts, length(all_dates))
      case_by_age1[, grep("age_bin", colnames(case_by_age1))] <- 
        matrix(c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", 
                 "70-79", "80-89", "90+"), nrow = nrow(case_by_age1), ncol = 10, byrow = T)
      # Set the new_confirmed columns to 0
      case_by_age1[, grep("new_confirmed_age_", colnames(case_by_age1))] <- 0
      
      # In each age group / date / region, add the number of reported cases from case_by_age
      for(i in seq_len(10) - 1){
        col_i <- paste0("new_confirmed_age_", i)
        case_by_age1[, col_i] <- case_by_age1[, col_i] + 
          case_by_age[paste(case_by_age1$date, case_by_age1$location_key, i, sep = "_"), n]
        case_by_age1[is.na(case_by_age1[, col_i]), col_i] <- 0
      }
      case_by_age1 <- as.data.table(case_by_age1)
      # Change location_key code
      case_by_age1[, location_key := paste0("CZ_", substr(location_key, 4, 5))]
    } 
    ## Select required columns
    cols <- c("date","location_key",names(case_by_age1)[grep("_confirmed|_tested|age_bin_",names(case_by_age1))])
    case_by_age1 <- case_by_age1[,..cols]
    regs <- case_by_age1[,unique(location_key)]
    nb_regs <- length(regs)
    nms <- names(case_by_age1)
    
    # Fill cases for any missing date for each region with 0s
    dates <- as.character(seq.Date(case_by_age1[,as.Date(min(date))],
                                   case_by_age1[,as.Date(max(date))], by=1))
    base <- CJ(location_key = regs,date = dates)
    case_by_age1 <- merge(base, case_by_age1, by = c("location_key","date"), all.x = T)
    cols <- nms[grep("age_bin_", nms)]
    age_labels <- case_by_age1[!is.na(age_bin_0), ..cols][1,]
    case_by_age1[is.na(age_bin_0), (cols) := age_labels]
    setnafill(case_by_age1, fill = 0,
              cols = nms[grep("new_confirmed_",nms)])
    
    ## Add column with week of confirmation
    case_by_age1[, week_day := format(as.Date(date), "%Y-%W")]
    
    # Combine data for 80-89 and 90+ age groups
    case_by_age1[, age_bin_8 := "80+"]
    case_by_age1[, age_bin_9 := NULL]
    case_by_age1[, cumulative_tested_age_8 := cumulative_tested_age_8 + cumulative_tested_age_9]
    case_by_age1[, cumulative_tested_age_9 := NULL]
    case_by_age1[, new_tested_age_8 := new_tested_age_8 + new_tested_age_9]
    case_by_age1[, new_tested_age_9 := NULL]
    case_by_age1[, cumulative_confirmed_age_8 := cumulative_confirmed_age_8 + cumulative_confirmed_age_9]
    case_by_age1[, cumulative_confirmed_age_9 := NULL]
    case_by_age1[, new_confirmed_age_8 := new_confirmed_age_8 + new_confirmed_age_9]
    case_by_age1[, new_confirmed_age_9 := NULL]
    
    if(total){
      case1 <- case_by_age1[, .(location_key, date, 
                                new_confirmed = new_confirmed_age_0 + new_confirmed_age_1 + new_confirmed_age_2 + new_confirmed_age_3 + 
                                  new_confirmed_age_4 + new_confirmed_age_5 + new_confirmed_age_6 + new_confirmed_age_7 + new_confirmed_age_8, 
                                new_tested = new_tested_age_0 + new_tested_age_1 + new_tested_age_2 + new_tested_age_3 + 
                                  new_tested_age_4 + new_tested_age_5 + new_tested_age_6 + new_tested_age_7 + new_tested_age_8, 
                                cumulative_confirmed = cumulative_confirmed_age_0 + cumulative_confirmed_age_1 + cumulative_confirmed_age_2 + cumulative_confirmed_age_3 + 
                                  cumulative_confirmed_age_4 + cumulative_confirmed_age_5 + cumulative_confirmed_age_6 + cumulative_confirmed_age_7 + cumulative_confirmed_age_8, 
                                cumulative_tested = cumulative_tested_age_0 + cumulative_tested_age_1 + cumulative_tested_age_2 + cumulative_tested_age_3 + 
                                  cumulative_tested_age_4 + cumulative_tested_age_5 + cumulative_tested_age_6 + cumulative_tested_age_7 + cumulative_tested_age_8, 
                                week_day)]
      return(case1)
    }
    return(case_by_age1)    
  }
}

#### Import death data 
import_death <- function(country, download, data_file, index, total = F){
  if (total & !(country %in% c("FR", "CZ"))){
    # Download and import death file, if you want to import a local version, set download to FALSE
    ## TAG COUNTRY:
    ## NON-AGE-STRATIFIED MODEL:
    ## Import the death data from a country-specific URL (containing the number of daily deaths)
    ## The output of this if section should be a data frame named "death1".
    ## Process the data so that case1 contains the following variables:
    # location_key: region code (NUTS-2/NUTS-1 depending on data availability)
    # date: date in format YYYY-MM-DD
    # new_deceased: daily number of new deaths
    if (country == "IT"){
      if (download){
        data_file <- "https://github.com/pcm-dpc/COVID-19/raw/master/dati-regioni/dpc-covid19-ita-regioni.csv"
      }
      death <- as.data.table(read.csv2(data_file, sep = ","))
      regs <- unique(death[codice_nuts_2 != "", codice_nuts_2])
      names(regs) <- unique(death[codice_nuts_2 != "", codice_regione])
      death[, codice_nuts_2 := regs[as.character(codice_regione)]]
      
      death[, date := as.character(as.Date(data))]
      
      death1 <- death[, .(date, codice_nuts_2, deceduti)]
      colnames(death1) <- c("date", "location_key", "cumulative_death")
      for(i in unique(death1$location_key)){
        # Compute number of deaths per day
        death1[location_key == i, new_deceased := c(0, diff(cumulative_death))]
      }
      ## Overwrite negative deaths with 0s 
      death1[new_deceased < 0, new_deceased := 0]
    }
    
    return(death1)
  } else {
    ## TAG COUNTRY:
    ## AGE-STRATIFIED MODEL:
    ## Import the death data from a country-specific URL (containing the number 
    ## of daily deaths, per 10 year age band).
    ## The output of this if section should be a data frame named "death_by_age1".
    ## Process the data so that death_by_age1 contains the following variables:
    # date: date in format YYYY-MM-DD
    # location_key: region code (NUTS-2/NUTS-1 depending on data availability)
    # 10 columns "new_deceased_age_[0 to 9]: daily number of new confirmed deaths in each age group
    # 10 columns "age_bin_[0 to 9]": 1 value per column for all rows, corresponding to the age group in this bin (age_bin_0 = "0-9" ; age_bin_1 = "10-19" etc..).   
    if (country == "FR"){
      if (download){
        data_file <- "https://www.data.gouv.fr/fr/datasets/r/08c18e08-6780-452d-9b8c-ae244ad529b3"
      }
      death_by_age <- fread(data_file)
      setnames(death_by_age, c("cl_age90","jour","dc"), c("age_bin","date","cumulative_deceased"))
      death_by_age[, date := as.character(date)]
      ## Calculate new deaths
      death_by_age[, new_deceased_age := diff(c(0, cumulative_deceased)), by = .(reg, age_bin)]
      death_by_age[, cumulative_deceased := NULL]
      ## Overwrite negative deaths with 0s
      death_by_age[new_deceased_age < 0, new_deceased_age := 0]
      ## Remove total row
      death_by_age <- death_by_age[age_bin != 0]
      ## Change age_bin
      death_by_age[,age_bin := fcase(age_bin < 90, paste0(age_bin - 9, "-", age_bin),
                                     age_bin == 90, "90-")]
      death_by_age[, grp := .GRP - 1, by = .(age_bin)]
      death_by_age <- dcast(death_by_age, reg + date ~ grp, value.var = c("new_deceased_age","age_bin"))
      ## Change INSEE region names to Google region names
      reg_nbs <- fread("Data/reg2016.txt", encoding = "Latin-1")
      reg_nbs[, reg_nm := NCCENR]
      reg_nbs[reg_nm == "Guyane", reg_nm := "French Guiana"]
      reg_nbs[reg_nm == "Normandie", reg_nm := "Normandy"]
      reg_nbs[reg_nm == "Nord-Pas-de-Calais-Picardie", reg_nm := "Hauts-de-France"]
      reg_nbs[reg_nm == "Alsace-Champagne-Ardenne-Lorraine", reg_nm := "Grand Est"]
      reg_nbs[reg_nm == "Bretagne", reg_nm := "Brittany"]
      reg_nbs[reg_nm == "Aquitaine-Limousin-Poitou-Charentes", reg_nm := "Nouvelle-Aquitaine"]
      reg_nbs[reg_nm == "Languedoc-Roussillon-Midi-Pyrénées", reg_nm := "Occitanie"]
      reg_nbs[reg_nm == "Corse", reg_nm := "Corsica"]
      
      ## Add Google location key for regions 
      reg_nbs[, location_key := index[match(reg_nbs[,reg_nm], subregion1_name), location_key]]
      death_by_age[, location_key := reg_nbs[match(death_by_age[,reg],REGION), location_key]]
      
      ## Remove overseas territories
      death_by_age <- death_by_age[!(reg %in% c(1,2,3,4,6))]
      ## Remove Corsica
      death_by_age1 <- death_by_age[reg != 94]
    } else if (country == "CZ"){
      if (download){
        data_file <- "https://onemocneni-aktualne.mzcr.cz/api/v2/covid-19/umrti.csv"
      }
      death_by_age <- fread(data_file)
      colnames(death_by_age) <- c("id", "date", "age", "gender", "nuts", "lau")
      # Remove entries where age and nuts is not reported
      death_by_age <- death_by_age[nuts != "" & !is.na(age), .(date, age, nuts)]
      # Set age group
      death_by_age[, age_group := floor(age / 10)]
      death_by_age[, n := 1]
      # Compute the number of cases per day, nuts, and age group
      death_by_age <- death_by_age[, lapply(.SD, sum), by = .(date, nuts, age_group), .SDcols = "n"]
      death_by_age[, id := paste(date, nuts, age_group,sep = "_")]
      setkey(death_by_age, id)
      
      ## Create output matrix death_by_age1
      # Generate a vector containing all dates
      min_date <- min(death_by_age$date)
      max_date <- max(death_by_age$date)
      all_dates <- seq(min_date, max_date, "day")
      # Extract all nuts regions
      all_nuts <- unique(death_by_age$nuts)
      # Create death_by_age1, containing the number of new and cumulative cases, and 
      # new and cumulative tests (set to NA in Czechia)
      death_by_age1 <- as.data.frame(
        matrix(NA, ncol = 22, nrow = length(all_dates) * length(all_nuts)))
      colnames(death_by_age1) <- 
        c("date", "location_key", paste0("new_deceased_age_", 0:9), paste0("age_bin_", 0:9))
      # Set the dates / nuts / age_bin columns
      death_by_age1$date <- as.character(rep(all_dates, each = length(all_nuts)))
      death_by_age1$location_key <- rep(all_nuts, length(all_dates))
      death_by_age1[, grep("age_bin", colnames(death_by_age1))] <- 
        matrix(c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", 
                 "70-79", "80-89", "90+"), nrow = nrow(death_by_age1), ncol = 10, byrow = T)
      # Set the new_confirmed columns to 0
      death_by_age1[, grep("new_deceased_age_", colnames(death_by_age1))] <- 0
      
      # In each age group / date / region, add the number of reported cases from death_by_age
      for(i in seq_len(10) - 1){
        col_i <- paste0("new_deceased_age_", i)
        death_by_age1[, col_i] <- death_by_age1[, col_i] + 
          death_by_age[paste(death_by_age1$date, death_by_age1$location_key, i, sep = "_"), n]
        death_by_age1[is.na(death_by_age1[, col_i]), col_i] <- 0
      }
      death_by_age1 <- as.data.table(death_by_age1)
      # Change location_key code
      death_by_age1[, location_key := paste0("CZ_", substr(location_key, 4, 5))]
    }
    ## Select required columns
    cols <- c("date","location_key",names(death_by_age1)[grep("new_deceased|age_bin_",names(death_by_age1))])
    death_by_age1 <- death_by_age1[,..cols]
    regs <- death_by_age1[,unique(location_key)]
    nb_regs <- length(regs)
    nms <- names(death_by_age1)
    
    
    # Fill cases for any missing date for each region with 0s
    dates <- as.character(seq.Date(death_by_age1[,as.Date(min(date))],
                                   death_by_age1[,as.Date(max(date))], by=1))
    base <- CJ(location_key = regs,date = dates)
    death_by_age1 <- merge(base, death_by_age1, by = c("location_key","date"), all.x = T)
    cols <- nms[grep("age_bin_", nms)]
    age_labels <- death_by_age1[!is.na(age_bin_0), ..cols][1,]
    death_by_age1[is.na(age_bin_0), (cols) := age_labels]
    setnafill(death_by_age1, fill = 0,
              cols = nms[grep("new_deceased_",nms)])
    
    # Combine data for 80-89 and 90+ age groups
    death_by_age1[, age_bin_8 := "80+"]
    death_by_age1[, age_bin_9 := NULL]
    death_by_age1[, new_deceased_age_8 := new_deceased_age_8 + new_deceased_age_9]
    death_by_age1[, new_deceased_age_9 := NULL]
    if(total){
      death1 <- death_by_age1[, .(location_key, date, 
                                  new_deceased = new_deceased_age_0 + new_deceased_age_1 + new_deceased_age_2 + new_deceased_age_3 + 
                                    new_deceased_age_4 + new_deceased_age_5 + new_deceased_age_6 + new_deceased_age_7 + new_deceased_age_8)]
      return(death1)
    }
  }
  return(death_by_age1)
}

#### Import all cases from Europe:
import_eu_case <- function(country_exclude){
  # import WHO case data
  data <- read.csv("https://covid19.who.int/WHO-COVID-19-global-data.csv", na.strings = "")
  # Remove country of interest
  dt_eu <- as.data.table(data[data$WHO_region == "EURO" & data$Country_code != country_exclude, c("Country_code", "Date_reported", "New_cases")])
  colnames(dt_eu) <- c("country", "date", "cases")
  # Change date format
  dt_eu[,date := as.Date(date)]
  # Compute number of daily cases reported in Europe
  dt_eu_sum <- dt_eu[, lapply(.SD, sum), by = date, .SDcols = "cases"]
  dt_eu_sum <- dt_eu_sum[!is.na(cases),]
  return(dt_eu_sum)
}

#### Import vaccine data
### If national == T: Overall population; if national == F: By NUTS-3 region 
### If total == T: Merge all age groups; if total == F: Stratified by age
### If vacc_ecdc == T: Use vaccination ECDC data; if vacc_ecdc == F: Use national databases
import_vacc <- function(country, vacc_ecdc, total = F, national = F){
  temp <- tempfile()
  if (vacc_ecdc){
    # Download ECDC vaccine data
    download.file("https://opendata.ecdc.europa.eu/covid19/vaccine_tracker/csv/data.csv",
                  temp)
    vacc_eu <- as.data.table(read.csv2(temp, sep = ","))
    unlink(temp)
    # Select entries from regions of country of interest
    vacc <- vacc_eu[ReportingCountry == country,]
    # Remove the column "Vaccine", which is always "UNK"
    # Remove the column "unknown dose", which is always 0
    vacc[, Vaccine := NULL]
    vacc[, FirstDoseRefused := NULL]
    vacc[, UnknownDose := NULL]
    vacc[, Population := NULL]
    vacc[, NumberDosesExported := NULL]
    
    nb_dose <- sum(grepl("DoseAdditional", colnames(vacc))) + 2
    colnames(vacc) <- c("date_week", "Country", "Population", "Numberdoses", 
                        paste0("dose", seq_len(nb_dose)), "NUTS", "Age")
    for(i in seq_len(nb_dose)[-c(1, 2, 3)]){
      vacc[, dose3 := dose3 + eval(parse(text = paste0("dose", i)))]
      vacc[, c(paste0("dose", i)) := NULL]
    }

    # Convert date_week variable to date
    vacc[,date_week := as.character(ISOweek2date(paste0(date_week,"-7")))]
    
    # Sum doses over different vaccines
    cols <- c("dose1", "dose2", "dose3")
    vacc <- vacc[, lapply(.SD,function(x) sum(x, na.rm = T)), .SDcols = cols, 
                 by = .(date_week, Country, Population, NUTS, Age)]
    
    # Calculate proportions vaccinated with different doses
    vacc[, dose1_prop := dose1 / Population]
    vacc[, dose2_prop := dose2 / Population]
    vacc[, dose3_prop := dose3 / Population]

    # Clean age group names
    vacc[, Age := sub("_","-",sub("Age","", Age))]
    vacc[Age == "<18", Age := "0-17"]
    
    # Drop doses for HCW and LTCF 
    vacc <- vacc[!(Age %in% c("HCW","LTCF"))]
    
    # Drop doses split by age group below 18 as not all countries have these
    # Can look at separating countries with and without these later
    vacc <- vacc[!(Age %in% c("0-4","5-9","10-14","15-17"))]
    ## TAG COUNTRY:
    ## If any entry has to be removed (e.g. overseas territory), remove the rows here
    ## If national, exclude doses coming from removed regions. If no region should be 
    ## removed, there's no need to change anything.
    if (country == "FR"){
      # Remove overseas territories
      vacc <- vacc[!is.element(NUTS, c("FRM","FRY1", "FRY2", "FRY3", "FRY4", "FRY5")),]
      if (national){
        vacc_regions <- vacc[NUTS != "FR", lapply(.SD,function(x) sum(x, na.rm = T)), .SD = cols, by = .(date_week, Country, Age)]
        vacc_regions[, NUTS := "FR"]
      } else{
        vacc_regions <- vacc[NUTS != country]  
      }
    } else {
      if (national){
        vacc_regions <- vacc[NUTS == country]  
      } else{
        vacc_regions <- vacc[NUTS != country]  
      }
    }
    if(total){
      # Extract total (non-age-stratified) vaccinations (child + adult). 
      # In Italy and Czechia, ALL does not include children
      if(country %in% c("IT", "CZ")){
        vacc_regions <- vacc_regions[Age %in% c("0-17","ALL")]
      } else vacc_regions <- vacc_regions[Age == "ALL"]
      vacc_regions[, Age := "Total"]
      vacc_regions <- vacc_regions[, lapply(.SD, sum), .SDcols = cols, by = .(date_week, NUTS, Age)]
      colnames(vacc_regions) <- c("date_week", "number", "Age", "dose1", "dose2", "dose3")
    }
    ## In ECDC data, the geographic resolution is different in vacc (NUTS 2) 
    ## and case_by_age1 (NUTS 3)
  } else {
    ## TAG COUNTRY:
    ## ALL MODEL WITH NATIONAL DATABASES:
    ## Import the vaccine data from a country-specific URL (containing the number 
    ## of daily cases, per 10 year age band).
    ## The output of this if() section should be a data frame named "vacc_regions".
    ## IF TOTAL == T: Process the data so that vacc_regions contains the following variables:
    # date_week: date in format YYYY-MM-DD
    # number: Region code
    # Population: Number of inhabitants in this region + age group (If not reported, can be set to NA, for now it's only used in France)
    # Age: Age group described in this row
    # dose1: Number of 1st dose distributed on that day / region / age group
    # dose2: Number of 2nd dose distributed on that day / region / age group
    # dose3: Number of 3rd dose distributed on that day / region / age group
    if (country == "FR"){
      # Import publicly available French data on coverage at NUTS3 level
      download.file("https://datavaccin-covid.ameli.fr/explore/dataset/donnees-vaccination-par-tranche-dage-type-de-vaccin-et-departement/download/?format=csv&timezone=Europe/London&lang=fr&use_labels_for_header=true&csv_separator=%3B",
                    temp)
      vacc <- as.data.table(read.csv2(temp, sep = ";"))
      # Weeks 2020-53 and 2021-01 are the same
      vacc <- vacc[semaine_injection != "2020-53",]
      # Remove NUTS2 regions
      vacc <- vacc[departement_residence != "Tout département" | 
                     Region_residence == "Toute région",]
      # Select column of interest
      vacc <- vacc[, .(date, departement_residence,
                       population_insee, classe_age, type_vaccin,
                       effectif_1_inj, effectif_termine, Effectif_rappel)]
      # Select entries grouping all vaccines
      vacc <- vacc[type_vaccin == "Tout vaccin",]
      vacc <- vacc[, type_vaccin := NULL]
      # Set NA values to 0
      vacc[is.na(effectif_1_inj), effectif_1_inj := 0]
      vacc[is.na(effectif_termine), effectif_termine := 0]
      vacc[is.na(Effectif_rappel), Effectif_rappel := 0]
      
      # Remove overseas territories
      vacc <- vacc[!is.element(departement_residence, c(971, 972, 973, 974, 976, 999)),]
      # Remove Corsica:
      vacc <- vacc[!is.element(departement_residence, c("2A", "2B")),]
      # Rename columns
      colnames(vacc) <- c("date_week", "number", "Population", "Age",
                          "dose1", "dose2", "dose3")
      vacc[Age == "75 et +", Age := "75+"]
      vacc[Age == "TOUT_AGE", Age := "Total"]
      
      # Split in two datasets: vaccine per region, and national vaccine
      vacc_regions <- vacc[number != "Tout département",]
      
      vacc_regions[, ID := paste0(date_week, "-", number, "-", Age)]
      setkey(vacc_regions, ID)
      
      if (total){
        # Extract total regional vaccinations
        vacc_regions <- vacc_regions[Age == "Total"]  
      }
      unlink(temp)
    } else if (country == "CZ"){
      # Import publicly available Czech data on coverage at NUTS3 level
      download.file("https://onemocneni-aktualne.mzcr.cz/api/v2/covid-19/ockovani.csv", temp)
      vacc <- as.data.table(read.csv2(temp, sep = ","))
      setnames(vacc,
               c("datum","vakcina","kraj_nuts_kod","kraj_nazev","vekova_skupina",
                 "prvnich_davek","druhych_davek","celkem_davek"),
               c("date_week","vaccine","number","region","age","dose1","dose2","total"))
      
      # Calculate booster doses
      vacc[, dose3 := total - (dose1 + dose2)]
      
      # Change number from NUTS code to just region number
      vacc[, number := substr(number, 4, 5)]
      
      # Sum over vaccine types
      cols <- c("dose1","dose2","dose3","total")
      vacc <- vacc[, lapply(.SD,sum), .SDcols = cols, by = .(date_week, number, region, age)]
      
      # Make data table of vaccinations for all dates, regions and age groups
      dates <- vacc[,unique(date_week)]
      reg_nbs <- vacc[,unique(number)]
      age_groups_vax <- vacc[,sort(unique(age))]
      base <- CJ(date_week = dates, number = reg_nbs, age = age_groups_vax)
      vacc_regions <- merge(base, vacc[, !"region"], by = c("date_week", "number", "age"), all.x = T)
      vacc_regions[, region := vacc[match(vacc_regions[,number], number), region]]
      
      # Fill missing values with 0s
      setnafill(vacc_regions, fill = 0, cols = cols)
      
      # Add NA population column
      vacc_regions[, Population := NA]
      
      vacc_regions[, ID := paste0(date_week, "-", number, "-", age)]
      setkey(vacc_regions, ID)
      
      if (total){
        # Sum vaccinations for each region over age groups
        vacc_regions <- vacc_regions[, lapply(.SD, sum), .SDcols = cols, by = .(date_week, number, region)]
      }
      unlink(temp)
    }

  }
  
  return(vacc_regions)
}

#### Proportion of variants 
import_variant <- function(cntry, min_percent = 10){
  temp <- tempfile()
  # Download ECDC variant dataset
  download.file("https://opendata.ecdc.europa.eu/covid19/virusvariant/csv/data.csv",
                temp)
  variant_eu <- as.data.table(read.csv2(temp, sep = ","))
  unlink(temp)
  # Check if variant data contains Alpha and Delta, and use static version of dataset if not
  if (!all(c("B.1.1.7", "B.1.617.2") %in% variant_eu[, unique(variant)])){
    variant_eu <- as.data.table(read.csv2("Data/variant.csv", sep = ","))
  }
  

  # Select entries in country of interest
  variant <- variant_eu[country_code == cntry,]
  # Check if variant data contains number_sequenced_known_variant, and use number_sequenced if not
  if(any(colnames(variant) == "number_sequenced_known_variant")){
    variant[is.na(number_sequenced_known_variant), number_sequenced_known_variant := number_sequenced]
    variant <- 
      variant[, .(country_code, year_week, variant, source, percent_variant,
                  number_detections_variant, number_sequenced_known_variant)]
  } else
    variant <- 
    variant[, .(country_code, year_week, variant, source, percent_variant,
                number_detections_variant, number_sequenced)]
  
  colnames(variant) <- c("country_code", "year_week", "variant", "source", 
                         "percent_variant", "nb_detection", "nb_sequenced")
  colnames(variant)[1] <- "country"
  variant[percent_variant == "", percent_variant := 0]
  variant[, percent_variant := as.numeric(percent_variant)]
  # Find maximum value of each variant across the time frame. Remove variants
  # where percent_variant was never above 10% 
  max_prop_variant <- variant[!is.na(percent_variant),lapply(.SD, max), 
                              by = variant, .SDcols = "percent_variant"]
  # Remove variants that were never reported in more than "min_percent"% of the sequenced cases
  variant <- variant[
    is.element(variant, max_prop_variant[percent_variant > min_percent, variant])]
  # Select and rename columns of interest
  variant <- variant[variant != "Other",]
  # Rename variants
  variant[variant %in% c("B.1.1.529","BA.1","BA.2", "BA.4", "BA.5"), variant := "Omicron"]
  variant[variant == "B.1.617.2", variant := "Delta"]
  variant[variant == "B.1.1.7", variant := "Alpha"]
  variant[variant == "B.1.351", variant := "Beta"]
  variant[variant == "P.1", variant := "Gamma"]
  cols <- c("percent_variant","nb_detection")
  variant <- variant[,lapply(.SD,sum),.SDcols = cols, 
                     by = .(country,year_week,variant,source,nb_sequenced)]
  # variant[, percent_variant := nb_detection/nb_sequenced]
  # Merge TESSy and GISAID data
  dates_two_sources <- unique(variant[source == "TESSy", year_week])
  variants_two_sources <- unique(variant[source == "TESSy", variant])
  for(var in variants_two_sources){
    tessy_var <- 
      variant[is.element(year_week, dates_two_sources) & source == "TESSy" & 
                variant == var, .(year_week, nb_detection, nb_sequenced)]
    
    variant[is.element(year_week, tessy_var$year_week) & source == "GISAID" & 
              variant == var, nb_detection := nb_detection + tessy_var$nb_detection]
    variant[is.element(year_week, tessy_var$year_week) & source == "GISAID" & 
              variant == var, nb_sequenced := nb_sequenced + tessy_var$nb_sequenced]
    if (all(variant[variant == var, source] == "TESSy"))
      variant[variant == var, source := "GISAID"]
  }
  variant <- variant[source == "GISAID", ]
  variant[, source := NULL]
  # Recompute percent_variant after merging
  variant[, percent_variant := nb_detection / nb_sequenced]
  variant[nb_sequenced == 0, percent_variant := 0]
  variant[, date := ISOweek2date(paste0(sub("-","-W",year_week),"-1"))]
  variant[variant == "Omicron" & date < "2021-08-01", percent_variant := 0]
  variant[variant == "Delta" & date < "2021-01-01", percent_variant := 0]
  variant <- variant[order(date),]
  return(variant)
}

#### Import keys, codes and names for each region from Google covid-19 open data
import_index <- function(){
  temp <- tempfile()
  download.file("https://storage.googleapis.com/covid19-open-data/v3/index.csv",temp)
  index <- as.data.table(read.csv2(temp, sep = ",", encoding = "UTF-8"))
  unlink(temp)
  return(index)
}

#### Number of inhabitants: 
### If national == T: Overall population; if national == F: By NUTS-3 region 
### If total == T: Merge all age groups; if total == F: Stratified by age
### If google == T: Use Google covid-19 open data; if google == F: Use Eurostat
import_pop <- function(country, index, total = F, national = F, google = F){
  if (national){ # use Eurostat national population data by year of age
    # Read in population data at NUTS2 level
    pop <- as.data.table(read.csv2("Data/demo_r_d2jan_1_Data.csv", sep = ",", encoding = "latin1"))
    
    names(pop) <- tolower(names(pop))
    
    # Drop unneeded columns
    pop <- pop[, .(number = geo, region = geo_label, age, population = value)]
    
    pop[, population := suppressWarnings(as.numeric(gsub(",","",population)))]
    
    # Drop population with unknown age as it's very small
    pop <- pop[age != "UNK"]
    pop[, age := sub("_","",sub("Y","",age))]
    pop[age == "LT1", age := "0"]
    pop[age == "OPEN", age := "100"]
    ## TAG COUNTRY
    ## If any entry has to be removed (e.g. overseas territory), remove the rows 
    ## here (e.g. for France, exclude overseas territories). If no region should be 
    ## removed, there's no need to change anything.
    if (country == "FR"){
      pop <- pop[substr(number,1,2) == "FR" & nchar(number)==4]
      pop <- pop[!(number %in% c("FRM0","FRY1","FRY2","FRY3","FRY4","FRY5"))]
      ## Compute the number of inhabitants in the country, excluding overseas territory
      pop <- pop[,.(population = sum(population)), by = .(age)]
      pop[, `:=`(number = "FR", region = "France")]
      setcolorder(pop, c("number", "region", "age", "population"))
    } else {
      pop <- pop[number == country]
    }
    
    if (total){
      pop <- pop[age == "TOTAL"]
    } else {
      pop <- pop[age != "TOTAL"]
      pop <- pop[, age := as.integer(age)]
    }    
    return(pop)
  } else {
    ## TAG COUNTRY
    ## Import the population data from a country-specific URL (containing the number 
    ## of inhabitants per age group / NUTS3 region). The else{} option contain code
    ## if google or eurostat data is usable.
    ## The output of this if section depends on the value of total.
    ## WHEN TOTAL == TRUE: Return a data frame named "pop", containing the following variables
    # number: Region id in a format that can be matched with location_key in the case data table, via the get_reg_nb() function in function_utils.R
    # population: number of inhabitants per region
    ## WHEN TOTAL == FALSE: Return a data frame named "pop_long" containing the following variables:
    # number: Region id in a format that can be matched with location_key in the case data table, via the get_reg_nb() function in function_utils.R
    # region: region name
    # age: Age group
    # population: Number of inhabitants per region / age group
    if (country == "FR"){
      ## Downloaded from: "https://www.insee.fr/fr/statistiques/fichier/1893198/estim-pop-dep-sexe-aq-1975-2022.xlsx"
      ## Converted to a csv file with data from 2020
      pop <- as.data.table(read.csv2("Data/pop_struct.csv", sep = ","))
      colnames(pop) <- c("Number", "region", "0-4", "5-9", "10-14", "15-19", 
                         "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", 
                         "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", 
                         "80-84", "85-89", "90-94", "95+", "Total")
      if (total){
        pop <- pop[,.(number = Number, region = region, population = Total)]
        
        # Remove France total and Corsica
        pop <- pop[!(region == "France" | number %in% c("2A", "2B"))]
        return(pop)
      } else {
        # Take number of cases as reference age group distribution and match vaccine and pop to it.
        pop_long <- 
          tidyr::pivot_longer(pop, c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
                                     "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", 
                                     "60-64", "65-69", "70-74", "75-79", "80-84", "85-89",
                                     "90-94", "95+", "Total"), values_to = "Population", 
                              names_to = "Age") %>% as.data.table
        names(pop_long) <- tolower(names(pop_long))

        # Drop total populations and national population
        pop_long <- pop_long[!(age == "Total" | region == "France")]
        # Remove Corsica
        pop_long <- pop_long[!is.element(number, c("2A", "2B"))]
        return(pop_long)    
      }
    } else {
      # use Google open data (N.B. These are "the most recently reported" populations (e.g. from a census), so not all 2020)
      if (google){ 
        temp <- tempfile()
        download.file("https://storage.googleapis.com/covid19-open-data/v3/demographics.csv", temp)
        pop <- as.data.table(read.csv2(temp, sep = ","))
        # Merge with Google index of place names
        pop <- merge(pop, index[, .(location_key, subregion1_name)], by = "location_key", all.x = T)
        
        # Extract data for desired country
        pop <- pop[substr(location_key,1,2) == country, ]
        
        ## TAG COUNTRY
        ## If using google data, and certain regions are missing, add the missing rows here
        # Subset to subnational data
        if (country == "IT"){
          # Add rows with region "IT_88_SU", "IT_52_MS", "IT_72_NA", which are missing in pop
          new_rows <- data.table("location_key" = c("IT_88_SU", "IT_52_MS", "IT_72_NA"))
          pop <- rbindlist(list(pop, new_rows), fill = T) # Row-Wisely combine pop and new_rows
          # Compute the number of inhabitants in IT_88_SU, which is equal to the population in
          # IT_88 minus the sum of the number of inhabitants in all NUTS-3 regions in IT_88
          pop[location_key == "IT_88_SU", 
              population := pop[location_key == "IT_88", population] - 
                sum(pop[substr(location_key, 1, 5) == "IT_88" & location_key != "IT_88", population], 
                    na.rm = T)]
          pop[location_key == "IT_88_SU", subregion1_name := pop[location_key == "IT_88", subregion1_name]] 
          # Compute the number of inhabitants in IT_52_MS, which is equal to the population in
          # IT_52 minus the sum of the number of inhabitants in all NUTS-3 regions in IT_52
          pop[location_key == "IT_52_MS", 
              population := pop[location_key == "IT_52", population] - 
                sum(pop[substr(location_key, 1, 5) == "IT_52" & location_key != "IT_52", population], 
                    na.rm = T)]
          pop[location_key == "IT_52_MS", subregion1_name := pop[location_key == "IT_52", subregion1_name]] 
          # Compute the number of inhabitants in IT_72_NA, which is equal to the population in
          # IT_72 minus the sum of the number of inhabitants in all NUTS-3 regions in IT_72
          pop[location_key == "IT_72_NA", 
              population := pop[location_key == "IT_72", population] - 
                sum(pop[substr(location_key, 1, 5) == "IT_72" & location_key != "IT_72", population], 
                    na.rm = T)]
          pop[location_key == "IT_72_NA", subregion1_name := pop[location_key == "IT_72", subregion1_name]] 
          # Select all NUTS3 regions
          pop <- pop[nchar(location_key) > 6 | is.element(location_key, c("IT_23", "IT_BZ", "IT_TN"))]
          
          # Write correspondence between NUTS2 and NUTS3 regions
          convert <- c(
            "ITF35", "ITF43", "ITF44", "ITF45", "ITF46", "ITF47", "ITF48", "ITF51", "ITF52",
            "ITF61", "ITC4C", "ITC4D", "ITF11", "ITF12", "ITF13", "ITF14", "ITF21", "ITF22",
            "ITF31", "ITF32", "ITF33", "ITF34", "ITF62", "ITF63", "ITF64", "ITF65", "ITG11",
            "ITG12", "ITG13", "ITG14", "ITG15", "ITG16", "ITG17", "ITG18", "ITG19", "ITG2D",
            "ITG2E", "ITG2H", "ITG2G", "ITH10", "ITH20", "ITH31", "ITH32", "ITH33", "ITH34",
            "ITH35", "ITH36", "ITH37", "ITH41", "ITH42", "ITH43", "ITH44", "ITH51", "ITH52",
            "ITH53", "ITH54", "ITH55", "ITH56", "ITH57", "ITH58", "ITH59", "ITI11", "ITI12",
            "ITI13", "ITI14", "ITI15", "ITI16", "ITI17", "ITI18", "ITI19", "ITI1A", "ITI21",
            "ITI22", "ITI31", "ITI32", "ITI33", "ITI34", "ITI35", "ITI41", "ITI42", "ITI43",
            "ITI44", "ITI45", "ITC11", "ITC12", "ITC13", "ITC14", "ITC15", "ITC16", "ITC17",
            "ITC18", "ITC20", "ITC31", "ITC32", "ITC33", "ITC34", "ITC41", "ITC42", "ITC43",
            "ITC44", "ITC46", "ITC47", "ITC48", "ITC49", "ITC4A", "ITC4B", "ITG2F")
          names(convert) <- c(
            "IT_72_SA", "IT_75_TA", "IT_75_BR", "IT_75_LE", "IT_75_FG", "IT_75_BA",
            "IT_75_BT", "IT_77_PZ", "IT_77_MT", "IT_78_CS", "IT_25_MI", "IT_25_MB",
            "IT_65_AQ", "IT_65_TE", "IT_65_PE", "IT_65_CH", "IT_67_IS", "IT_67_CB",
            "IT_72_CE", "IT_72_BN", "IT_72_NA", "IT_72_AV", "IT_78_KR", "IT_78_CZ",   
            "IT_78_VV", "IT_78_RC", "IT_82_TP", "IT_82_PA", "IT_82_ME", "IT_82_AG",
            "IT_82_CL", "IT_82_EN", "IT_82_CT", "IT_82_RG", "IT_82_SR", "IT_88_SS",
            "IT_88_NU", "IT_88_SU", "IT_88_OR", "IT_BZ", "IT_TN", "IT_34_VR",   
            "IT_34_VI", "IT_34_BL", "IT_34_TV", "IT_34_VE", "IT_34_PD", "IT_34_RO",
            "IT_36_PN", "IT_36_UD", "IT_36_GO", "IT_36_TS", "IT_45_PC", "IT_45_PR",
            "IT_45_RE", "IT_45_MO", "IT_45_BO", "IT_45_FE", "IT_45_RA", "IT_45_FC",
            "IT_45_RN", "IT_52_MS", "IT_52_LU", "IT_52_PT", "IT_52_FI", "IT_52_PO",
            "IT_52_LI", "IT_52_PI", "IT_52_AR", "IT_52_SI", "IT_52_GR", "IT_55_PG",
            "IT_55_TR", "IT_57_PU", "IT_57_AN", "IT_57_MC", "IT_57_AP", "IT_57_FM",
            "IT_62_VT", "IT_62_RI", "IT_62_RM", "IT_62_LT", "IT_62_FR", "IT_21_TO",
            "IT_21_VC", "IT_21_BI", "IT_21_VB", "IT_21_NO", "IT_21_CN", "IT_21_AT",
            "IT_21_AL", "IT_23", "IT_42_IM", "IT_42_SV", "IT_42_GE", "IT_42_SP",   
            "IT_25_VA", "IT_25_CO", "IT_25_LC", "IT_25_SO", "IT_25_BG",   "IT_25_BS",
            "IT_25_PV", "IT_25_LO", "IT_25_CR", "IT_25_MN", "IT_88_CA")
          pop[, location_key := convert[location_key]]
          
        } else{
          pop <- pop[between(nchar(location_key), 4, 6)]
        }
        
        if (total){
          # Select column containing location key and number of inhabitants
          pop <- pop[, .(location_key, population)]
          colnames(pop) <- c("number", "population")
          return(pop)
        } else {
          pop_long <- melt(pop, id.vars = c("location_key","subregion1_name"),
                           measure.vars = patterns("population_age_"),
                           variable.name = "age",
                           value.name = "population")
          pop_long[, age := sub("_","-",sub("population_age_","",sub("_and_older","+",age)))]
          pop_long[age == "00-09", age := "0-9"]
          return(pop_long)
        }
        unlink(temp)
      } else { # use Eurostat
        # Read in population data at NUTS3 level
        pop <- as.data.table(read.csv2("Data/demo_r_pjangrp3_1_Data.csv", sep = ",", encoding = "latin1"))
        
        names(pop) <- tolower(names(pop))
        
        # Drop unneeded columns
        pop <- pop[, .(number = geo, region = geo_label, age, population = value)]
        
        pop[, population := as.numeric(gsub(",","",population))]
        # Drop population with unknown age as it's very small
        pop <- pop[age != "UNK"]
        pop[, age := sub("_","",sub("Y","",age))]
        pop <- pop[age != "GE85"]
        pop[age == "LT5", age := "0-4"]
        pop[age == "GE90", age := "90+"]
        
        # Select data for desired country
        pop <- pop[substr(number,1,2) == country & nchar(number) > 4]
        
        # Get region number
        pop[, number := get_reg_nb(number, country)]
        if (total){
          pop <- pop[age == "TOTAL"]
        } else {
          pop <- pop[age != "TOTAL"]
        }
        
        return(pop)
      }
    }
  }
}

#### Import map
import_map <- function(country, url = F, index, dt_incidence){
  dir <- tempdir()
  if(url == FALSE){
    temp <- "Data/NUTS_RG_20M_2021_3035.shp.zip"
    path <- paste0(dir, "/",sub(".*/","",sub(".zip","",temp)))
  } else{
    temp <- tempfile()
    download.file(url, temp)
    path <- paste0(dir, "/",sub(".*/","",sub(".zip","",url)))
  }
  
  # Extract the shapefile into a temporary folder
  unzip(temp, exdir = dir)
  # Read the shapefile
  map <- sf::st_read(path,quiet = T)
  # Transform map to WGS84 coordinate reference system
  map <- sf::st_transform(map,"WGS84")
  # Subset map to country of interest
  map1 <- map[map$CNTR_CODE == country,]
  ## Add the region number to map1
  ## TAG COUNTRY
  ## Create column reg_nb in map1, which corresponds to the reference of the region
  ## in the "index" database. To do so, first remove overseas territories (if any),
  ## Select the entries of interest in index, set a key in index, and add subregion code 
  ## to map1 (in a column named "reg_nb")
  if (country == "FR"){
    # Exclude French overseas territories
    map1 <- map1[-(grep("FRY", map1$NUTS_ID)),]
    # Remove Corsica
    map1 <- map1[-(grep("FRM", map1$NUTS_ID)),]
    
    # Select entries from subregion2 level in Google data
    index <- index[country_code == country & nchar(location_key) > 6, 
                   .(subregion2_name, subregion2_code)] 
    index$subregion2_name <- stringi::stri_trans_general(
      index$subregion2_name, id = "Latin-ASCII")
    map1$NUTS_NAME <- stringi::stri_trans_general(map1$NUTS_NAME, id = "Latin-ASCII")
    setkey(index, subregion2_name)
    # Add department number to map1
    map1$reg_nb <- index[map1$NUTS_NAME, subregion2_code]
  } else if (country == "CZ"){
    # Select entries from subregion1 level in Google data
    index <- index[country_code == country & nchar(location_key) == 5,
                   .(NUTS_ID = sub("nuts/","",datacommons_id), subregion1_code)]
    setkey(index, NUTS_ID)
    # Add region number to map1
    map1$reg_nb <- index[map1$NUTS_ID, subregion1_code]
  } else if (country == "IT"){
    map1$reg_nb <- map1$FID
  }
  
  # Delete the temporary file
  if(url != FALSE) unlink(temp)
  unlink(dir)
  key_incidence <- unique(dt_incidence$location_key)
  names(key_incidence) <- get_reg_nb(key_incidence, country)
  map1$key <- key_incidence[map1$reg_nb]
  return(map1)
}

#### Contact data, from socialmixr or Prem et al (2021) 
import_contact <- function(country, index, age_groups, pop_age, pop){
  pop_by_age <- colSums(pop_age)
  min_ages <- get_min_age(age_groups)
  nb_age <- length(age_groups)
  
  ## TAG COUNTRY:
  ## AGE-STRATIFIED MODEL:
  ## Import a contact matrix from a country-specific survey, or multi-country analysis.
  ## Must return a matrix named C (dimensions 9 * 9), containing the number of contacts
  ## between the different age groups (0-10 years old; 10-20;...; 80+).
  if (country == "FR"){
    # Import survey
    survey <- get_survey("https://doi.org/10.5281/zenodo.1157918")
    
    if (nb_age > 9) max_age <- 9 else max_age <- nb_age
    # Calculate contact matrix
    m <- contact_matrix(survey,
                        countries = "France",
                        age.limits = min_ages[seq_len(max_age)],
                        symmetric = T,
                        missing.participant.age = "remove",
                        missing.contact.age = "remove")
    C <- m$matrix
  } else if (country == "CZ"){
    # Read in synthetic contact matrices (Prem et al, 2021)
    contact_matrices <- fread("Data/synthetic_contacts_2020.csv")
    # Get contact matrix as data table for required country
    contact <- contact_matrices[iso3c == index[match(country,index[,country_code]),iso_3166_1_alpha_3] &
                                  setting == "overall" & location_contact == "all"]
    # Correct variable name
    setnames(contact,"age_cotactee","age_contactee")
    
    # Merge with population
    age_cols <- names(contact)[grep("age", names(contact))]
    
    contact[, (paste0("min_",age_cols)):=lapply(.SD,function(x) as.numeric(sub("\\+","",sub(" to.*","",x)))),.SDcols = age_cols]
    
    # Sum mean numbers of contacts over contact age groups being aggregated
    contact[, age_contactee := cut(min_age_contactee, c(min_ages, Inf), labels = age_groups, right = F)]
    contact1 <- contact[,.(contacts = sum(mean_number_of_contacts)),by = .(age_contactor,age_contactee,min_age_contactor)]
    
    min_ages_contact <- contact[,unique(min_age_contactor)]
    age_groups_contact <- contact[,unique(age_contactor)]
    
    pop_contact <- copy(pop)
    pop_contact[,age_group_contact := cut(age,c(min_ages_contact,Inf),labels = age_groups_contact,right = F)]
    pop_contact <- pop_contact[,.(population = sum(population)),by = .(age_group_contact)]
    
    contact1 <- merge(contact1,pop_contact,by.x = "age_contactor",by.y = "age_group_contact")
    contact1[,age_contactor := cut(min_age_contactor,c(min_ages,Inf),labels = age_groups,right = F)]
    
    # Take population-weighted average of mean number of contacts over contactor age groups being aggregated
    contact2 <- contact1[,.(contacts = sum(contacts * population)/sum(population)), by = .(age_contactor,age_contactee)]
    
    # Convert data table to matrix
    nb_age1 <- contact2[,length(unique(age_contactor))]
    C <- matrix(contact2[,contacts], nrow = nb_age1, ncol = nb_age1, byrow = T) # N.B. contactors in rows, contactees in columns
  }
  
  missing_age <- nb_age - nrow(C)
  sum_pop_missing_age <- sum(pop_by_age[seq(nb_age - missing_age, nb_age)])
  # Split contacts in final age group across final and missing age groups
  if (missing_age != 0){
    C <- C[c(1:(nb_age-1),rep(nb_age-1,missing_age)),]
    C <- C[,c(1:(nb_age-1),rep(nb_age-1,missing_age))]
    for(i in seq(nb_age - missing_age, nb_age)){
      C[, i] <- C[, i] * (pop_by_age[i] / sum_pop_missing_age)
    }
  }
  
  # Rename rows and columns
  colnames(C) <- age_groups
  rownames(C) <- age_groups
  
  return(C)
}

#### Import urban/rural status of each NUTS 3 region from Eurostat data 
import_urban_rural <- function(country){
  urban <- as.data.table(read_xlsx("Data/NUTS2021.xlsx", sheet = "Urban-rural remoteness"))
  
  urban <- urban[, .(`NUTS ID`,`NUTS LABEL`,`CATEGORY CODE`,`CATEGORY LABEL`)]
  names(urban) <- c("NUTS","name","code","meaning")
  
  # Select country of interest
  urban <- urban[substr(NUTS,1,2) == country]
  
  return(urban)
}

#### Import test data where it's not available in the Google data
import_test <- function(cntry, pop, age_groups = NULL){
  ## TAG COUNTRY
  ## !! If importing a country-specific test database, not needed otherwise
  ## download the file and process
  ## the data to create a data frame named "test", containing the following variables:
  # date: date in format YYYY-MM-DD
  # age: age group (10 year age band)
  # population: Number of inhabitants
  # nb_tests: Number of tests on that date / age group
  if (cntry == "IT"){
    temp <- tempfile()
    download.file("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv", temp)
    test <- as.data.table(read.csv2(temp, sep = ","))
    test <- test[, .(date = as.Date(data), age = "Total", nb_tests = c(0, diff(tamponi)), 
                     population = pop$population)]
    
    # Fill in numbers of tests in any missing weeks with average of preceding and following weeks (if available)
    mssng_dates <- test[is.na(nb_tests), date]
    for (i in seq_along(mssng_dates)){
      mssng_date <- mssng_dates[i]
      if ((mssng_date - 7) %in% test[,date] & (mssng_date + 7) %in% test[,date]){
        test[date == mssng_date, nb_tests := (test[date == mssng_date - 7, nb_tests] + test[date == mssng_date + 7, nb_tests])/2]
      } else if ((mssng_date - 7) %in% test[,date] & !((mssng_date + 7) %in% test[,date])){
        test[date == mssng_date, nb_tests := test[date == mssng_date - 7, nb_tests]]
      } else {
        test[date == mssng_date, nb_tests := test[date == mssng_date + 7, nb_tests]]
      }
    }
    
    # Divide tests by number of days
    test[, nb_tests := nb_tests/.N, by = .(date)]
    
    unlink(temp)
  } else if (cntry == "CZ"){
    temp <- tempfile()
    download.file("https://onemocneni-aktualne.mzcr.cz/api/v2/covid-19/nakazeni-hospitalizace-testy.csv", temp)
    test <- as.data.table(read.csv2(temp, sep = ","))
    
    # Drop unneeded columns
    test[, c("id", "nove_hospitalizace", "nove_jip") := NULL]
    
    # Rename columns
    setnames(test, 
             c("datum", "vekova_kategorie", "indikace_testu", "dokoncene_ockovani",
               "reinfekce", "provedene_testy", "potvrzene_pripady"), 
             c("date", "age_group", "test_indication", "full_vaccination",
               "reinfection", "tests", "confirmed_cases"))
    
    # Sum tests by age group
    test <- test[, .(tests = sum(tests, na.rm = T)), by = .(date, age_group)]
    
    # Remove tests with unknown age as they are a very small proportion (0.03%)
    test[age_group == "-", sum(tests)]/test[, sum(tests)] #0.0003023779
    test <- test[age_group != "-"]
    if(all(age_groups == "Total")){
      base <- CJ(date = test[,unique(date)], age = "TOTAL")
      base <- merge(base, pop[, .(age, population)], by = "age")
      base[, age_group := "TOTAL"]
      test <- test[, lapply(.SD, sum), by = date, .SDcols = "tests"]
      test[, age_group := "TOTAL"]
      # N.B. This duplicates tests in the same age group on each date
      test <- merge(base, test, by = c("date", "age_group"), all.x = T)
      test[, age_group := NULL]
      setnames(test, c("tests"), c("nb_tests"))
    } else{
      ## Redistribute tests by model age groups
      # Make base data table with all ages and dates for merging test data into
      base <- CJ(date = test[,unique(date)], age = 0:100)
      # Merge population data
      base <- merge(base, pop[, .(age, population)], by = "age")
      
      # Merge test data
      age_groups_test <- test[, sort(unique(age_group))]
      min_ages_test <- get_min_age(age_groups_test)
      base[, age_group := cut(age, c(min_ages_test, Inf), label = age_groups_test, right = F)]
      # N.B. This duplicates tests in the same age group on each date
      test <- merge(base, test, by = c("date", "age_group"), all.x = T)
      # Change 'tests' to numeric type
      test[, tests := as.numeric(tests)]
      # Divide tests by population fraction
      test[, tests := tests * population/sum(population), by = .(date, age_group)]
      
      # Change age groups
      min_ages <- get_min_age(age_groups)
      test[, age_group := cut(age, c(min_ages, Inf), labels = age_groups, right = F)]
      # Sum tests by model age groups
      cols <- c("population", "tests")
      test <- test[, lapply(.SD, sum), .SDcols = cols, by = .(date, age_group)]
      # Change column names
      setnames(test, c("age_group","tests"), c("age","nb_tests"))
    }
    unlink(temp)
  } else{
    temp <- tempfile()
    download.file("https://opendata.ecdc.europa.eu/covid19/testing/csv/data.csv", temp)
    test <- as.data.table(read.csv2(temp, sep = ","))
    
    # Select data for country of interest
    test <- test[country_code == cntry & level == "national", 
                 .(country = country_code, date_week = year_week, nb_tests = tests_done)]
    # Create date variable
    test[, date := ISOweek2date(paste0(date_week,"-1"))]
    
    # Create data table with all dates
    dates <- seq.Date(test[, min(date)], test[, max(date) + 6], by = 1)
    base <- data.table(date = dates, date_week = ISOweek(dates))
    # Merge test data
    # N.B. This duplicates doses for dates within the same ISO week
    test <- merge(base, test[,!"date"], by = "date_week", all.x = T)
    
    # Fill in numbers of tests in any missing weeks with average of preceding and following weeks (if available)
    test[, nb_tests := as.numeric(nb_tests)]
    mssng_dates <- test[is.na(nb_tests), date]
    for (i in seq_along(mssng_dates)){
      mssng_date <- mssng_dates[i]
      if ((mssng_date - 7) %in% test[,date] & (mssng_date + 7) %in% test[,date]){
        test[date == mssng_date, nb_tests := (test[date == mssng_date - 7, nb_tests] + test[date == mssng_date + 7, nb_tests])/2]
      } else if ((mssng_date - 7) %in% test[,date] & !((mssng_date + 7) %in% test[,date])){
        test[date == mssng_date, nb_tests := test[date == mssng_date - 7, nb_tests]]
      } else {
        test[date == mssng_date, nb_tests := test[date == mssng_date + 7, nb_tests]]
      }
    }
    
    # Divide tests by number of days
    test[, nb_tests := nb_tests/.N, by = .(date_week)]
    
    # Add dummy age group column
    test[, age := "Total"]
    
    unlink(temp)
  }
  return(test)
}
