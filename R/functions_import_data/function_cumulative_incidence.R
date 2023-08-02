#### Cumulated incidence

### Compute the cumulative incidence per age group, per region, at each wave
### We considered three different waves: Wild Type + alpha ; Delta ; Omicron
cumu_incidence <- function(country, dt_variant, dt_case, dt_pop, prop_delta = 0.3, 
                           prop_omicron = 0.3, total = F){
  ## Define the time window of each wave
  start_delta <- dt_variant[variant == "Delta" & percent_variant >= prop_delta, 
                            year_week][1]
  start_omicron <- dt_variant[variant == "Omicron" & percent_variant >= prop_omicron, 
                              year_week][1]
  
  if (total){
    ## Select columns date, location, week of onset, Age, and number of cases
    case_nb <- dt_case[, .(date, location_key, week_day, nb_cases = new_confirmed, nb_tests = new_tested)]
    
    ## Use the same location format as the population and vaccine datasets 
    case_nb[, region_nb := get_reg_nb(location_key, country)]
    
    ## Add the column "variant", which describes what wave each date belongs to.
    case_nb[week_day < start_delta, variant := "WT+alpha"]
    case_nb[week_day >= start_delta & week_day < start_omicron, variant := "delta"]
    case_nb[week_day >= start_omicron, variant := "omicron"]
    
    ## Add the column "cumu_cases": the cumulative number of cases by age, location, and wave
    case_nb[, cumu_cases := cumsum(nb_cases), by = .(location_key, variant)]
    
    # Compute cumulative incidence per region, age group, and variant
    case_nb <- merge(case_nb, dt_pop, by.x = "region_nb", by.y = "number", all.x=T)
    case_nb[, cumu_incidence := 1e5 * cumu_cases / population]
    setnames(case_nb,"region_nb","region")
    
    ## Add dummy age group for reshaping into case array
    case_nb[,age := "Total"]
    
    return(case_nb)
  } else {
    ## Convert dt_case into long format:
    long_case_nb <- 
      tidyr::pivot_longer(dt_case, grep("new_confirmed", colnames(dt_case)), 
                          names_to  = "age_group", values_to = "nb_cases") %>% as.data.table
    long_test_nb <- 
      tidyr::pivot_longer(dt_case, grep("new_tested", colnames(dt_case)), 
                          names_to  = "age_group", values_to = "nb_tests") %>% as.data.table
    long_case_nb[, nb_tests := long_test_nb$nb_tests]
    
    ## Select columns date, location, week of onset, Age, and number of cases
    long_case_nb <- long_case_nb[, .(date, location_key, week_day, age_group, nb_cases, nb_tests)]
    
    ## Use the same location format as the population and vaccine datasets 
    long_case_nb[, region_nb := get_reg_nb(location_key, country)]
    
    ## Currently, the "age" variable contains the age bin reference from dt_case, 
    ## we convert this column into the actual age group using the "age_bin" columns in dt_case
    age_convert <- as.character(dt_case[1, grep("age_bin", colnames(dt_case)), with = F])
    names(age_convert) <- unique(long_case_nb$age_group)
    long_case_nb[, age_group := age_convert[age_group]]
    
    ## Add the column "variant", which describes what wave each date belongs to.
    long_case_nb[week_day < start_delta, variant := "WT+alpha"]
    long_case_nb[week_day >= start_delta & week_day < start_omicron, variant := "delta"]
    long_case_nb[week_day >= start_omicron, variant := "omicron"]
    
    ## Add the column "cumu_cases": the cumulative number of cases by age, location, and wave
    long_case_nb[, cumu_cases := cumsum(nb_cases), by = .(age_group, location_key, variant)]
    
    ## Aggregate population in dt_pop into 10-year age groups in dt_case
    min_ages <- get_min_age(age_convert)
    dt_pop1 <- copy(dt_pop)
    dt_pop1[,age_group:=cut(get_min_age(age),c(min_ages,Inf),labels=age_convert,right=F)]
    dt_pop_agg <- dt_pop1[,.(population=sum(population)),by=.(region_nb=number,age_group)]
    
    # Compute cumulative incidence per region, age group, and variant
    long_case_nb <- merge(long_case_nb,dt_pop_agg,by=c("region_nb","age_group"),all.x=T)
    long_case_nb[, cumu_incidence := 1e5 * cumu_cases / population]
    long_case_nb[age_group == "90-", age_group := "90+"]
    setnames(long_case_nb,c("region_nb","age_group"),c("region","age"))
    
    return(long_case_nb)
  }
}
