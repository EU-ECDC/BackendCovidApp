#### Match age structure ####

### Summary:
## Match the age structure of population and vaccine coverage to the incidence data
## Example with French data:
## Local population in vacc_FR does not exactly match pop (<2% difference when comparable), 
## we use pop to match the age group in vacc_FR to case_by_age_FR
## Correspondence between new (left) and old (right) age groups:
## Main assumption: Homogeneous coverage with age group (e.g vacc_FR0 = vacc_FR1 = .. = vacc_FR0-11) 
# vacc_FR0-10 = vacc_FR0-11
# vacc_FR10-20 = vacc_FR12-17
# vacc_FR20-30 = (vacc_FR18-24 * pop20-24 + vacc_FR25-29 * pop25-29) / (pop20-24 + pop25-29)
# vacc_FR30-39 = vacc_FR25-39
# vacc_FR40-49 = vacc_FR40-54
# vacc_FR50-59 = (vacc_FR40-54 * pop50-55 + vacc_FR55-64 * pop55-59) / (pop50-54 + pop55-59)
# vacc_FR60-69 = (vacc_FR55-64 * pop60-64 + vacc_FR65-74 * pop65-69) / (pop60-64 + pop65-59)
# vacc_FR70-79 = (vacc_FR65-74 * pop70-74 + vacc_FR75+ * pop75-79) / (pop70-74 + pop75-79)
# vacc_FR80-79 = vacc_FR75+
# vacc_FR90+ = vacc_FR75+

adjust_age_group <- function(vacc_reg, pop_reg, age_groups = NULL, vacc_ecdc = F, 
                             total = F, national = F, country = NULL){
  if (is.null(vacc_reg)){
    # If input is NULL, set output to NULL
    vacc_long <- NULL
  }
  ## Set ID in pop_reg
  pop_df <- copy(pop_reg)
  if(any(colnames(pop_df) == "age")){
    pop_df[, id := paste0(number, "-", age)]
    setkey(pop_df, id)
  }
  
  ## If national level, divide the nationwide number of doses by the number of inhabitants
  ## If total, compute across all age groups, otherwise stratify by age groups
  ## if vacc_ecdc, the vaccine data was taken from the ECDC vaccine tracker. 
  
  ## All national data were taken from ECDC.
  ## No local age stratified data in ECDC vaccine database (equivalent to vacc_ecdc == T, 
  ## total == F and national == F)
  if (national){
    if (total){
      vacc_output <- copy(vacc_reg)
      
      # Calculate proportion of population vaccinated for each dose
      cols <- c("dose1", "dose2", "dose3")
      vacc_output[, (cols) := lapply(.SD, function(x) x/pop_df[,population]), .SDcols = cols]
      
      # Melt to long format
      vacc_long <- melt(vacc_output, measure.vars = cols, variable.name = "dose", value.name = "prop_dose")
      
      # Rename column
      setnames(vacc_long, "Age", "age")
      
      # Calculate coverage
      vacc_long[, cumu_dose := cumsum(prop_dose), by = .(age, dose)]
    } else {
      min_ages <- get_min_age(age_groups)
      
      # Select columns of interest in population data table
      pop_df <- pop_df[,.(age,population)]

      ## data frame containing all combination of date_week and age
      base <- CJ(date_week = vacc_reg[,unique(date_week)], age = pop_df[,age])
      base <- merge(base, pop_df, by = "age")
      
      ## Remove entries of vacc_reg with unknown age
      vacc_reg <- vacc_reg[Age != "ALL"]
      vacc_reg <- vacc_reg[Age != "UNK"]
      age_groups_vax <- vacc_reg[, sort(unique(Age))]
      min_ages_vax <- get_min_age(age_groups_vax)
      
      base[,age_group := cut(age, c(min_ages_vax,Inf), labels = age_groups_vax, right = F)]
      
      # Merge with base data table
      vacc_output <- merge(base, vacc_reg[,.(date_week,Age,dose1,dose2,dose3)], by.x = c("date_week", "age_group"), by.y = c("date_week", "Age"), all.x = T)
      # Fill missing values with 0s
      cols <- c("dose1","dose2","dose3")
      setnafill(vacc_output, fill = 0, cols = cols)
      # Convert dose columns to numeric
      vacc_output[,(cols) := lapply(.SD, as.numeric), .SDcols = cols]
      # Split doses by population fraction
      vacc_output[,(cols) := lapply(.SD, function(x) x*population/sum(population)), .SDcols = cols, by = .(date_week,age_group)]
      
      # Change age groups
      vacc_output[, age_group := cut(age, c(min_ages,Inf), labels = age_groups, right = F)]
      
      # Sum doses over new age groups
      vacc_output <- vacc_output[, lapply(.SD, sum), .SDcols = c("population",cols), by = .(date_week, age_group)]
      vacc_output[, (cols) := lapply(.SD, function(x) x/population), .SDcols = cols]
      
      vacc_long <- melt(vacc_output, measure.vars = cols, variable.name = "dose", value.name = "prop_dose")
      
      # Rename column
      setnames(vacc_long, "age_group", "age")
      
      # Calculate coverage
      vacc_long[, cumu_dose := cumsum(prop_dose), by = .(age, dose)]
    }
  } else{
    if (vacc_ecdc){
      if (total){
        # Group population by NUTS2:
        pop_nuts2 <- aggregate(pop_df$population, list(substr(pop_df$number, 1, 4)), sum) %>% 
          as.data.table
        colnames(pop_nuts2) <- colnames(pop_df)
        setkey(pop_nuts2, "number")
        
        vacc_output <- copy(vacc_reg)
        
        # Calculate proportion of population vaccinated for each dose
        cols <- c("dose1", "dose2", "dose3")
        vacc_output[, (cols) := lapply(.SD, function(x) x/pop_nuts2[number,population]), 
                    .SDcols = cols]
        # Convert to NUTS3
        vacc_output_nuts3 <- data.table(
          date_week = rep(unique(vacc_output$date_week), each = length(pop_df$number)),
          region_nb = rep(pop_df$number, length(unique(vacc_output$date_week))),
          Age = "Total", dose1 = 0, dose2 = 0, dose3 = 0
        )
        
        vacc_output[, ID := paste0(date_week, "_", number)]
        vacc_output_nuts3[, ID := paste0(date_week, "_", substr(region_nb, 1, 4))]
        setkey(vacc_output, "ID")
        vacc_output_nuts3[, dose1 := vacc_output[ID, dose1]]
        vacc_output_nuts3[, dose2 := vacc_output[ID, dose2]]
        vacc_output_nuts3[, dose3 := vacc_output[ID, dose3]]
        vacc_output_nuts3[is.na(vacc_output_nuts3)] <- 0
        vacc_output_nuts3[, ID := NULL]
        
        # Melt to long format
        vacc_long <- melt(vacc_output_nuts3, measure.vars = cols, variable.name = "dose", value.name = "prop_dose")
        
        # Rename column
        setnames(vacc_long, "Age", "age")
        # Calculate coverage
        vacc_long[, cumu_dose := cumsum(prop_dose), by = .(region_nb, age, dose)]
        colnames(vacc_long) <- c("date_week", "region_nb", "age", "dose", 
                                 "prop_dose", "cumu_dose")
      } else{
        stop("Local age-stratified ECDC vaccine data is not implemented")
      }
    } else{
      ## TAG COUNTRY
      ## AGE STRATIFIED MODELS
      ## Create a data table named vacc_long, which contains the proportion of the population 
      ## per region / age group who received each dose of vaccine (on that day and cumulative).
      ## vacc_long contains the following variables:
      # date_week: date in format YYYY-MM-DD
      # region_nb: Region code
      # age: Age group described in this row
      # dose: Dose number
      # prop_dose: Proportion of the population vaccinated that day
      # cumu_dose: Cumulative proportion of the population vaccinated
      if (country == "FR"){
        vacc_reg_prop <- copy(vacc_reg)
        # Calculate proportions vaccinated with different doses
        vacc_reg_prop[, dose1_prop := dose1 / Population]
        vacc_reg_prop[, dose2_prop := dose2 / Population]
        vacc_reg_prop[, dose3_prop := dose3 / Population]
        
        if (total){
          vacc_output <- vacc_reg_prop[, .(date_week, number, dose1_prop, dose2_prop, 
                                           dose3_prop)] 
          
          # Rename columns
          colnames(vacc_output) <- c("date_week", "region_nb", "dose1", "dose2", 
                                     "dose3")
          # Convert to long format
          vacc_long <- tidyr::pivot_longer(vacc_output, c("dose1", "dose2", "dose3"), 
                                           names_to = "dose", values_to = "prop_dose") %>% 
            as.data.table
          # Calculate vaccine coverage
          vacc_long[, cumu_dose := cumsum(prop_dose), by = .(region_nb, dose)]        
        } else{
          ### Data initialisation
          ## We create a new vaccine data table which will contain the coverage in new age groups
          vacc_output <- data.table()
          
          ## One row per new age group, region, and date
          nb_age <- length(age_groups)
          nb_dep <- length(unique(vacc_reg_prop$number))
          nb_dates <- length(unique(vacc_reg_prop$date_week))
          
          ## Initialise the values of dates, region and age
          vacc_output[, date_week := rep(unique(vacc_reg_prop$date_week), each = nb_dep * nb_age)]
          vacc_output[, dep := rep(unique(vacc_reg_prop$number), nb_age * nb_dates)]
          vacc_output[, age_group := rep(rep(age_groups, each = nb_dep), nb_dates)]
          ## Create a unique ID for each row, and use it to create a key
          vacc_output[, ID := paste0(date_week, "-", dep, "-", age_group)]
          setkey(vacc_output, ID)
          ## age group from old vaccine uptake data and population used to compute the
          ## new age groups. If only one age group is needed (e.g. 0-11 becomes 0-5), 
          ## vacc_age2, pop_age1, and pop_age2 are set to NA, since only the old vaccine uptake
          ## will be needed
          vacc_age1 <- c("00-11", "12-17", "18-24", "25-39", "40-54", "40-54", "55-64",
                         "65-74", "75+")#, "75+")
          vacc_age2 <- c(NA, NA, "25-39", NA, NA, "55-64", "65-74", "75+", NA)#, NA)
          pop_age1 <- c(NA, NA, "20-24", NA, NA, "50-54", "60-64", "70-74", NA)#, NA)
          pop_age2 <- c(NA, NA, "25-29", NA, NA, "55-59", "65-69", "75-79", NA)#, NA)
          # Link vacc_age1, vacc_age2, pop_age1, and pop_age2 to the new age groups
          names(vacc_age1) <- names(vacc_age2) <- names(pop_age1) <- names(pop_age2) <- 
            age_groups
          ### Computing the new coverage:
          ## Add columns with the number of inhabitants corresponding to pop_age1 and pop_age2
          ## in each region, using the key in pop_df
          vacc_output[, pop1 := pop_df[paste(dep, pop_age1[age_group], sep = "-"), population]]
          vacc_output[, pop2 := pop_df[paste(dep,  pop_age2[age_group], sep = "-"), population]]
          
          ## Add columns with the vaccine coverage (for each dose), corresponding to vacc_age1 and 
          ## vacc_age2, in each region at each date, using the key in vacc_reg_prop
          vacc_output[, id1:=paste(date_week, dep, vacc_age1[age_group], sep = "-")]
          vacc_output[, c("vacc_old1d1", "vacc_old1d2", "vacc_old1d3") := 
                        vacc_reg_prop[id1,.(dose1_prop, dose2_prop, dose3_prop)]]
          vacc_output[, id2:=paste(date_week, dep, vacc_age2[age_group], sep = "-")]
          vacc_output[, c("vacc_old2d1", "vacc_old2d2", "vacc_old2d3") := 
                        vacc_reg_prop[id2,.(dose1_prop, dose2_prop, dose3_prop)]]

          ## Compute the new coverage in each region, using the equations from the summary
          vacc_output[, dose1 := (vacc_old1d1 * pop1 + vacc_old2d1 * pop2)/ (pop1 + pop2)]
          vacc_output[, dose2 := (vacc_old1d2 * pop1 + vacc_old2d2 * pop2)/ (pop1 + pop2)]
          vacc_output[, dose3 := (vacc_old1d3 * pop1 + vacc_old2d3 * pop2)/ (pop1 + pop2)]
          ## If only one age group was needed to compute dose1, recalculate the new coverage 
          vacc_output[is.na(vacc_age2[age_group]), dose1 := vacc_old1d1]
          vacc_output[is.na(vacc_age2[age_group]), dose2 := vacc_old1d2]
          vacc_output[is.na(vacc_age2[age_group]), dose3 := vacc_old1d3]

          
          ### Solve problem with the 10-19 age group
          ## In the current calculation, the coverage at 10-19 is equal to the coverage at
          ## 12-17, which is a poor estimation (10-11, 12-17 and 18-19 had access to vaccine)
          ## at different dates. We implement a different method to compute the coverage at 
          ## 10-19 year old:
          ## nbvacc1011 = (2/5 * pop1014 * nbvacc011) / (pop0004 + pop0509 + 2/5 * pop1014)
          ## nbvacc1819 = (2/5 * pop1519 * nbvacc1824) / (2/5 * pop1519 + pop2024)
          ## vacc1019 = (nbvacc1011 + nbvacc1217 + nbvacc1819) / (pop1014 + pop1519)
          ## Assumption: homogeneous number of inhabitants within the age groups
          
          ## First, we create a data table that will contain the number of inhabitants per age group
          pop1019 <- data.table()
          # pop0111 = (pop0004 + pop0509 + 2/5 * pop1014)
          pop1019[, pop0111 := (pop_df[number != 0 & age == "0-4", population] +
                                  pop_df[number != 0 & age == "5-9", population] +
                                  2/5 * pop_df[number != 0 & age == "10-14", population])]
          # pop1011 = 2/5 * pop1014
          pop1019[, pop1011 := 2/5 * pop_df[number != 0 & age == "10-14", population]]
          # prop1011 = (2/5 * pop1014) / (pop0004 + pop0509 + 2/5 * pop1014)
          pop1019[, prop1011 := pop1011 / pop0111]
          # pop1824 = (2/5 * pop1519 + pop2024)
          pop1019[, pop1824 := (2/5 * pop_df[number != 0 & age == "15-19", population] +
                                  pop_df[number != 0 & age == "20-24", population])]
          # pop1819 = 2/5 * pop1519
          pop1019[, pop1819 := 2/5 * pop_df[number != 0 & age == "15-19", population]]
          # prop1819 = (2/5 * pop1519) / (2/5 * pop1519 + pop2024)
          pop1019[, prop1819 := pop1819 / pop1824]
          # pop1019 = pop1014 + pop1519
          pop1019[, pop1019 := pop_df[number != 0 & age == "10-14", population] +
                    pop_df[number != 0 & age == "15-19", population]]
          ## Set an ID for the new dataset, using the region
          pop1019[, ID := pop_df[number != 0 & age == "15-19", number]]
          setkey(pop1019, ID)
          
          ## Compute the vacc_reg_prop ID required for each age group we need (0-11, 12-17, and 18-24)
          ## id_1019 is a matrix, containing the ID per age group
          id_1019 <- t(apply(vacc_output[age_group == "10-19", .(date_week, dep)], 1, 
                             function(X) paste0(paste(X, collapse = "-"), 
                                                c("-00-11", "-12-17", "-18-24"))))
          # Compute the coverage for each dose
          vacc_output[age_group == "10-19", c("dose1", "dose2", "dose3") := 
                        (vacc_reg_prop[id_1019[,1], .(dose1, dose2, dose3)] * pop1019[dep, prop1011] +
                           vacc_reg_prop[id_1019[,2], .(dose1, dose2, dose3)] +
                           vacc_reg_prop[id_1019[,3], .(dose1, dose2, dose3)] * pop1019[dep, prop1819])/
                        pop1019[dep, pop1019]]
          
          # Select columns of interest
          vacc_output <- vacc_output[, .(date_week, dep, age_group, dose1, dose2, dose3)] 
          # Rename columns
          colnames(vacc_output) <- c("date_week", "region_nb", "age", "dose1", "dose2", 
                                     "dose3")
          vacc_long <- tidyr::pivot_longer(vacc_output, c("dose1", "dose2", "dose3"), 
                                           names_to = "dose", values_to = "prop_dose") %>% 
            as.data.table
          vacc_long[, cumu_dose := cumsum(prop_dose), by = .(age, region_nb, dose)]        
        }
      } else if (country == "CZ"){
        if (total){
          # Select columns of interest
          vacc_output <- vacc_reg[, .(date_week, number, dose1, dose2, dose3)]
          # Merge with population data
          vacc_output <- merge(vacc_output, pop_df, by = "number", all.x = T)
          
          # Calculate proportion vaccinated on each day
          cols <- c("dose1","dose2","dose3")
          vacc_output[, (cols) := lapply(.SD, function(x) x/population), .SDcols = cols]
          
          # Convert to long format
          vacc_long <- melt(vacc_output, measure.vars = c("dose1", "dose2", "dose3"),
                            variable.name = "dose", value.name = "prop_dose")
          
          # Calculate vaccine coverage
          vacc_long[, cumu_dose := cumsum(prop_dose), by = .(number, dose)]
          
          setnames(vacc_long, c("date_week", "number"), c("date_week", "region_nb"))
        } else{
          vacc_output <- CJ(date_week = vacc_reg[,unique(date_week)], number = vacc_reg[,unique(number)], age1 = 0:90)
          age_groups_pop <- pop_df[,unique(age)]
          
          # Merge with population data
          min_ages_pop <- get_min_age(age_groups_pop)
          vacc_output[,age := cut(age1, c(min_ages_pop,Inf), labels = age_groups_pop, right = F)]
          vacc_output <- merge(vacc_output, pop_df[,!"region"], by = c("number","age"), all.x = T)
          vacc_output[, population := population/.N, by = .(number, age, date_week)]
          
          # Change age groups to vaccination age groups
          age_groups_vax <- vacc_reg[,unique(age)]
          min_ages_vax <- get_min_age(age_groups_vax)
          vacc_output[, age := cut(age1, c(min_ages_vax,Inf), labels = age_groups_vax, right = F)]
          
          # Merge with vaccination data
          vacc_output <- merge(vacc_output, vacc_reg[,!"total"], by = c("date_week","number","age"), all.x = T)
          
          cols <- c("dose1","dose2","dose3")
          
          vacc_output[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]
          vacc_output[, (cols) := lapply(.SD,function(x) x/.N), .SDcols = cols, by = .(date_week,number,age)]
          
          # Change age groups to model age groups
          min_ages <- get_min_age(age_groups)
          vacc_output[, age := cut(age1, c(min_ages,Inf), labels = age_groups, right = F)]
          # Sum doses over model age groups
          vacc_output <- vacc_output[,lapply(.SD, sum), .SDcols = c("population",cols), by = .(date_week,number,age)]
          
          # Calculate proportion vaccinated each day
          vacc_output[,(cols) := lapply(.SD, function(x) x/population), .SDcols = cols]
          
          # Melt to long format
          vacc_long <- melt(vacc_output, measure.vars = cols, variable.name = "dose", value.name = "prop_dose")
          # Calculate coverage
          vacc_long[, cumu_dose := cumsum(prop_dose), by = .(age, number, dose)]
          vacc_long[, population := NULL]
          setnames(vacc_long, c("date_week", "number"), c("date_week", "region_nb"))
        }
      }
    }
  }
  return(vacc_long)
}
