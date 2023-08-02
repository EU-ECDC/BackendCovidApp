if (!is.element("Data", dir())) stop("Set working directory to Backend_covid_ECDC")

# Increase the timeout time for the download for slower connections
options(timeout = max(600, getOption("timeout")))

# Import all libraries
source("R/functions_import_data/import_libraries.R")
# Import all scripts
source("R/functions_import_data/import_scripts.R")
# Set the interpretControl function to our customised version
assignInNamespace("interpretControl", my.interpretControl, ns = "hhh4addon", 
                  pos = "package:hhh4addon")
# Set the simulate function to our customised version
assignInNamespace("simulate.hhh4lag", my.simulate.hhh4lag, ns = "hhh4addon", 
                  pos = "package:hhh4addon")

download <- T
data_file <- c("")
nsim <- 100

### Generate output file for each country
## TAG COUNTRY
## Add new countries to the loop (using the 2 letter ISO code), and set whether to use 
## total = T (age-stratified) or total = F (non-age-stratified)
for(country in c("IT", "CZ", "FR")){
  print( paste0( "Running ", country, " ... " ) )
  if (country == "IT") total <- T else total <- F
  # Define the prediction date
  pred_date <- Sys.Date()
  range_dates <- c("2020-09-01", pred_date)
  
  # Generate list_output object
  list_output <- generate_list_output(country, range_dates, download, data_file, 
                                      total, last_date = TRUE, nsim = nsim)
  
  # Save list_output
  if (!total){
    saveRDS(list_output, paste0("Output/output_model_", country, ".RDS"))
  } else {
    list_output$pop_age <- list_output$pop * 0
    list_output$pop_age[] <- 1
    saveRDS(list_output, paste0("Output/output_model_", country, "_total.RDS"))
  }
}
