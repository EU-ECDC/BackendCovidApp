#### Import hhh4 functions / data ####

if (!is.element("Data", dir())) stop("Set working directory to Backend_covid_ECDC")

# Increase the timeout time for the download for slower connections
options(timeout = max(600, getOption("timeout")))

# Import all libraries
source("R/functions_import_data/import_libraries.R")
# Import all scripts
source("R/functions_import_data/import_scripts.R")
# Import functions for calibration
source("R/calibration/function_calib.R")
# Import functions for figures
source("R/calibration/figures_calib.R")
# Import functions for main figures
source("R/calibration/function_main_figure.R")

# Set the interpretControl function to our customised version
assignInNamespace("interpretControl", my.interpretControl, ns = "hhh4addon", 
                  pos = "package:hhh4addon")
# Set the simulate function to our customised version
assignInNamespace("simulate.hhh4lag", my.simulate.hhh4lag, ns = "hhh4addon", 
                  pos = "package:hhh4addon")

country <- "all"
# country <- "CZ"
# If run == TRUE: Run the calibration models; otherwise: just generate the figures
run <- T #FALSE

# Define the prediction date
pred_date <- as.character(Sys.Date() - 3)

# Run calibration:
# Full age-stratified model
run_calib(country, all_total = F, empty = F, run, pred_date)
# Empty age-stratified model (without covariates or seasonality)
run_calib(country, all_total = F, empty = T, run, pred_date)
# Full non-age-stratified model
run_calib(country, all_total = T, empty = F, run, pred_date)
# Empty non-age-stratified model (without covariates or seasonality)
run_calib(country, all_total = T, empty = T, run, pred_date)

# Make tables of prediction scores for empty model (no covariates or seasonality) and full model
table_scores(all_total = F, "case")
table_scores(all_total = F, "death")
table_scores(all_total = T, "case")
table_scores(all_total = T, "death")

