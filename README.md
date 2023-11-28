# Backend code for RShiny App: Analysis of COVID-19 outbreak risk at subnational level in the vaccine era
This repository contains the script and functions to fit data on COVID-19 outbreaks in various EU countries (currently France, Czechia, and Italy), generate 28-day case and death forecasts, and simulate the impact of changes in transmission on outbreak risks. The model is implemented using a distributed-lag version of the Endemic-Epidemic model (as described in the [surveillance](https://cran.r-project.org/web/packages/surveillance/index.html) and [hhh4addon](https://github.com/jbracher/hhh4addon) R packages). In France and Czechia, the model implemented is age-stratified (following the approach of the [hhh4contacts](https://cran.r-project.org/web/packages/hhh4contacts/index.html) package). The specifics of the model in each country are detailed in the file [backend_code_methods.pdf](backend_code_methods.pdf).

## Installation
Clone/download this project onto your machine.

The following R packages are required to run the code:
```
* data.table
* qs
* readxl
* surveillance
* remotes
* hhh4addon
* hhh4contacts
* sf
* socialmixr
* ggplot2
* spdep
* abind
* ISOweek
* dplyr
* tidyr
```
and can be installed in R by running:
```R
install.packages(c("data.table","qs","readxl",”surveillance”,”remotes”,”hhh4contacts”,”sf”,”socialmixr”,”ggplot2”,”spdep”,”abind”,”ISOweek”,”dplyr”,”tidyr”))
```
`hhh4addon` is installed when the code is run if it has not already been installed.

## Data
The analysis uses COVID-19 and demographic data from various sources. Some datasets, such as the case data, vaccination data and testing data, are by default downloaded in real-time as the code runs (see Table 1 in [backend_code_methods.pdf](backend_code_methods.pdf) for details of the sources for each country). Other datasets, such as the population data and contact data, are saved as static versions in the [Data](Data) directory. It contains the following files:

* [NUTS2021.xlsx](Data/NUTS2021.xlsx): File with urban/rural status of NUTS-3 regions from [Eurostat rural development methodology](https://ec.europa.eu/eurostat/web/rural-development/methodology)
* [demo_r_d2jan_1_Data.csv](Data/demo_r_d2jan_1_Data.csv): Population of EU countries in 2020 by year-of-age at different NUTS levels from [Eurostat demographic data](https://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=demo_r_d2jan&lang=en)
* [demo_r_pjangrp3_1_Data.csv](Data/demo_r_pjangrp3_1_Data.csv): Population of EU countries in 2020 in 5-year age groups at different NUTS levels downloaded from [Eurostat demographic data](https://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=demo_r_pjangrp3&lang=en)
* [index.csv](Data/index.csv): Index file of names and location codes of different places in [Google COVID-19 Open Data](https://github.com/GoogleCloudPlatform/covid-19-open-data)
* [pop_struct.csv](Data/pop_struct.csv): 2020 NUTS-3-level population data for France downloaded from [INSEE](https://www.insee.fr/fr/statistiques/1893198)
* [synthetic_contacts_2020.csv](Data/synthetic_contacts_2020.csv): Synthetic contact data for 172 countries from [Prem et al 2021](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009098) downloaded from [here](https://github.com/kieshaprem/synthetic-contact-matrices/tree/6e0eebc2524db277ae8110431b187ac589884dd9/output/syntheticcontactmatrices2020).
* [variant.csv](Data/variant.csv): National-level variant frequency data for EU countries downloaded from [Data on SARS-CoV-2 variants in the EU/EEA](https://www.ecdc.europa.eu/en/publications-data/data-virus-variants-covid-19-eueea). This dataset is a static backup version used only if there are errors in the version downloaded when the code runs.
* [File NUTS_RG_20M_2021_3035.shp](Data/NUTS_RG_20M_2021_3035.shp): contains all shapefiles, downloaded from [the Eurostat website](https://gisco-services.ec.europa.eu/distribution/v2/nuts/shp/NUTS_RG_20M_2021_3035.shp.zip), Copied in this folder because downloading the file in R can sometimes trigger errors.


## Running the code
To run the model and generate forecasts in the different countries: set the working directory to this repository, and run the command 
```R
source(“R/script_model/script_import_to_pred.R”)
```
The overall run time should be approximately 5 hours on a standard laptop with a 3.0 GHz processor and 32 GB RAM. This will generate the output files contained in the [Output](Output) directory. These files should then be moved to the [RShiny GitHub repository](https://github.com/alxsrobert/RShiny_covid_ECDC). Change the value of `pred_date` ([script_import_to_pred.R](R/script_model/script_import_to_pred.R)) to fit the model up to a more recent prediction date. 

The predictions on the [RShiny GitHub repository](https://github.com/alxsrobert/RShiny_covid_ECDC) are automatically updated weekly using the scripts in the [workflow](R/workflow) folder and [GitHub Actions](.github/workflows). The workflow scripts are similar to [script_import_to_pred.R](R/script_model/script_import_to_pred.R), but are separated by the type of forecasts ([workflow_fit](R/workflow/workflow_fit.R) to generate forecasts at previous and current dates, [workflow_sim](R/workflow/workflow_sim.R) to generate forecasts according to different scenarios). They were separated to avoid going over the GitHub Action time limit.

To run the calibration analysis over the last six months of data, run the command:
```R
source(“R/calibration/script_calib.R”)
```
This will generate various figures and tables summarising the calibration results, and comparing the predictions from our model to forecasts generated by the European COVID-19 Forecast Hub ensemble model. The figures and tables will be generated in different folders of the [Output](Output) directory using the full model, an empty age-stratifed, a non-age-stratified version and an empty non-stratified version of the model. This highlights the impact of added complexity on the forecasting ablities of the model.

The framework developed in this repository can be applied to countries that report daily subnational case data (at NUTS-3 level) and subnational death data (at NUTS1-NUTS2 level), along with local vaccine and testing data. If the data sources are stratified by age groups, then an age-stratified version of the model can be implemented.

## Modifying the code
Instructions for how to integrate a new country in the analysis, how to add a new variant to the fitting section, and how to add alternative NPIs in the Scenario forecasts can be found in the sections below. We also briefly present what is contained in every script.
Please contact us via email (<alexis.robert@lshtm.ac.uk>; <l.chapman4@lancaster.ac.uk>) or open an issue if you need assistance with editing the repository.

### How to add a new country to the backend code (Tag)
As most EU countries only have non-age-stratified subnational data available, we only provide detailed instructions for adding non-age-stratified forecasts here. All countries are identified by their two-letter ISO country code as defined in the `country_code` variable in the [index  file](Data/index.csv) of the [Google COVID-19 Open Data](https://github.com/GoogleCloudPlatform/covid-19-open-data), e.g. `FR` for France. 

The country-specific sections of code are indicated by a tag (`## TAG COUNTRY:`). Use `CTRL-SHIFT-F` (or `CMD-SHIFT-F` on a Mac) to find all occurrences of the tag in the repository. The description, expected input and output format of these sections is described in the comment lines following the tag.  

Modify the scripts and functions in the [R](R) directory as follows:
* [functions_import_data/function_import_data.R](R/functions_import_data/function_import_data.R)
  * `import_case()`: 
    * Process the data so that it contains the following variables, by adding `else if` statements where `case` and `case1` are defined:
      * `location_key`: NUTS-3 region code
      * `date`: date in format YYYY-MM-DD
      * `new_confirmed`: daily number of new confirmed cases
      * `new_tested`: daily number of new individuals tested; this should be `NA` if the testing data is missing from the case data source
      * `cumulative_confirmed`: daily cumulative number of confirmed cases
      * `cumulative_tested`: daily cumulative number of individuals tested; this should be `NA` if the testing data is missing from the case data source
      * `week_day`: year and week number in the format YYYY-WW where WW is the week of the year as a decimal number (00-53), with the first Monday of the year as day 1 of week 1.
  * `import_death()`: 
    * Import and process the data on the number of death per region so that it contains the following variables:
      * `location_key`: NUTS-1/NUTS-2 region code
      * `date`: date in format YYYY-MM-DD
      * `new_deceased`: daily number of new deaths
  * `import_vacc()`:
    * Import the vaccine data from a country-specific URL (containing the number  of daily cases, per 10 year age band), the output of this `if()` section should be a data frame called "vacc_regions", which contains the following variables:
      * `date_week`: date in format YYYY-MM-DD.
      * `number`: Region code.
      * `Population`: Number of inhabitants in this region + age group (If not reported, can be set to NA, for now it is only used in France).
      * `Age`: Age group described in this row.
      * `dose1`: Number of 1st dose distributed on that day / region / age group.
      * `dose2`: Number of 2nd dose distributed on that day / region / age group.
      * `dose3`: Number of booster dose distributed on that day / region / age group.
  * `import_pop()`:
    * If the population data source is country-specific, add an else if statement after `if (country == “FR”)` (i.e. after the Tag).
    * If the data source is Eurostat, add the two-letter country code to the vector of country codes  (in `else if (country %in% c(“CZ”))`).
    * Process the data so that it has the following variables:
      * `number`: part or whole NUTS-3 region code, in a format that can be matched with `location_key` in the `case` data table, via the `get_reg_nb()` function in `function_utils.R`
      * `population`: population of NUTS-3 region
    * If `total == FALSE`, the data frame returned by the function should also contain the following variables (in addition to the two listed above):
      * `region`: region name
      * `age`: age group
  * `import_test()`:
    * If the testing data is not already present in the case data table, import it via the `import_test()` function.
    * If the data source for the testing data is country-specific, add an `if` statement in `import_test()` to import and process the testing data, otherwise the national-level total testing data from the [ECDC testing database](https://www.ecdc.europa.eu/en/publications-data/covid-19-testing) is used by default. The dataframe created in this if section should contain the following variables
      * `date`: date in format YYYY-MM-DD
      * `age`: age group (10 year age band)
      * `population`: Number of inhabitants
      * `nb_tests`: Number of tests on that date / age group
  * `import_contact()`:
    * Import a contact matrix from a country-specific survey (like in France), or a multi-country analysis (like in Slovakia).
    * The `if()` section must define a matrix named `C` (dimensions 9*9), which contains the number of contacts between the different age groups (0-10 years old; 10-20;...; 80+).
  * `import_map()`:
    * Create column `reg_nb` in map1, which corresponds to the reference of each subnational area in the `index` database (created in `import_index()`). To do so, first remove overseas territories (if any), select the entries of interest in `index`, set a key in `index`, and add the subregion code to `map1` (in a column called `reg_nb`).
  * `import_all_files()`:
    * If the population data is taken from the google database, set `google <- T`.
    * If vaccine data is taken from the [ECDC vaccine database](https://www.ecdc.europa.eu/en/publications-data/data-covid-19-vaccination-eu-eea), set `vacc_ecdc <- T`.
    * If the testing data is included in the case database (imported with `import_case()`), set `test = NULL` (e.g. France).
    * Otherwise, if the testing data comes from a national database (age-stratified), use `import_test()` with `total = F` (e.g. Czechia).
    * Otherwise, to use the ECDC testing database, run `import_test()` with `total = T` (e.g. Italy)
* [functions_model_sim/function_covariates.R](R/functions_model_sim/function_covariates.R):
  * `dow_cov()`:
    * Define the vector `bank_holiday`, which contains the dates of all bank holidays in the country in 2021 and 2022.
  * `all_covariates()`:
    * When computing the proportion of the population tested in the past two weeks (calling `test_by_age_cov`), if the testing data is included in the incidence dataset, use `dt_incidence` as an argument when calling the function `test_by_age_cov`, otherwise use `dt_test`.
* [functions_import_data/function_adjust_age_group.R](R/functions_import_data/function_adjust_age_group.R):
  * `adjust_age_group()`:
    * Create a data table named `vacc_long`, which contains the proportion of the population 
per region / age group who received each dose of vaccine (on that day and cumulative). `vacc_long` contains the following variables:
      * `date_week`: date in format YYYY-MM-DD
      * `region_nb`: Region code
      * `age`: Age group described in this row
      * `dose`: Dose number
      * `prop_dose`: Proportion of the population vaccinated that day
      * `cumu_dose`: Cumulative proportion of the population vaccinated
* [function_utils.R](R/function_utils.R):
  * Add an `if()` section to functions `get_reg_nb`, `get_reg_nb`, `get_age`, `get_nuts2_reg`, and `get_nuts3_reg`. These sections are used to extract the region number, region-age, age, NUTS-2 code, and NUTS-3 code from Google COVID-19 Open Data `location_key`.
* [functions_model_sim/function_calculate_cfr.R](R/functions_model_sim/function_calculate_cfr.R):
  * `calculate_cfr()`: If the spatial granularity of case and death data differs, aggregate the case data (object `case`) at NUTS-2 or NUTS-1 level (example for Italy and France).


### How to add a new variant
A new variant, `X` say, can be added to the model by:
* Adding a line in the `import_variant()` function in [function_import_data.R](R/functions_import_data/function_import_data.R) to set which PANGO variant codes in the ECDC variant data should be coded as `variant = X` in the `variant` data table.
* Adding a binary matrix `X_mat` in the `variant_covariate()` function in [function_covariates.R](R/functions_model_sim/function_covariates.R) with 1s when `X` was above a certain proportion of cases, `prop_X`, and adding `prop_X` as an argument of `variant_covariate()` with a default value.
* `X_mat` should be a matrix with the same dimensions as `observed` in the `sts` object with identical columns for the different region(-age) groups (as the variant data is only available at a national level).
* Adding `X` to the epidemic terms `ne_terms` in the model equations (in [function_generate_list_output.R](R/script_model/function_generate_list_output.R)).

### How to add a new targeted age group

New scenarios can be generated by adding new elements to the vector `target` in the script [function_generate_list_output.R](`R/script_model/function_generate_list_output.R`). The age group(s) included in this new element should be defined in the following `for` loop, using the vector `cols` (similar to the other elements of the `target` vector). For instance: 

```R
if (target[l] == "children") {
cols <- c(grep("0-9", colnames(model_fit$stsObj@observed)), grep("10-19", colnames(model_fit$stsObj@observed)))
}
```

Bear in mind that adding scenarios will make the output dataset larger, eventually slowing down the first display of the figures in the Shiny App.

### Presentation of each script
Directory structure:


```
├── .github
│   ├── workflows
│   │   ├── copy_shiny.yml
│   │   └── run.yml
├── Data
│   ├── NUTS2021.xlsx
│   ├── NUTS_RG_20M_2021_3035.shp.zip
│   ├── demo_r_d2jan_1_Data.csv
│   ├── demo_r_pjangrp3_1_Data.csv
│   ├── index.csv
│   ├── pop_struct.csv
│   ├── source_NUTS.txt
│   ├── synthetic_contacts_2020.csv
│   └── variant.csv
├── LICENSE
├── Output
│   ├── output_model_CZ.RDS
│   ├── output_model_FR.RDS
│   └── output_model_IT_total.RDS
├── R
│   ├── calibration
│   │   ├── figures_calib.R
│   │   ├── function_main_figure.R
│   │   ├── script_calib.R
│   │   └── function_calib.R
│   ├── function_utils.R
│   ├── functions_custom_hhh4
│   │   ├── function_interpretControl.R
│   │   └── function_simulate.R
│   ├── functions_import_data
│   │   ├── function_adjust_age_group.R
│   │   ├── function_cumulative_incidence.R
│   │   ├── function_import_data.R
│   │   ├── import_libraries.R
│   │   └── import_scripts.R
│   ├── functions_model_sim
│   │   ├── function_calculate_cfr.R
│   │   ├── function_covariates.R
│   │   ├── function_create_sts.R
│   │   ├── function_format_pred.R
│   │   ├── function_models.R
│   │   ├── function_n_step_ahead.R
│   │   └── function_predict_covariates.R
│   ├── script_model
│   │   ├── function_generate_list_output.R
│   │   └── script_import_to_pred.R
│   └── workflow
│       ├── workflow_fit.R
│       └── workflow_sim.R
├── README.md
└── backend_code_methods.pdf
```
The scripts in the [R](R) folder contain the following:
* `function_utils.R`: General utility functions used by multiple functions/scripts, e.g. functions for extracting the region numbers from the Google `location_key`
* `functions_custom_hhh4/`
  * `function_interpretControl.R`: Custom version of `interpretControl()` function from [hhh4addon](https://github.com/jbracher/hhh4addon) package that enables simulation of age-structured distributed-lag models with separate autoregressive and neighbourhood components
  * `function_simulate.R`: Custom version of `simulate.hhh4lag()` function from [hhh4addon](https://github.com/jbracher/hhh4addon) package to enable simulation of distributed-lag age-stratified models with separate autoregressive and neighbourhood components and to account for cumulative incidence being a covariate in autoregressive and neighbourhood predictors
* `functions_import_data/`
  * `function_adjust_age_group.R`: Function for calculating vaccination coverage in age groups in case data and population data
  * `function_cumulative_incidence.R`: Function for calculating the cumulative incidence per age group (if model is age-stratified), per region, during each variant wave (wild-type + Alpha, Delta, Omicron)
  * `function_import_data.R`: Functions for importing the datasets required to fit the model
  * `import_libraries.R`: Script that loads the packages required for the analysis
  * `import_scripts.R`: Script that sources all the files with functions required for the analysis
* `functions_model_sim/`
  * `function_calculate_cfr.R`: Functions to compute Case Fatality Ratios from recent weeks ; Create a data frame containing the number of weekly deaths per region ; Implement regression models to predict changes in Case fatality Ratios ; and generate up to four-week-ahead forecasts on the number of weekly deaths.
  * `function_covariates.R`: Functions for constructing the covariate matrices used in fitting the model
  * `function_create_sts.R`: Function for creating the `sts` object used in fitting the model that accounts for whether the model is age-structured or not, and creates the `neighbourhood` entry of the `sts` object accordingly. Adapted from the `noroBE()` function in [hhh4contacts](https://cran.r-project.org/web/packages/hhh4contacts/index.html)
  * `function_format_pred.R`: Function that generates a data frame of the forecasts in the `list_pred` object output from `nStepAhead()`
  * `function_models.R`: Function for fitting the (age-stratified) distributed-lag `hhh4` model that creates the `control` object from the model equations and covariate data
  * `function_n_step_ahead.R`: Function for forecasting cases `n` days ahead with predicted values of covariates and different scenarios for variant transmissibility and importation level
  * `function_predict_covariates.R`: Functions for projecting covariates for case forecasts with different projection types
* `calibration/`
  * `function_calib.R`: Functions to run the calibration analysis for some (or all) countries at a set of dates, and generate figures and tables describing the calibration results.
  * `figures_calib.R`: Functions to produce the calibration figures and tables per country. 
  * `function_main_figure.R`: Functions to produce the calibration figures and tables comparing scores and performance in each country. 
  * `script_calib.R`: Top-level script to run the calibration analysis and generate the figures and tables, for a set of models, dates and countries. The calibration files, figures, and tables will be added in the [Output](Output) directory. 
* `workflow/`
  * `workflow_fit.R`: Top-level script used in GitHub actions to produce case forecasts at current and previous dates.
  * `workflow_sim.R`: Top-level script used in GitHub actions to produce case forecasts in different transmission scenarios.
* `script_model/`
  * `function_generate_list_output.R`: Function for generating case forecasts for a specific country over a certain date range
  * `script_import_to_pred.R`: Top-level script for setting options and running code to produce case forecasts

## Authors
* Alexis Robert: <alexis.robert@lshtm.ac.uk>
* Lloyd Chapman: <l.chapman4@lancaster.ac.uk>
