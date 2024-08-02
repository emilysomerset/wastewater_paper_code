
# You first need to obtain a personal access token (PAT) by following the instruction via 
# this link: https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token.
# Enter the token as a string to the object PAT then run the following code to install the package:

# PAT <- " " # enter your PAT token here
# devtools::install_github("gqlNU/publicWW", auth_token = PAT)

## Will have problems downloading this with R version 4.4.1 since some dependencies 
## I first ran this code with a older version of R. 
## were removed from the CRAN. 

library(tibble) #tibble_3.2.1
library(dplyr) #dplyr_1.1.4
library(publicWW)    

#Functions
source("./data/make_daily_data.R")  # from publicWW modified by me

####################
###  user input
#####################


## All of this code is taken from https://github.com/gqlNU/publicWW
## Slight modification was made to get daily data instead of weekly data. 

#  region is included as fixed ('fixed') or random effect ('re') or not included at all ('none')
region_effect <- 'none'

#  LOD threshold
L <- 160

#  this is by default set to FALSE (excluding the genomic covariates) as we do not have permission to make the data public
include_genomic_covars <- FALSE

#  covariates to be included
covars <- tibble(covariates =  c('IMD_score','old_prop','young_prop','bame_proportion','population_density','f_industrial','pct_genome_coverage','num_snp'), #  covariates to be included
                 aggregation= c('wmean','wmean','wmean','wmean','wmean','wmean','none','none') # how to aggregate from LSOA to catchment, taking weighted mean or sum, none means that the covariate is spatio temporal and should note be aggregated.
)
if (!include_genomic_covars) {
  covars <- tibble(covariates =  c('IMD_score','old_prop','young_prop','bame_proportion','population_density','f_industrial'), #  covariates to be included
                   aggregation= c('wmean','wmean','wmean','wmean','wmean','wmean') # how to aggregate from LSOA to catchment, taking weighted mean or sum, none means that the covariate is spatio temporal and should note be aggregated.
  )
}

covars <- NULL
analysis_unit <- list(time='week',space='uwwCode')  #  space and time units of analysis (they should be columns in the dataset)

# Defaults to syears = NULL, smonths = NULL. To select specific time regions change these values.
syears  <- NULL #2021  #  which year of data to be modelled
smonths <- NULL #6:9

# construct catchment level covariates from LSOA values?
construct_catchment_covariates <- FALSE 
#####################
###  end user input
#####################

#  read wastewater data in
dat <- read_public_data('30Mar2022')
colnames(dat)[5:ncol(dat)] %>% length()
# 01/06/2021 to 30/03/2022

#  gather and format data for fitting
ww.weekly <- make_weekly_data(dat = dat, analysis_unit = analysis_unit, covars = covars, 
                              LOD = L, syears = syears, smonths = smonths,
                              construct_catchment_covariates=construct_catchment_covariates) # no covariates


ww.daily <- make_daily_data(dat = dat, analysis_unit = analysis_unit, covars = covars, 
                            LOD = L, syears = syears, smonths = smonths,
                            construct_catchment_covariates=construct_catchment_covariates)# no covariates


save(file=".data/work_d.RData", list=c("ww.daily","ww.weekly"))

