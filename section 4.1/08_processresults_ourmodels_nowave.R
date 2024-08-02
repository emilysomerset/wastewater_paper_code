# > sessionInfo()
# R version 4.4.1 (2024-06-14) -- "Race for Your Life"
# Platform: aarch64-apple-darwin20
# Running under: macOS Monterey 12.4

# To install OSplines package
# install.packages("remotes")
# remotes::install_github("AgueroZZ/OSplines")

library(OSplines) # OSplines_0.1.1
library(aghq) # aghq_0.4.1
library(dplyr) # dplyr_1.1.4
library(Matrix) # Matrix_1.7-0 
library(TMB) # TMB_1.9.14 

load("./data/work_d.RData")
source('process_results.R')

work_d <- work_d %>% 
  rename(y = value)

for (i in 1:103){
  load(paste0("./nowcast_nowave/model",i,'.RData'))
  
  ## Have to do a bit of extra work because I didn't save the true y values for the prediction stuff. 
  ## Need to fill them in so I can get quantile estimates. 
  df_full <- df_full %>% 
    dplyr::select(-'y') %>% # as mentioned above
    left_join(work_d %>% dplyr::select(y, sample_date, site_id), by = c("sample_date","site_id")) # as mentioned above
  
  
  results <-  process_results(df_full, 
                              tmbdat,
                              samps1,
                              polyOrder=3,  
                              id_group=1, 
                              id_group_name = NULL)
  
  write.csv(results$df_full, file = paste0('./nowcast_nowave/results_model',i,'.csv'))
  write.csv(results$station_ave_df, file = paste0('./nowcast_nowave/CANresults_model',i,'.csv'))
}
