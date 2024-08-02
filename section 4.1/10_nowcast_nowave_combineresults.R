# > sessionInfo()
# R version 4.4.1 (2024-06-14) -- "Race for Your Life"
# Platform: aarch64-apple-darwin20
# Running under: macOS Monterey 12.4

library(dplyr) # dplyr_1.1.4
library(Matrix) # Matrix_1.7-0 

#####################################################
# Combine the results for our model
#####################################################

df_full <- read.csv('./nowcast_firstwave/results_fullmodel.csv')
df_full_can <- read.csv('./nowcast_firstwave/CANresults_fullmodel.csv')
df_full <- df_full %>% 
  mutate(sample_date = ymd(sample_date))
df_full_can <- df_full_can %>% 
  mutate(sample_date = ymd(sample_date)) %>% 
  mutate(model = "fullmodel")

dates_with_newdata <- df_full %>% 
  filter(sample_date >= ymd('2022-12-01') & sample_date <= ymd('2023-03-31')) %>% 
  filter(!is.na(y))%$% sample_date %>% unique()

dates_totry = dates_with_newdata


df_full <- df_full %>% dplyr::select(- 'X')
df_full_can <- df_full_can %>% dplyr::select(- 'X')

for (i in 1:length(dates_totry)){
tmp <- read.csv(paste0('./nowcast_nowave/results_model',i,'.csv'))
tmp2 <- read.csv(paste0('./nowcast_nowave/CANresults_model',i,'.csv'))

tmp <- tmp %>% 
  mutate(sample_date = ymd(sample_date)) %>% 
  filter(sample_date <= dates_totry[i]) %>% 
  select(-'X')

tmp2 <- tmp2 %>% 
  mutate(sample_date = ymd(sample_date)) %>% 
  filter(sample_date <= dates_totry[i]) %>% 
  select(-'X') %>% 
  mutate(model = paste0("model",i))

df_full <- full_join(df_full, tmp)
df_full_can <- full_join(df_full_can, tmp2)
}

df_full_nowave <- df_full
df_full_can_nowave <- df_full_can

save(list='df_full_nowave', file=paste0('./nowcast_nowave/results_allmodels.RData'))
save(list='df_full_can_nowave', file=paste0('./nowcast_nowave/CANresults_allmodels.RData'))

rm(list=ls()[!(ls()%in% 'dates_totry')])





#####################################################
# Combine the results for the comparison model. 
#####################################################

df_full_xiaotian <- read.csv('./nowcast_firstwave_xiaotian/summary_fullmodel.csv')
df_full_can_xiaotian <- read.csv('./nowcast_firstwave_xiaotian/CAN_summary_fullmodel.csv')
df_full_xiaotian <- df_full_xiaotian %>% 
  rename(y = value.raw,
         smoothed = mu_total,
         sample_date = date) %>% 
  mutate(sample_date = ymd(sample_date),
         exp_smoothed = exp(smoothed)) %>% 
  dplyr::select(-'X') %>% 
  mutate(model = "fullmodel") 

df_full_can_xiaotian <- df_full_can_xiaotian %>% 
  rename(smoothed = mu_total,
         sample_date = date) %>% 
  mutate(sample_date = ymd(sample_date),
         exp_smoothed = exp(smoothed)) %>% 
  dplyr::select(-'X') %>% 
  mutate(model = "fullmodel")  

for (i in 1:length(dates_totry)){
  tmp <- read.csv(paste0('./nowcast_nowave_morebasis_xiaotian/summary_model',i,'.csv'))
  tmp2 <- read.csv(paste0('./nowcast_nowave_morebasis_xiaotian/CAN_summary_model',i,'.csv'))
  
  tmp <- tmp %>% 
    mutate(date = ymd(date)) %>% 
    rename(sample_date = date,
           y = value.raw,
           smoothed = mu_total) %>% 
    mutate(model = paste0("model",i)) %>% 
    mutate(exp_smoothed = exp(smoothed)) %>% 
    dplyr::select(-'X') %>% 
    filter(sample_date <= dates_totry[i]) 
  
  tmp2 <- tmp2 %>% 
    mutate(date = ymd(date)) %>% 
    rename(sample_date = date,
           smoothed = mu_total) %>% 
    mutate(model = paste0("model",i)) %>% 
    mutate(exp_smoothed = exp(smoothed)) %>% 
    dplyr::select(-'X')%>% 
    filter(sample_date <= dates_totry[i]) 
  
  df_full_xiaotian <- rbind(df_full_xiaotian, tmp) 
  df_full_can_xiaotian <- rbind(df_full_can_xiaotian, tmp2) 
  }

df_full_nowave_xiaotian <- df_full_xiaotian
df_full_can_nowave_xiaotian <- df_full_can_xiaotian

save(list='df_full_nowave_xiaotian', file=paste0('./nowcast_nowave_morebasis_xiaotian/results_allmodels_xiaotian.RData'))
save(list = 'df_full_can_nowave_xiaotian', file=paste0('./nowcast_nowave_morebasis_xiaotian/CANresults_allmodels_xiaotian.RData'))


