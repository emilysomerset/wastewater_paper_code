# R version 4.4.1 (2024-06-14)
# Platform: aarch64-apple-darwin20
# Running under: macOS Monterey 12.4

library(dplyr)      # dplyr_1.1.4
library(data.table) # data.table_1.15.4
library(refund)     # refund_0.1-35 
library(rstan)      # rstan_2.32.6
library(lubridate)  # lubridate_1.9.3
library(forecast)   # forecast_8.23.0

source("WWmodel_meanconcentrate.R")
source("WW_forecast.R")

load('./data/work_d.RData')

modeldata2 = work_d

# just renaming some stuff so I can use their functions
modeldata2 <- modeldata2 %>% 
  rename("date" = sample_date,
         "Location" = site_id) %>% 
  mutate(log.value.raw = log(value))


date_range <-modeldata2$date %>% range()

modeldata2 <- data.frame(date = rep(seq(date_range[1], date_range[2],1),each=15),
                         Location = rep(unique(modeldata2$Location),length(seq(date_range[1], date_range[2],1)) )) %>% 
  left_join(data.frame(modeldata2), by = c("Location","date")) %>% 
  mutate(target = "N2") %>% 
  as.data.table(.)

dates_with_newdata <- modeldata2 %>% 
  filter(date >= ymd('2022-12-01') & date <= ymd('2023-03-31')) %>% 
  filter(!is.na(value))%$% date %>% unique()


dates_totry = dates_with_newdata


for (i in 1:length(dates_totry)) {
  df <- modeldata2 %>% 
    filter(date <= dates_totry[i])
  
  
  trainingdata= df 
  testingdata <- modeldata2 %>% 
    filter(date > dates_totry[i] & date <= (dates_totry[i] + 7))
  
  set.seed(2)

  res <- WWmodel_meanconcentrate(modeldata = trainingdata, 
                                 ID = c("Location", "target"), 
                                 date = "date", 
                                 value = "log.value.raw", 
                                 covariate = NULL, 
                                 iteration = 5000, 
                                 burnin = 2500, 
                                 cores = 4,
                                 chains=4,
                                 nbasis = 50)
  
  set.seed(3)
  
  pred_df = fcst_trend(h.ahead=7, 
                       modeldata = trainingdata, 
                       model_res=res, 
                       ID=c("Location", "target"),
                       date ="date", 
                       value = "log.value.raw", 
                       covariate = NULL, 
                       iteration = 5000, 
                       burnin = 2500,
                       nbasis = 50,
                       fpca.fit = res$fpca.fit)
  
  res_df <- res$df
  
  # this is how I saved the results. 
  # save(file=paste0("./nowcast_nowave_morebasis_xiaotian/model",i,".RData"), list = c("trainingdata","testingdata", "res_df","pred_df"))
}






