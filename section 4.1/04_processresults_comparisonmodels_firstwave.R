# > sessionInfo()
# R version 4.4.1 (2024-06-14) -- "Race for Your Life"
# Platform: aarch64-apple-darwin20
# Running under: macOS Monterey 12.4

library(dplyr) # dplyr_1.1.4
library(Matrix) # Matrix_1.7-0 


for (i in 1:54){
  load(paste0('./nowcast_firstwave_xiaotian/model',i,'.RData'))
  
  df = full_join(res_df %>% mutate(type = "train"), 
                 pred_df %>% dplyr::select(-'time.ahead') %>% mutate(type = "test"))
  
  df = df %>% 
    group_by(Location, date) %>% 
    mutate(vv = 1:2501) %>% # change this to be the number of iterations - burnin + 1
    group_by(Location,vv) %>% 
    arrange(date) %>% 
    mutate(diff1 = c(NA,diff(value, lag = 1)),
           diff2 = c(rep(NA,2), diff(value, lag = 2)),
           diff7 = c(rep(NA,7), diff(value, lag = 7))) 
  
  data_full = full_join(trainingdata, testingdata) %>% 
    rename(value.raw = value)
  
  df <- df %>% 
    left_join(data_full, by = c("date","Location"))
  
  df_summary = df %>% 
    group_by(Location, date, value.raw) %>% 
    summarise(f1_IS_fixed = median(value),
              f1_IS_fixed_upr = quantile(value, 0.975),
              f1_IS_fixed_lwr = quantile(value, 0.025),
              exp_f1_IS_fixed = median(exp(value)),
              exp_f1_IS_fixed_upr = quantile(exp(value), 0.975),
              exp_f1_IS_fixed_lwr = quantile(exp(value), 0.025),
              pred_mm = median(pred),
              pred_lwr = quantile(pred, 0.975),
              pred_upr = quantile(pred, 0.025),
              mu_total = median(mu_total),
              p1_increase = length(which(diff1>0))/length(na.omit(diff1)),
              p2_increase = length(which(diff2>0))/length(na.omit(diff2)),
              p7_increase = length(which(diff7>0))/length(na.omit(diff7)),
              quantile_obs = ifelse(!is.na(value.raw[1]), length(which(pred <= value.raw[1]))/length(pred), NA))
  
  # There is a lot here that was not used for the final manuscript
  df_summary_canada = df %>% 
    group_by(date, vv, mu_total) %>%
    summarise(mean_log = mean(value),
              mean_exp = mean(exp(value)),
              mean_exp_corr = mean(exp(value+sigma2/2)),
              pred = mean(pred)) %>% 
    group_by(vv) %>% 
    arrange(date) %>% 
    mutate(diff1 = c(NA,diff(mean_log, lag = 1)),
           diff2 = c(rep(NA,2), diff(mean_log, lag = 2)),
           diff7 = c(rep(NA,7), diff(mean_log, lag = 7))) %>% 
    mutate(diff1_exp = c(NA,diff(mean_exp, lag = 1)),
           diff2_exp = c(rep(NA,2), diff(mean_exp, lag = 2)),
           diff7_exp = c(rep(NA,7), diff(mean_exp, lag = 7))) %>% 
    mutate(diff1_exp_corr = c(NA,diff(mean_exp_corr, lag = 1)),
           diff2_exp_corr = c(rep(NA,2), diff(mean_exp_corr, lag = 2)),
           diff7_exp_corr = c(rep(NA,7), diff(mean_exp_corr, lag = 7))) %>% 
    group_by(date) %>% 
    summarise(ave_f1_IS_fixed = median(mean_log),
              ave_f1_IS_fixed_upr = quantile(mean_log, 0.975),
              ave_f1_IS_fixed_lwr = quantile(mean_log, 0.025),
              ave_exp_f1_IS_fixed = median(mean_exp),
              ave_exp_f1_IS_fixed_upr = quantile(mean_exp, 0.975),
              ave_exp_f1_IS_fixed_lwr = quantile(mean_exp, 0.025),
              ave_exp_corr_f1_IS_fixed = median(mean_exp_corr),
              ave_exp_corr_f1_IS_fixed_upr = quantile(mean_exp_corr, 0.975),
              ave_exp_corr_f1_IS_fixed_lwr = quantile(mean_exp_corr, 0.025),
              pred_mm = median(pred),
              pred_lwr = quantile(pred, 0.975),
              pred_upr = quantile(pred, 0.025),
              mu_total = median(mu_total),
              p1_increase = length(which(diff1>0))/length(na.omit(diff1)),
              p2_increase = length(which(diff2>0))/length(na.omit(diff2)),
              p7_increase = length(which(diff7>0))/length(na.omit(diff7)),
              p1_exp_increase = length(which(diff1_exp>0))/length(na.omit(diff1_exp)),
              p2_exp_increase = length(which(diff2_exp>0))/length(na.omit(diff2_exp)),
              p7_exp_increase = length(which(diff7_exp>0))/length(na.omit(diff7_exp)),
              p1_exp_corr_increase = length(which(diff1_exp_corr>0))/length(na.omit(diff1_exp_corr)),
              p2_exp_corr_increase = length(which(diff2_exp_corr>0))/length(na.omit(diff2_exp_corr)),
              p7_exp_corr_increase = length(which(diff7_exp_corr>0))/length(na.omit(diff7_exp_corr)))
  

# this was done to save space on my computer after I got what I needed from the model.  
  if (file.exists(paste0('./nowcast_firstwave_xiaotian/model',i,'.RData'))) {
    #Delete file if it exists
    file.remove(paste0('./nowcast_firstwave_xiaotian/model',i,'.RData'))
  }
  
  write.csv(df_summary, file = paste0('./nowcast_firstwave_xiaotian/summary_model',i,'.csv'))
  write.csv(df_summary_canada, file = paste0('./nowcast_firstwave_xiaotian/CAN_summary_model',i,'.csv'))
}


## We also want to process the full model in the same way: 

load(paste0('./nowcast_firstwave_xiaotian/fullmodel.RData'))

df = res_df 

df = df %>% 
  group_by(Location, date) %>% 
  mutate(vv = 1:2501) %>% # change this to be the number of iterations - burnin + 1
  group_by(Location,vv) %>% 
  arrange(date) %>% 
  mutate(diff1 = c(NA,diff(value, lag = 1)),
         diff2 = c(rep(NA,2), diff(value, lag = 2)),
         diff7 = c(rep(NA,7), diff(value, lag = 7))) 

data_full <- data_full %>% 
  rename(value.raw = value)

df <- df %>% 
  left_join(data_full, by = c("date","Location"))

df_summary = df %>% 
  group_by(Location, date, value.raw) %>% 
  summarise(f1_IS_fixed = median(value),
            f1_IS_fixed_upr = quantile(value, 0.975),
            f1_IS_fixed_lwr = quantile(value, 0.025),
            exp_f1_IS_fixed = median(exp(value)),
            exp_f1_IS_fixed_upr = quantile(exp(value), 0.975),
            exp_f1_IS_fixed_lwr = quantile(exp(value), 0.025),
            pred_mm = median(pred),
            pred_lwr = quantile(pred, 0.975),
            pred_upr = quantile(pred, 0.025),
            mu_total = median(mu_total),
            p1_increase = length(which(diff1>0))/length(na.omit(diff1)),
            p2_increase = length(which(diff2>0))/length(na.omit(diff2)),
            p7_increase = length(which(diff7>0))/length(na.omit(diff7)),
            quantile_obs = ifelse(!is.na(value.raw[1]), length(which(pred <= value.raw[1]))/length(pred), NA))

df_summary_canada = df %>% 
  group_by(date, vv, mu_total) %>%
  summarise(mean_log = mean(value),
            mean_exp = mean(exp(value)),
            mean_exp_corr = mean(exp(value+sigma2/2)),
            pred = mean(pred)) %>% 
  group_by(vv) %>% 
  arrange(date) %>% 
  mutate(diff1 = c(NA,diff(mean_log, lag = 1)),
         diff2 = c(rep(NA,2), diff(mean_log, lag = 2)),
         diff7 = c(rep(NA,7), diff(mean_log, lag = 7))) %>% 
  mutate(diff1_exp = c(NA,diff(mean_exp, lag = 1)),
         diff2_exp = c(rep(NA,2), diff(mean_exp, lag = 2)),
         diff7_exp = c(rep(NA,7), diff(mean_exp, lag = 7))) %>% 
  mutate(diff1_exp_corr = c(NA,diff(mean_exp_corr, lag = 1)),
         diff2_exp_corr = c(rep(NA,2), diff(mean_exp_corr, lag = 2)),
         diff7_exp_corr = c(rep(NA,7), diff(mean_exp_corr, lag = 7))) %>% 
  group_by(date) %>% 
  summarise(ave_f1_IS_fixed = median(mean_log),
            ave_f1_IS_fixed_upr = quantile(mean_log, 0.975),
            ave_f1_IS_fixed_lwr = quantile(mean_log, 0.025),
            ave_exp_f1_IS_fixed = median(mean_exp),
            ave_exp_f1_IS_fixed_upr = quantile(mean_exp, 0.975),
            ave_exp_f1_IS_fixed_lwr = quantile(mean_exp, 0.025),
            ave_exp_corr_f1_IS_fixed = median(mean_exp_corr),
            ave_exp_corr_f1_IS_fixed_upr = quantile(mean_exp_corr, 0.975),
            ave_exp_corr_f1_IS_fixed_lwr = quantile(mean_exp_corr, 0.025),
            pred_mm = median(pred),
            pred_lwr = quantile(pred, 0.975),
            pred_upr = quantile(pred, 0.025),
            mu_total = median(mu_total),
            p1_increase = length(which(diff1>0))/length(na.omit(diff1)),
            p2_increase = length(which(diff2>0))/length(na.omit(diff2)),
            p7_increase = length(which(diff7>0))/length(na.omit(diff7)),
            p1_exp_increase = length(which(diff1_exp>0))/length(na.omit(diff1_exp)),
            p2_exp_increase = length(which(diff2_exp>0))/length(na.omit(diff2_exp)),
            p7_exp_increase = length(which(diff7_exp>0))/length(na.omit(diff7_exp)),
            p1_exp_corr_increase = length(which(diff1_exp_corr>0))/length(na.omit(diff1_exp_corr)),
            p2_exp_corr_increase = length(which(diff2_exp_corr>0))/length(na.omit(diff2_exp_corr)),
            p7_exp_corr_increase = length(which(diff7_exp_corr>0))/length(na.omit(diff7_exp_corr)))


write.csv(df_summary, file = paste0('./nowcast_firstwave_xiaotian/summary_fullmodel.csv'))
write.csv(df_summary_canada, file = paste0('./nowcast_firstwave_xiaotian/CAN_summary_fullmodel.csv'))


