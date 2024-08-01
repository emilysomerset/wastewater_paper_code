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

## Have to do a bit of extra work because I didn't save the true y values for the prediction stuff. 
## Need to fill them in so I can get quantile estimates. 
load("./data/work_d.RData")

work_d <- work_d %>% 
  rename(y = value)

for (i in 1:54){
  load(paste0("./nowcast_firstwave/model",i,'.RData'))
  
  df_full <- df_full %>% 
    dplyr::select(-'y') %>% # as mentioned above
    left_join(work_d %>% dplyr::select(y, sample_date, site_id), by = c("sample_date","site_id")) # as mentioned above
  
  
  ## The remainder of this code is not as neat as done in Section 3.1 and 3.2 
  ## Could wrap all of this into a function like I did for those sections
  ## There's also not necessary computations that were never used for the paper. 
  
  P = as.matrix(tmbdat$P)
  X = as.matrix(tmbdat$X)
  Xfstat = as.matrix(tmbdat$Xfstat)
  daily = as.matrix(tmbdat$daily)
  obs = as.matrix(tmbdat$obs)
  knots = tmbdat$knots
  
  coefsamps1 <- samps1$samps[1:ncol(P),]
  global_samps1 <- samps1$samps[(ncol(P) + 1):(ncol(P) + ncol(X)-ncol(Xfstat)),]
  Xf_samps1 <- samps1$samps[(ncol(P) + ncol(X)-ncol(Xfstat) + 1):(ncol(P) + ncol(X)),]
  daily_samps1 <- samps1$samps[(ncol(P) + ncol(X) + 1):(ncol(P) + ncol(X) + ncol(daily)),]
  IS_samps1 <- samps1$samps[(ncol(P) + ncol(X) + ncol(daily)+1):(ncol(P) + ncol(X) + ncol(daily)+ncol(obs)),]
  IS_deriv_samps1 <- samps1$samps[(ncol(P) + ncol(X) + ncol(daily)+ncol(obs)+1):(nrow(samps1$samps)),]
  
  f1 <- compute_post_fun(samps = coefsamps1, global_samps = global_samps1, 
                         knots = knots, 
                         refined_x = df_full$t,
                         p = polyOrder, degree = 0)
  
  f1_obs <- compute_post_fun(samps = coefsamps1, global_samps = global_samps1, 
                             knots = knots, 
                             refined_x = df$t,
                             p = polyOrder, degree = 0)
  
  f1deriv <- compute_post_fun(samps = coefsamps1, global_samps = global_samps1, 
                              knots = knots, 
                              refined_x = df_full$t,
                              p = polyOrder, degree = 1)
  
  f1deriv_obs <- compute_post_fun(samps = coefsamps1, global_samps = global_samps1, 
                                  knots = knots, 
                                  refined_x = df$t,
                                  p = polyOrder, degree = 1)
  
  Xfstat_full <-model.matrix(~site_id, 
                             data = df_full %>% mutate(site_id = factor(site_id)), 
                             contrasts.arg = list(site_id = "contr.sum"))[,-1]

  
  f1_IS_fixed <- f1[,-1]+ IS_samps1 + Xfstat_full%*%Xf_samps1
  f1_IS_fixed_deriv <- f1deriv[,-1]+ IS_deriv_samps1
  
  f1_IS_fixed_obs <- f1_obs[,-1]+ obs%*%IS_samps1 + Xfstat%*%Xf_samps1
  f1_IS_fixed_deriv_obs <- f1deriv_obs[,-1]+ obs%*%IS_deriv_samps1
  
  ## Nominal AR2 + ospline+ fixed effects
  df_full$exp_f1_IS_fixed <- as.numeric(apply(exp(f1_IS_fixed), MARGIN=1,median))
  df_full$exp_f1_IS_fixed_lwr <- as.numeric(apply(exp(f1_IS_fixed), MARGIN=1,quantile, 0.025))
  df_full$exp_f1_IS_fixed_upr <- as.numeric(apply(exp(f1_IS_fixed), MARGIN=1,quantile, 0.975))
  
  df_full$exp_f1_IS_fixed_deriv <- as.numeric(apply(exp(f1_IS_fixed)*(f1_IS_fixed_deriv), MARGIN=1,median))
  df_full$exp_f1_IS_fixed_deriv_lwr <- as.numeric(apply(exp(f1_IS_fixed)*(f1_IS_fixed_deriv), MARGIN=1,quantile, 0.025))
  df_full$exp_f1_IS_fixed_deriv_upr <- as.numeric(apply(exp(f1_IS_fixed)*(f1_IS_fixed_deriv), MARGIN=1,quantile, 0.975))
  
  ## Posterior prob increase 
  
  df_full$post_prob_f1_IS <- as.numeric(apply(f1_IS_fixed_deriv,MARGIN=1, FUN= function(x){length(which(x>0))/length(x)}))
  
  ## Log AR2 + ospline+ fixed effects
  df_full$f1_IS_fixed <- as.numeric(apply(f1_IS_fixed, MARGIN=1,median))
  df_full$f1_IS_fixed_lwr <- as.numeric(apply(f1_IS_fixed, MARGIN=1,quantile, 0.025))
  df_full$f1_IS_fixed_upr <- as.numeric(apply(f1_IS_fixed, MARGIN=1,quantile, 0.975))
  
  df_full$f1_IS_fixed_deriv <- as.numeric(apply(f1_IS_fixed_deriv, MARGIN=1,median))
  df_full$f1_IS_fixed_deriv_lwr <- as.numeric(apply(f1_IS_fixed_deriv, MARGIN=1,quantile, 0.025))
  df_full$f1_IS_fixed_deriv_upr <- as.numeric(apply(f1_IS_fixed_deriv, MARGIN=1,quantile, 0.975))
  
  ## Full data
  f1_IS_fixed_daily <- f1[,-1]+ IS_samps1 + Xfstat_full%*%Xf_samps1 + daily%*%daily_samps1
  
  ## Log AR2 + ospline+ fixed effects + daily
  df_full$f1_IS_fixed_daily <- as.numeric(apply(f1_IS_fixed_daily, MARGIN=1,median))
  df_full$f1_IS_fixed_daily_lwr <- as.numeric(apply(f1_IS_fixed_daily, MARGIN=1,quantile, 0.025))
  df_full$f1_IS_fixed_daily_upr <- as.numeric(apply(f1_IS_fixed_daily, MARGIN=1,quantile, 0.975))
  
  ## AR2 + ospline+ fixed effects + daily
  df_full$exp_f1_IS_fixed_daily <- as.numeric(apply(exp(f1_IS_fixed_daily), MARGIN=1,median))
  df_full$exp_f1_IS_fixed_daily_lwr <- as.numeric(apply(exp(f1_IS_fixed_daily), MARGIN=1,quantile, 0.025))
  df_full$exp_f1_IS_fixed_daily_upr <- as.numeric(apply(exp(f1_IS_fixed_daily), MARGIN=1,quantile, 0.975))
  
  ## Posterior prediction
  set.seed(2)
  cov_samps = exp(samps1$thetasamples[[3]])
  preds <- lapply(as.list(1:nrow(df_full)), function(i){
    eta = f1_IS_fixed_daily[i,]%>% t() %>% as.vector()
    pred_values <- rgamma(n=3000, shape = 1/cov_samps^2, scale = exp(eta)*cov_samps^2)
    log_pred_values = log(pred_values)
    obs_y = df_full$y[i]
    if(!is.na(obs_y)){quantile_obs = length(which(pred_values <= obs_y))/length(pred_values)}else{quantile_obs = NA}
    summary_exp = c(median(pred_values), quantile(pred_values, 0.025), quantile(pred_values, 0.975))
    summary_log = c(median(log_pred_values), quantile(log_pred_values, 0.025), quantile(log_pred_values, 0.975))
    data.frame(exp_pred = summary_exp[1],
               exp_pred_lwr = summary_exp[2],
               exp_pred_upr = summary_exp[3],
               pred = summary_log[1],
               pred_lwr = summary_log[2],
               pred_upr = summary_log[3],
               quantile_obs = quantile_obs)
  })
  
  
  preds <- Reduce(rbind, preds)
  
  df_full <- cbind(df_full, preds)
  
  post_samps_df <- df_full %>% 
    dplyr::select(1:8) %>% 
    cbind(as.data.frame(f1_IS_fixed)) %>% 
    melt(id.vars = 1:8)
  
  post_samps_df_deriv <- df_full %>% 
    dplyr::select(1:8) %>% 
    cbind(as.data.frame(f1_IS_fixed_deriv)) %>% 
    melt(id.vars = 1:8)
  
  ## Log Ospline
  g_result <- extract_mean_interval_given_samps(f1)
  
  df_full$smoothed <- g_result$mean
  df_full$smoothed_upper <- g_result$pupper
  df_full$smoothed_lower <- g_result$plower
  
  df_full$derivsmoothed <- as.numeric(apply(f1deriv[,-1], MARGIN=1,median))
  df_full$derivsmoothed_lower <- as.numeric(apply(f1deriv[,-1], MARGIN=1,quantile, 0.025))
  df_full$derivsmoothed_upper<- as.numeric(apply(f1deriv[,-1], MARGIN=1,quantile, 0.975))
  
  df_full$post_prob_f1 <- as.numeric(apply(f1deriv[,-1],MARGIN=1, FUN= function(x){length(which(x>0))/length(x)}))
  
  
  ## Ospline
  df_full$exp_smoothed <- as.numeric(apply(exp(f1[,-1]), MARGIN=1,median))
  df_full$exp_smoothed_upper <- as.numeric(apply(exp(f1[,-1]), MARGIN=1,quantile, p = 0.975))
  df_full$exp_smoothed_lower <- as.numeric(apply(exp(f1[,-1]), MARGIN=1,quantile, p = 0.025))
  
  df_full$exp_derivsmoothed <- as.numeric(apply((exp(f1[,-1]) * f1deriv[,-1]), MARGIN=1,median))
  df_full$exp_derivsmoothed_upper <- as.numeric(apply((exp(f1[,-1]) * f1deriv[,-1]), MARGIN=1,quantile, 0.975))
  df_full$exp_derivsmoothed_lower<- as.numeric(apply((exp(f1[,-1]) * f1deriv[,-1]), MARGIN=1,quantile, 0.025))
  
  post_samps_df_ospline <- df_full %>% 
    dplyr::select(1:8) %>% 
    cbind(as.data.frame(f1[,-1])) %>% 
    melt(id.vars = 1:8)
  
  post_samps_df_ospline_deriv <- df_full %>% 
    dplyr::select(1:8) %>% 
    cbind(as.data.frame(f1deriv[,-1])) %>% 
    melt(id.vars = 1:8)
  
  output_postsamps = post_samps_df %>% 
    rename("f1_IS_fixed" = value) %>% 
    cbind(post_samps_df_deriv%>% dplyr::select(value) %>% rename("f1_IS_fixed_deriv"= value) ) %>% 
    cbind(post_samps_df_ospline%>% dplyr::select(value) %>% rename("smoothed"= value) ) %>% 
    cbind(post_samps_df_ospline_deriv%>% dplyr::select(value) %>% rename("derivsmoothed"= value))
  
  df_full <- df_full %>% 
    mutate(model = paste0("model",i)) %>% 
    dplyr::select(-c("censored_y","denom","obs","nindex","t"))
  
  tmp<- output_postsamps  %>% 
    group_by(variable, sample_date) %>% 
    summarise(mean_of_logs = mean(f1_IS_fixed),
              mean_of_exps = mean(exp(f1_IS_fixed)),
              mean_of_logs_deriv = mean(f1_IS_fixed_deriv),
              mean_of_exps_deriv = mean(f1_IS_fixed_deriv*exp(f1_IS_fixed))) %>% 
    ungroup() %>% 
    group_by(sample_date) %>% 
    summarise(ave_f1_IS_fixed = median(mean_of_logs),
              ave_f1_IS_fixed_upr = quantile(mean_of_logs, 0.975),
              ave_f1_IS_fixed_lwr = quantile(mean_of_logs, 0.025),
              ave_exp_f1_IS_fixed = median(mean_of_exps),
              ave_exp_f1_IS_fixed_upr = quantile(mean_of_exps, 0.975),
              ave_exp_f1_IS_fixed_lwr = quantile(mean_of_exps, 0.025),
              ave_f1_IS_fixed_deriv = median(mean_of_logs_deriv),
              ave_f1_IS_fixed_deriv_upr = quantile(mean_of_logs_deriv, 0.975),
              ave_f1_IS_fixed_deriv_lwr = quantile(mean_of_logs_deriv, 0.025),
              ave_exp_f1_IS_fixed_deriv = median(mean_of_exps_deriv),
              ave_exp_f1_IS_fixed_deriv_upr = quantile(mean_of_exps_deriv, 0.975),
              ave_exp_f1_IS_fixed_deriv_lwr = quantile(mean_of_exps_deriv, 0.025),
              post_prob_ave_exp_f1_IS_fixed_deriv = length(which(mean_of_exps_deriv>0))/length(mean_of_exps_deriv),
              post_prob_ave_f1_IS_fixed_deriv = length(which(mean_of_logs_deriv>0))/length(mean_of_logs_deriv))
  
  
  write.csv(df_full, file = paste0('./nowcast_firstwave/results_model',i,'.csv'))
  write.csv(tmp, file = paste0('./nowcast_firstwave/CANresults_model',i,'.csv'))
}


