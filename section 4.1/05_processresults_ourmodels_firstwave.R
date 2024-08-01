rm(list=ls())

# library(fanetc)
library(splines)
library(cowplot)
# library(OSplines)
library(TMB)
library(aghq)
library(bayesplot)
library(lemon)
library(dplyr)
library(magrittr)
library(lubridate)
library(reshape2)

compute_weights_precision <- function(x){
  d <- diff(x)
  Precweights <- diag(d)
  Precweights
}
local_poly_helper <- function (knots, refined_x, p = 2) {
  if (min(knots) >= 0) {
    dif <- diff(knots)
    nn <- length(refined_x)
    n <- length(knots)
    D <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x[j] <= knots[i]) {
          D[j, i] <- 0
        }
        else if (refined_x[j] <= knots[i + 1] & refined_x[j] >= 
                 knots[i]) {
          D[j, i] <- (1/factorial(p)) * (refined_x[j] - 
                                           knots[i])^p
        }
        else {
          k <- 1:p
          D[j, i] <- sum((dif[i]^k) * ((refined_x[j] - 
                                          knots[i + 1])^(p - k))/(factorial(k) * factorial(p - 
                                                                                             k)))
        }
      }
    }
  }
  else if (max(knots) <= 0) {
    refined_x_neg <- refined_x
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    dif <- diff(knots_neg)
    nn <- length(refined_x_neg)
    n <- length(knots_neg)
    D <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_neg[j] <= knots_neg[i]) {
          D[j, i] <- 0
        }
        else if (refined_x_neg[j] <= knots_neg[i + 1] & 
                 refined_x_neg[j] >= knots_neg[i]) {
          D[j, i] <- (1/factorial(p)) * (refined_x_neg[j] - 
                                           knots_neg[i])^p
        }
        else {
          k <- 1:p
          D[j, i] <- sum((dif[i]^k) * ((refined_x_neg[j] - 
                                          knots_neg[i + 1])^(p - k))/(factorial(k) * 
                                                                        factorial(p - k)))
        }
      }
    }
  }
  else {
    refined_x_neg <- refined_x
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    dif <- diff(knots_neg)
    nn <- length(refined_x_neg)
    n <- length(knots_neg)
    D1 <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_neg[j] <= knots_neg[i]) {
          D1[j, i] <- 0
        }
        else if (refined_x_neg[j] <= knots_neg[i + 1] & 
                 refined_x_neg[j] >= knots_neg[i]) {
          D1[j, i] <- (1/factorial(p)) * (refined_x_neg[j] - 
                                            knots_neg[i])^p
        }
        else {
          k <- 1:p
          D1[j, i] <- sum((dif[i]^k) * ((refined_x_neg[j] - 
                                           knots_neg[i + 1])^(p - k))/(factorial(k) * 
                                                                         factorial(p - k)))
        }
      }
    }
    refined_x_pos <- refined_x
    refined_x_pos <- ifelse(refined_x > 0, refined_x, 0)
    knots_pos <- knots
    knots_pos <- unique(sort(ifelse(knots > 0, knots, 0)))
    dif <- diff(knots_pos)
    nn <- length(refined_x_pos)
    n <- length(knots_pos)
    D2 <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_pos[j] <= knots_pos[i]) {
          D2[j, i] <- 0
        }
        else if (refined_x_pos[j] <= knots_pos[i + 1] & 
                 refined_x_pos[j] >= knots_pos[i]) {
          D2[j, i] <- (1/factorial(p)) * (refined_x_pos[j] - 
                                            knots_pos[i])^p
        }
        else {
          k <- 1:p
          D2[j, i] <- sum((dif[i]^k) * ((refined_x_pos[j] - 
                                           knots_pos[i + 1])^(p - k))/(factorial(k) * 
                                                                         factorial(p - k)))
        }
      }
    }
    D <- cbind(D1, D2)
  }
  D
}
global_poly_helper <- function (x, p = 2) {
  result <- NULL
  for (i in 1:p) {
    result <- cbind(result, x^(i - 1))
  }
  result
}
prior_conversion <- function (d, prior, p) {
  Cp <- (d^((2 * p) - 1))/(((2 * p) - 1) * (factorial(p - 1)^2))
  prior_q <- list(a = prior$a, u = (prior$u * (1/sqrt(Cp))))
  prior_q
}
compute_post_fun <- function (samps, global_samps = NULL, knots, refined_x, p, degree = 0) {
  if (p <= degree) {
    return(message("Error: The degree of derivative to compute is not defined. Should consider higher order smoothing model or lower order of the derivative degree."))
  }
  if (is.null(global_samps)) {
    global_samps = matrix(0, nrow = p, ncol = ncol(samps))
  }
  if (nrow(global_samps) != p | nrow(samps) != (length(knots) - 
                                                1)) {
    return(message("Error: Incorrect dimension of global_samps or samps. Check whether the choice of p or the choice of knots are consistent with the fitted model."))
  }
  if (ncol(samps) != ncol(global_samps)) {
    return(message("Error: The numbers of posterior samples do not match between the O-splines and global polynomials."))
  }
  X = global_poly(refined_x, p = p)
  X <- as.matrix(X[, 1:(p - degree)])
  for (i in 1:ncol(X)) {
    X[, i] <- (factorial(i + degree - 1)/factorial(i - 1)) * 
      X[, i]
  }
  B = as(local_poly(knots, refined_x = refined_x, p = (p - 
                                                         degree)), "dgTMatrix")
  fitted_samps_deriv <- X %*% global_samps[(1 + degree):p, 
  ] + B %*% samps
  result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_deriv)))
  result
}
global_poly <- function (x, p = 2) {
  result <- NULL
  for (i in 1:p) {
    result <- cbind(result, x^(i - 1))
  }
  result
}
local_poly <- function (knots, refined_x, p = 2) {
  dif <- diff(knots)
  nn <- length(refined_x)
  n <- length(knots)
  D <- matrix(0, nrow = nn, ncol = n - 1)
  for (j in 1:nn) {
    for (i in 1:(n - 1)) {
      if (refined_x[j] <= knots[i]) {
        D[j, i] <- 0
      }
      else if (refined_x[j] <= knots[i + 1] & refined_x[j] >= 
               knots[i]) {
        D[j, i] <- (1/factorial(p)) * (refined_x[j] - 
                                         knots[i])^p
      }
      else {
        k <- 1:p
        D[j, i] <- sum((dif[i]^k) * ((refined_x[j] - 
                                        knots[i + 1])^(p - k))/(factorial(k) * factorial(p - 
                                                                                           k)))
      }
    }
  }
  D
}
extract_mean_interval_given_samps <- function (samps, level = 0.95) {
  x <- samps[, 1]
  samples <- samps[, -1]
  result <- data.frame(x = x)
  alpha <- 1 - level
  result$plower <- as.numeric(apply(samples, MARGIN = 1, quantile, 
                                    p = (alpha/2)))
  result$pupper <- as.numeric(apply(samples, MARGIN = 1, quantile, 
                                    p = (level + (alpha/2))))
  result$mean <- as.numeric(apply(samples, MARGIN = 1, mean))
  result
}
## I did something stupid and didn't save the true values of y for the prediction stuff. 
## Need to fill them in so I can get quantile estimates. 
fullmodel <- read.csv('~/Wastewater/toEnglish/xiaotian_phac/nowcast_firstwave/results_fullmodel.csv')
fullmodel <- fullmodel %>% 
  mutate(sample_date = ymd(sample_date)) %>% 
  dplyr::select(- 'X')

for (i in 10:11){
  load(paste0("~/Wastewater/toEnglish/xiaotian_phac/nowcast_firstwave/model",i,'.RData'))
  
  df_full <- df_full %>% 
    dplyr::select(-'y') %>% 
    left_join(fullmodel %>% dplyr::select(y, sample_date, site_id), by = c("sample_date","site_id")) # add the data (pooor planning)
  
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
  
  # norm <- Xfstat_full
  # norm[,2] <- -1
  # norm[,1]<- -1
  
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
  
  
  # ggplot(df_full %>% filter(sample_date <= "2021-05-23") %>% filter(!is.na(y)), aes(sample_date, exp_pred))+
  #   geom_ribbon(aes(ymax = exp_pred_upr, ymin = exp_pred_lwr),alpha = 0.2)+
  #   geom_line()+
  #   geom_point(aes(sample_date, y), size = 0.5)+
  #   facet_wrap(~site_id)+
  #   geom_vline(xintercept = ymd("2021-12-04"))
  # 
  # df_full %>% filter(sample_date <= "2021-12-03") %>% 
  #   mutate(in_range = y <= exp_pred_upr & y >= exp_pred_lwr) %>%
  #   filter(!is.na(y)) %>% 
  #   summarise(out = length(which(in_range == TRUE))/length(in_range))  
  
  
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
  
  # if (i %in% c(17,20,21)){
  #   save(file = paste0('./Xiaotian et al/phac/nowcast_firstwave2/postsamps_model',i,'.RData'), list = 'output_postsamps')
  # }
  
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
  
  
  write.csv(df_full, file = paste0('~/Wastewater/toEnglish/xiaotian_phac/nowcast_firstwave/results_model',i,'.csv'))
  write.csv(tmp, file = paste0('~/Wastewater/toEnglish/xiaotian_phac/nowcast_firstwave/CANresults_model',i,'.csv'))
}


