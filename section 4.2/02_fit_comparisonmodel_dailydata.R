# > sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: aarch64-apple-darwin20
# Running under: macOS Monterey 12.4

library(INLA)        #INLA_24.06.27 
library(dplyr)       #dplyr_1.1.4 
library(reshape2)    #reshape2_1.4.4
library(data.table)  #data.table_1.15.4


load("./data/work_d.RData")

ww.daily <- ww.daily %>%
  filter(RegionName == "London")
n.stat <- length(unique(ww.daily$uwwName))
n = length(unique(ww.daily$time_index))

data_tmp <- list(y = ww.daily$log_e_gc,
                 time_index = rep(1:n,n.stat), 
                 idx.rep = rep(1:n.stat, each = n),
                 idx.ts = rep(1:n, n.stat), 
                 day = rep(1:n, n.stat))


pcprior.rho = list(theta = list(prior = "pc.cor0", param = c(0.1,0.9)))
hyper.prec = list(prec = list(prior = "pc.prec", param = c(10,0.05)))
priors_fixed <- list(fixed = list(mean.intercept = 0, prec.intercept = 1e-04,
                                mean = 0, prec = 1e-04))

formula <- as.formula('y~ f(idx.ts, model = "ar1", hyper = pcprior.rho, replicate = idx.rep)+
                      f(time_index, model = "rw1", hyper = hyper.prec)+
                      f(day, model = "iid", hyper = hyper.prec) ')


outputs <- inla(formula, data = data_tmp, family = "gaussian", 
                control.fixed = priors_fixed$fixed, 
                control.compute = list(config = TRUE, return.marginals.predictor = TRUE),
                verbose = FALSE)

nsims =1000
samples <- INLA::inla.posterior.sample(n=nsims, outputs)


gsub(':[[:digit:]]+', '', rownames(samples[[1]]$latent)) %>% unique()
# "Predictor"  "idx.ts"     "time_index" "day"       "(Intercept)" 

a <- unlist(lapply(samples, function(s) s$latent)) # this is the components of the link function
a <- matrix(a, nrow = 1000, byrow = TRUE)
dim(a) #1000 X 5455


# just get the predictor for now
rr <- range(which(gsub(':[[:digit:]]+', '', rownames(samples[[1]]$latent))== "Predictor"))
a2 <- a[,1:rr[2]]
a2 <- t(a2)

pred <- apply(a2, 1, median)
pred_lwr <- apply(a2, 1, quantile, 0.025)
pred_upr <- apply(a2, 1, quantile, 0.975)


## full signal 
rr1 <- range(which(gsub(':[[:digit:]]+', '', rownames(samples[[1]]$latent))== "time_index"))
rr2 <- range(which(gsub(':[[:digit:]]+', '', rownames(samples[[1]]$latent))== "idx.ts"))
rr3 <- range(which(gsub(':[[:digit:]]+', '', rownames(samples[[1]]$latent))== "(Intercept)"))


dd <- model.matrix(~factor(rep(1:n.stat, each = n)))
dd2 <- model.matrix(~factor(rep(1:n, n.stat))-1)


signal <- (matrix(1, ncol = 1, nrow = nrow(dd2))%*% t(a[,rr3[1]])) +  t(a[,rr2[1]:rr2[2]]) + (dd2 %*% t(a[,rr1[1]:rr1[2]])) ## full_signal

full_signal_results <- ww.daily %>% 
  ungroup() %>% 
  cbind(signal) %>% 
  melt(measure.vars = paste0(1:1000)) %>% 
  group_by(variable, uwwName) %>% 
  mutate(value_diff = c(NA, diff(value)),
         value_diff_7 = c(rep(NA,7), diff(value,lag=7)),
         exp_value_diff = c(NA, diff(exp(value))),
         exp_value_diff_7 = c(rep(NA,7), diff(exp(value),lag=7))) %>% 
  # dplyr::select(uwwName, Date,time_index,variable, value_diff_7,value) %>% 
  mutate(which_pos = value_diff_7>0 & shift(value_diff_7,7)>0,
         exp_which_pos = exp_value_diff_7>0 & shift(exp_value_diff_7,7)>0) %>% 
  mutate(which_pos = ifelse(time_index <= 14,NA, which_pos),
         exp_which_pos = ifelse(time_index <= 14,NA, exp_which_pos)) %>% 
  group_by(Date, time_index,site_index, uwwName, RegionName, gc,log_e_gc) %>% 
  summarise(f1_IS_fixed = median(value),
            f1_IS_fixed_upr = quantile(value, 0.975),
            f1_IS_fixed_lwr = quantile(value, 0.025),
            exp_f1_IS_fixed = median(exp(value)),
            exp_f1_IS_fixed_upr = quantile(exp(value), 0.975),
            exp_f1_IS_fixed_lwr = quantile(exp(value), 0.025),
            f1_IS_fixed_deriv = median(value_diff, na.rm = TRUE),
            f1_IS_fixed_deriv_upr = quantile(value_diff, 0.975, na.rm = TRUE),
            f1_IS_fixed_deriv_lwr = quantile(value_diff, 0.025, na.rm = TRUE),
            exp_f1_IS_fixed_deriv = median(exp_value_diff, na.rm = TRUE),
            exp_f1_IS_fixed_deriv_upr = quantile(exp_value_diff, 0.975, na.rm = TRUE),
            exp_f1_IS_fixed_deriv_lwr = quantile(exp_value_diff, 0.025, na.rm = TRUE),
            f1_IS_fixed_deriv_7 = median(value_diff_7, na.rm = TRUE),
            f1_IS_fixed_deriv_7_upr = quantile(value_diff_7, 0.975, na.rm = TRUE),
            f1_IS_fixed_deriv_7_lwr = quantile(value_diff_7, 0.025, na.rm = TRUE),
            exp_f1_IS_fixed_deriv_7 = median(exp_value_diff_7, na.rm = TRUE),
            exp_f1_IS_fixed_deriv_7_upr = quantile(exp_value_diff_7, 0.975, na.rm = TRUE),
            exp_f1_IS_fixed_deriv_7_lwr = quantile(exp_value_diff_7, 0.025, na.rm = TRUE),
            f1_IS_fixed_which_pos = length(which(which_pos))/length(which_pos),
            exp_f1_IS_fixed_which_pos = length(which(exp_which_pos))/length(exp_which_pos))


signal <- (matrix(1, ncol = 1, nrow = nrow(dd2))%*% t(a[,rr3[1]])) + t(a[,rr2[1]:rr2[2]]) 

IS_fixed_signal_results <- ww.daily %>% 
  ungroup() %>% 
  cbind(signal) %>% 
  melt(measure.vars = paste0(1:1000)) %>% 
  group_by(variable, uwwName) %>% 
  mutate(value_diff = c(NA, diff(value)),
         value_diff_7 = c(rep(NA,7), diff(value,lag=7)),
         exp_value_diff = c(NA, diff(exp(value))),
         exp_value_diff_7 = c(rep(NA,7), diff(exp(value),lag=7))) %>% 
  # dplyr::select(uwwName, Date,time_index,variable, value_diff_7,value) %>% 
  mutate(which_pos = value_diff_7>0 & shift(value_diff_7,7)>0,
         exp_which_pos = exp_value_diff_7>0 & shift(exp_value_diff_7,7)>0) %>% 
  mutate(which_pos = ifelse(time_index <= 14,NA, which_pos),
         exp_which_pos = ifelse(time_index <= 14,NA, exp_which_pos)) %>% 
  group_by(Date, time_index,site_index, uwwName, RegionName, gc,log_e_gc) %>% 
  summarise(IS_fixed = median(value),
            IS_fixed_upr = quantile(value, 0.975),
            IS_fixed_lwr = quantile(value, 0.025),
            exp_IS_fixed = median(exp(value)),
            exp_IS_fixed_upr = quantile(exp(value), 0.975),
            exp_IS_fixed_lwr = quantile(exp(value), 0.025),
            IS_fixed_deriv = median(value_diff, na.rm = TRUE),
            IS_fixed_deriv_upr = quantile(value_diff, 0.975, na.rm = TRUE),
            IS_fixed_deriv_lwr = quantile(value_diff, 0.025, na.rm = TRUE),
            exp_IS_fixed_deriv = median(exp_value_diff, na.rm = TRUE),
            exp_IS_fixed_deriv_upr = quantile(exp_value_diff, 0.975, na.rm = TRUE),
            exp_IS_fixed_deriv_lwr = quantile(exp_value_diff, 0.025, na.rm = TRUE),
            IS_fixed_deriv_7 = median(value_diff_7, na.rm = TRUE),
            IS_fixed_deriv_7_upr = quantile(value_diff_7, 0.975, na.rm = TRUE),
            IS_fixed_deriv_7_lwr = quantile(value_diff_7, 0.025, na.rm = TRUE),
            exp_IS_fixed_deriv_7 = median(exp_value_diff_7, na.rm = TRUE),
            exp_IS_fixed_deriv_7_upr = quantile(exp_value_diff_7, 0.975, na.rm = TRUE),
            exp_IS_fixed_deriv_7_lwr = quantile(exp_value_diff_7, 0.025, na.rm = TRUE),
            IS_fixed_which_pos = length(which(which_pos))/length(which_pos),
            exp_IS_fixed_which_pos = length(which(exp_which_pos))/length(exp_which_pos))

signal <- t(a[,rr2[1]:rr2[2]]) 

IS_signal_results <- ww.daily %>% 
  ungroup() %>% 
  cbind(signal) %>% 
  melt(measure.vars = paste0(1:1000)) %>% 
  group_by(variable, uwwName) %>% 
  mutate(value_diff = c(NA, diff(value)),
         value_diff_7 = c(rep(NA,7), diff(value,lag=7)),
         exp_value_diff = c(NA, diff(exp(value))),
         exp_value_diff_7 = c(rep(NA,7), diff(exp(value),lag=7))) %>% 
  mutate(which_pos = value_diff_7>0 & shift(value_diff_7,7)>0,
         exp_which_pos = exp_value_diff_7>0 & shift(exp_value_diff_7,7)>0) %>% 
  mutate(which_pos = ifelse(time_index <= 14,NA, which_pos),
         exp_which_pos = ifelse(time_index <= 14,NA, exp_which_pos)) %>% 
  group_by(Date, time_index,site_index, uwwName, RegionName, gc,log_e_gc) %>% 
  summarise(IS = median(value),
            IS_upr = quantile(value, 0.975),
            IS_lwr = quantile(value, 0.025),
            exp_IS = median(exp(value)),
            exp_IS_upr = quantile(exp(value), 0.975),
            exp_IS_lwr = quantile(exp(value), 0.025),
            IS_deriv = median(value_diff, na.rm = TRUE),
            IS_deriv_upr = quantile(value_diff, 0.975, na.rm = TRUE),
            IS_deriv_lwr = quantile(value_diff, 0.025, na.rm = TRUE),
            exp_IS_deriv = median(exp_value_diff, na.rm = TRUE),
            exp_IS_deriv_upr = quantile(exp_value_diff, 0.975, na.rm = TRUE),
            exp_IS_deriv_lwr = quantile(exp_value_diff, 0.025, na.rm = TRUE),
            IS_deriv_7 = median(value_diff_7, na.rm = TRUE),
            IS_deriv_7_upr = quantile(value_diff_7, 0.975, na.rm = TRUE),
            IS_deriv_7_lwr = quantile(value_diff_7, 0.025, na.rm = TRUE),
            exp_IS_deriv_7 = median(exp_value_diff_7, na.rm = TRUE),
            exp_IS_deriv_7_upr = quantile(exp_value_diff_7, 0.975, na.rm = TRUE),
            exp_IS_deriv_7_lwr = quantile(exp_value_diff_7, 0.025, na.rm = TRUE),
            IS_which_pos = length(which(which_pos))/length(which_pos),
            exp_IS_which_pos = length(which(exp_which_pos))/length(exp_which_pos))




signal <- dd2 %*% t(a[,rr1[1]:rr1[2]]) 

f1_signal_results <- ww.daily %>% 
  ungroup() %>% 
  cbind(signal) %>% 
  melt(measure.vars = paste0(1:1000)) %>% 
  group_by(variable, uwwName) %>% 
  mutate(value_diff = c(NA, diff(value)),
         value_diff_7 = c(rep(NA,7), diff(value,lag=7)),
         exp_value_diff = c(NA, diff(exp(value))),
         exp_value_diff_7 = c(rep(NA,7), diff(exp(value),lag=7))) %>% 
  mutate(which_pos = value_diff_7>0 & shift(value_diff_7,7)>0,
         exp_which_pos = exp_value_diff_7>0 & shift(exp_value_diff_7,7)>0) %>% 
  mutate(which_pos = ifelse(time_index <= 14,NA, which_pos),
         exp_which_pos = ifelse(time_index <= 14,NA, exp_which_pos)) %>% 
  group_by(Date, time_index,site_index, uwwName, RegionName, gc,log_e_gc) %>% 
  summarise(f1 = median(value),
            f1_upr = quantile(value, 0.975),
            f1_lwr = quantile(value, 0.025),
            exp_f1 = median(exp(value)),
            exp_f1_upr = quantile(exp(value), 0.975),
            exp_f1_lwr = quantile(exp(value), 0.025),
            f1_deriv = median(value_diff, na.rm = TRUE),
            f1_deriv_upr = quantile(value_diff, 0.975, na.rm = TRUE),
            f1_deriv_lwr = quantile(value_diff, 0.025, na.rm = TRUE),
            exp_f1_deriv = median(exp_value_diff, na.rm = TRUE),
            exp_f1_deriv_upr = quantile(exp_value_diff, 0.975, na.rm = TRUE),
            exp_f1_deriv_lwr = quantile(exp_value_diff, 0.025, na.rm = TRUE),
            f1_deriv_7 = median(value_diff_7, na.rm = TRUE),
            f1_deriv_7_upr = quantile(value_diff_7, 0.975, na.rm = TRUE),
            f1_deriv_7_lwr = quantile(value_diff_7, 0.025, na.rm = TRUE),
            exp_f1_deriv_7 = median(exp_value_diff_7, na.rm = TRUE),
            exp_f1_deriv_7_upr = quantile(exp_value_diff_7, 0.975, na.rm = TRUE),
            exp_f1_deriv_7_lwr = quantile(exp_value_diff_7, 0.025, na.rm = TRUE),
            f1_which_pos = length(which(which_pos))/length(which_pos),
            exp_f1_which_pos = length(which(exp_which_pos))/length(exp_which_pos))

signal <- (matrix(1, ncol = 1, nrow = nrow(dd2))%*% t(a[,rr3[1]])) + (dd2 %*% t(a[,rr1[1]:rr1[2]])) ## common signal 2414 X 1000

f1_int_signal_results <- ww.daily %>% 
  ungroup() %>% 
  cbind(signal) %>% 
  melt(measure.vars = paste0(1:1000)) %>% 
  group_by(variable, uwwName) %>% 
  mutate(value_diff = c(NA, diff(value)),
         value_diff_7 = c(rep(NA,7), diff(value,lag=7)),
         exp_value_diff = c(NA, diff(exp(value))),
         exp_value_diff_7 = c(rep(NA,7), diff(exp(value),lag=7))) %>% 
  mutate(which_pos = value_diff_7>0 & shift(value_diff_7,7)>0,
         exp_which_pos = exp_value_diff_7>0 & shift(exp_value_diff_7,7)>0) %>% 
  mutate(which_pos = ifelse(time_index <= 14,NA, which_pos),
         exp_which_pos = ifelse(time_index <= 14,NA, exp_which_pos)) %>% 
  group_by(Date, time_index,site_index, uwwName, RegionName, gc,log_e_gc) %>% 
  summarise(f1_int = median(value),
            f1_int_upr = quantile(value, 0.975),
            f1_int_lwr = quantile(value, 0.025),
            exp_f1_int = median(exp(value)),
            exp_f1_int_upr = quantile(exp(value), 0.975),
            exp_f1_int_lwr = quantile(exp(value), 0.025),
            f1_int_deriv = median(value_diff, na.rm = TRUE),
            f1_int_deriv_upr = quantile(value_diff, 0.975, na.rm = TRUE),
            f1_int_deriv_lwr = quantile(value_diff, 0.025, na.rm = TRUE),
            exp_f1_int_deriv = median(exp_value_diff, na.rm = TRUE),
            exp_f1_int_deriv_upr = quantile(exp_value_diff, 0.975, na.rm = TRUE),
            exp_f1_int_deriv_lwr = quantile(exp_value_diff, 0.025, na.rm = TRUE),
            f1_int_deriv_7 = median(value_diff_7, na.rm = TRUE),
            f1_int_deriv_7_upr = quantile(value_diff_7, 0.975, na.rm = TRUE),
            f1_int_deriv_7_lwr = quantile(value_diff_7, 0.025, na.rm = TRUE),
            exp_f1_int_deriv_7 = median(exp_value_diff_7, na.rm = TRUE),
            exp_f1_int_deriv_7_upr = quantile(exp_value_diff_7, 0.975, na.rm = TRUE),
            exp_f1_int_deriv_7_lwr = quantile(exp_value_diff_7, 0.025, na.rm = TRUE),
            f1_int_which_pos = length(which(which_pos))/length(which_pos),
            exp_f1_int_which_pos = length(which(exp_which_pos))/length(exp_which_pos))


full_results <- full_signal_results %>% 
  left_join(IS_fixed_signal_results, by = c('Date', 'time_index','site_index', 'uwwName', 'RegionName', 'gc','log_e_gc')) %>% 
  left_join(IS_signal_results, by = c('Date', 'time_index','site_index', 'uwwName', 'RegionName', 'gc','log_e_gc')) %>% 
  left_join(f1_signal_results, by = c('Date', 'time_index','site_index', 'uwwName', 'RegionName', 'gc','log_e_gc')) %>% 
  left_join(f1_int_signal_results, by = c('Date', 'time_index','site_index', 'uwwName', 'RegionName', 'gc','log_e_gc'))


save(file = "./results_comparisonmodel.RData", full_results)
