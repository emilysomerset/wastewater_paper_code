rm(list=ls())
library(dplyr)
library(splines)
# library(OSplines)
library(TMB)
library(aghq)
library(bayesplot)
library(lemon)
library(openxlsx)
library(magrittr)

load("~/Wastewater/toEnglish/Data/work_d_xiaotian_phac.RData")

compile(file="~/Wastewater/toEnglish/cpp/model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper.cpp")
try(dyn.unload(dynlib("~/Wastewater/toEnglish/cpp/model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper")),silent = TRUE)
dyn.load(dynlib("~/Wastewater/toEnglish/cpp/model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper"))

modeldata <- modeldata2
rm(modeldata2)

raw_d <- modeldata %>% 
  rename(sample_date = date,
         site_id = Location) 

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

# clean column names
work_d <- raw_d %>% 
  janitor::clean_names() %>% 
  group_by(sample_date) %>% 
  mutate(sample_date = ymd(sample_date)) %>% 
  ungroup()




dates_with_newdata <- work_d %>% 
  filter(sample_date >= ymd('2022-12-01') & sample_date <= ymd('2023-03-31')) %>% 
  filter(!is.na(value_raw))%$% sample_date %>% unique()


# dates_totry = dates_with_newdata[seq(1, length(dates_with_newdata),7)]
dates_totry = dates_with_newdata

for (i in 102:length(dates_totry)){
  
  df <- work_d %>% 
    dplyr::mutate("y" = value_raw) 
  
  
  rr = df %>% filter(!is.na(y)) %$% sample_date %>% range()
  
  df <- expand.grid(sample_date = seq(rr[1],rr[2],1),
                    site_id = unique(df$site_id)) %>% 
    left_join(df, by = c("sample_date","site_id")) %>% 
    dplyr::select(sample_date, site_id,  y) %>% 
    mutate(censored_y = FALSE, 
           denom = 1) %>% 
    filter(sample_date <= (dates_totry[i] + 7)) %>% 
    mutate(y = ifelse(sample_date > dates_totry[i], NA_real_, y))
  
  
  df =df %>%  
    arrange(site_id,sample_date) %>% 
    ungroup() 
  
  df_full <- df %>% 
    mutate(obs = !is.na(y)) %>% 
    mutate(nindex = 1:nrow(.))
  
  df <- df %>% 
    filter(!is.na(y)) %>% 
    mutate(nindex2 = 1:nrow(.))
  
  df$t <- df$sample_date %>% as.numeric()
  df$t <- df$t - min(df$t)
  df_full$t <- df_full$sample_date %>% as.numeric()
  df_full$t <- df_full$t - min(df_full$t)
  
  stationsizes = df_full %>%  
    group_by(site_id) %>%  
    summarise(ll = length(site_id)) %$% ll
  
  n1 = df_full %>%  
    mutate(n1 = 0:(nrow(.)-1)) %>%  
    group_by(site_id) %>%  
    slice(1) %$% n1 
  
  n1 <- c(n1, (nrow(df_full)))
  
  # Setup the design and precision matrix 
  polyOrder = 3
  knots <- seq(0,(max(df$t)), length.out = ceiling((max(df$t))/7))
  
  B <- local_poly_helper(knots = knots, 
                         refined_x = df$t,
                         p = polyOrder)
  
  P <- compute_weights_precision(x=knots)
  
  X_global = global_poly_helper(x=df$t, p= polyOrder)
  Xfstat <- model.matrix(~site_id, 
                         data = df %>% mutate(site_id = factor(site_id)), 
                         contrasts.arg = list(site_id = "contr.sum"))[,-1]
  Xfstat <- as.matrix(Xfstat)
  X = cbind(X_global, Xfstat)
  
  # daily <- model.matrix(~t -1, 
  #                       data = df %>% mutate(t = factor(t)))
  
  daily <- model.matrix(~t -1,
                        data = df_full %>% mutate(t = factor(t)))
  
  
  obs <- model.matrix(~nindex -1, data = df_full %>% mutate(nindex=factor(nindex)) %>% filter(obs))
  station <- model.matrix(~factor(site_id) -1, data = df)
  y_ind_obs <- df %>% filter(censored_y==FALSE)%$%nindex2 -1
  y_ind_cens <- df %>% filter(censored_y==TRUE)%$%nindex2 - 1
  observed_y <- model.matrix(~nindex2 -1, data = df %>% mutate(nindex2=factor(nindex2)) %>% filter(censored_y==FALSE))
  censored_y <- model.matrix(~nindex2 -1, data = df %>% mutate(nindex2=factor(nindex2)) %>% filter(censored_y==TRUE))
  
  tmbdat <- list(
    # Design matrix
    Xfstat = as(Xfstat, 'dgTMatrix'),
    X_global = as(X_global, 'dgTMatrix'),
    X = as(X, 'dgTMatrix'),
    B = as(B, 'dgTMatrix'),
    P = as(P,'dgTMatrix'),
    daily = as(daily, 'dgTMatrix'),
    obs = as(obs, 'dgTMatrix'),
    logPdet = as.numeric(determinant(P,logarithm = T)$modulus),
    # Response
    y = df$y,
    # PC Prior params
    n1=n1,
    stationsizes = stationsizes,
    denom = df$denom,
    y_ind_obs=y_ind_obs,
    # y_ind_cens = y_ind_cens,
    # cens_dir = rep(1, length(y_ind_cens)),
    station = as(station, 'dgTMatrix'),
    knots = knots
  )
  
  
  prior_IWP <- prior_conversion(d=20, prior =list(u=log(2),alpha = 0.5),p=polyOrder)
  
  # Set other priors
  tmbdat$u1 = prior_IWP$u
  tmbdat$alpha1 = 0.5
  tmbdat$u2 = 0.5
  tmbdat$alpha2 = 0.5
  tmbdat$betaprec = 0.01
  tmbdat$lambda_phi = -log(0.5)/50
  tmbdat$lambda_tau = -log(0.5)/0.5
  tmbdat$lambda_cov = -log(0.5)/0.5
  
  set.seed(2)
  init_daily <- rnorm(ncol(tmbdat$daily), 0, 0.1);
  init_W <- rnorm(ncol(tmbdat$obs),0,0.1)
  tmbparams <- list(
    W = c(rep(0, (ncol(tmbdat$X)+ ncol(tmbdat$B))), init_daily, init_W, 0,diff(init_W)), # W = c(U,beta,Z); U = B-Spline coefficients
    theta1 = 10, # -2log(sigma)
    theta2 = 0,
    cov_log = 0,
    theta3 = 0,
    theta4 = 0
  )
  
  
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    random = "W",
    DLL = "model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper",
    silent = TRUE
  )
  
  aghq_k = 3
  
  mdl1 <- aghq::marginal_laplace_tmb(ff,k=aghq_k,startingvalue = c(10,0,0,0,0))
  
  samps1 <- aghq::sample_marginal(mdl1, M = 3000) 
  marginals <- mdl1$marginals
  save(file=paste0("~/Wastewater/toEnglish/xiaotian_phac/nowcast_nowave/model",i,".RData"), list = c("df", "df_full","marginals","samps1","tmbdat","polyOrder"))
  rm(list = c("mdl1","samps1","marginals"))
  print(i)
}
