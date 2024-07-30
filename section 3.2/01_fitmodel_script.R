# > sessionInfo()
# R version 4.4.1 (2024-06-14) -- "Race for Your Life"
# Platform: aarch64-apple-darwin20
# Running under: macOS Monterey 12.4

# To install OSplines package
# install.packages("remotes")
# remotes::install_github("AgueroZZ/OSplines")

library(OSplines) # OSplines_0.1.1
library(aghq) # aghq_0.4.1
library(readr) # readr_2.1.5
library(dplyr) # dplyr_1.1.4
library(magrittr) # magrittr_2.0.3 
library(lubridate) # lubridate_1.9.3 
library(Matrix) # Matrix_1.7-0 
library(TMB) #TMB_1.9.14 

compute_weights_precision <- function(x){
  d <- diff(x)
  Precweights <- diag(d)
  Precweights
}

setwd('./section 3.2/')

raw_d <- read_csv("./data/SCAN_AllPlants_SDR_7April23_rev.csv")


raw_d <- raw_d %>% 
  mutate(sample_date = ymd(paste0(Year,"-",Month,"-",Day))) %>% 
  rename(site_id = `Plant Abbr`) 

# clean column names
work_d <- raw_d %>% 
  janitor::clean_names() 


## outcome name: 'y'
## If no censoring, create column censored_y = FALSE
## If no normalizing value create column: denom = 1
## name the station column: "site_id"
## name the date: "sample_date"

# f a “0” appears, it means the assay was a non-detect. 
# The detection limit varies by sample depending on the amount of solids by dry weight included, 
# but is between 500–1000 cp/g dry weight.

df <- work_d %>% 
  dplyr::mutate("y" = rsv_gc_g_dry_weight,
                "denom" = pm_mo_v_gc_g_dry_weight) %>% # no censored denom
  mutate(censored_y = ifelse(y == 0, TRUE, FALSE), 
         y = ifelse(y==0, 1000, y)) 


rr = df %>% filter(!is.na(y)) %$% sample_date %>% range()

df <- expand.grid(sample_date = seq(rr[1],rr[2],1),
                  site_id = unique(df$site_id)) %>% 
  left_join(df, by = c("sample_date","site_id")) %>% 
  dplyr::select(sample_date, site_id,  y, denom, censored_y) 

summary(df$censored_y)


# df %>% filter(censored_y == FALSE) %>% ggplot(aes(sample_date, log(y), col = site_id))+geom_line()+geom_point()+ facet_wrap(~site_id, ncol = 5)

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
knots <- seq(0,(max(df$t)), length.out = 50)

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
  y_ind_cens = y_ind_cens,
  cens_dir = rep(1, length(y_ind_cens)),
  station = as(station, 'dgTMatrix'),
  knots = knots
)
# 
compile(file="./cpp/model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper_censored.cpp")
try(dyn.unload(dynlib("./cpp/model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper_censored")),silent = TRUE)
dyn.load(dynlib("./cpp/model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper_censored"))

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
  DLL = "model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper_censored",
  silent = TRUE
)

aghq_k = 3

mdl1 <- aghq::marginal_laplace_tmb(ff,k=aghq_k,startingvalue = c(10,0,0,0,0))

samps1 <- aghq::sample_marginal(mdl1, M = 3000) # this is not working for some reson
marginals <- mdl1$marginals
save(file="~/Wastewater/toEnglish/RSV/model2.RData", list = c("df", "df_full","marginals","samps1","tmbdat","polyOrder"))
