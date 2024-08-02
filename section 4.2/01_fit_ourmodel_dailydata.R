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

# load functions
source('prep_data.R')
source('process_results.R')

#compile C++ code
compile(file="./cpp/model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper_censored.cpp")
try(dyn.unload(dynlib("./cpp/model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper_censored")),silent = TRUE)
dyn.load(dynlib("./cpp/model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper_censored"))

# Filter London stations
ww.daily <- ww.daily %>% 
  filter(RegionName == "London") 

L = 200
ww.daily <- ww.daily %>% 
  mutate(denom = 1, 
         censored_y = ifelse(gc == L/2, TRUE, FALSE),
         gc = ifelse(gc == L/2, L, gc)) 


data_foranalysis <- prep_data(outcome_column_name = "gc",
                              denom_column_name = "denom",
                              site_id = "uwwName",   # column for wastewater station IDs
                              sample_date = "Date",  #date column. Should be year-month-day
                              data = ww.daily, # put the working data set
                              polyOrder =3,  #polynomial order for the Ospline (we recommend 3 for wastewater data)
                              pred_also = FALSE,
                              nknots = 40) 

tmbdat <- data_foranalysis$tmbdat
df_full <- data_foranalysis$df_full

polyOrder = 3
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
  theta1 = 0, # -2log(sigma)
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

samps1 <- aghq::sample_marginal(mdl1, M = 3000) #
marginals <- mdl1$marginals
save(file="model1.RData", list = c("df", "df_full","marginals","samps1","tmbdat","polyOrder"))


### Process and save the results. 

results <- process_results(df_full, tmbdat, samps1, polyOrder, id_group = NULL)
save(file = "results_ourmodel.RData", list = "results")
