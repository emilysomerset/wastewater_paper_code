rm(list=ls())
library(dplyr)
library(splines)
library(TMB)
library(aghq)
library(bayesplot)
library(lemon)
library(magrittr)
# library(fanetc)
library(data.table)
library(refund)
library(readr)
library(rstan)
library(lubridate)
source("~/Wastewater/toEnglish/WWmodel_meanconcentrate.R")
source("~/Wastewater/toEnglish/WW_forecast.R")

##########################################
#### Data preprocess
##########################################
raw_d <- readr::read_rds(file = '~/Wastewater/toEnglish/Data/ww-db-2021-11-18.rds')


RawData = as.data.table(raw_d)
RawData<- RawData %>% 
  dplyr::select("Location",  "date", "method", "target", "replicate", "labname", "value.raw","uqnd")
lab = "all"
city = "all"
model_target = "N1"

DataPrep = function(RawData, model_target, lab = "all", city = "all",covariate = NULL) {
  data = RawData %>%  
    filter(!(Location %in% c(NA, "OTW","TLV"))) 
  
  data$date = as.Date(data$date)
  data = data[labname == "NML", ]
  data = data[replicate <= 2, ]
  
  normdata = data[target != "PMMV" & !uqnd %in% c("UQ", "ND"), ]
  # normdata = data[target != "PMMV", ]  # Emily edit removed uqnd
  
  normdata = normdata[, c("Location", "date", "method", "target", "replicate", "labname", "value.raw")]
  
  widedata = reshape(normdata,
                     idvar = c("Location", "method", "target", "replicate", "labname"),
                     timevar = "date", direction = "wide")
  
  widedata = widedata[order(Location, target, replicate, method)]
  
  wdata = as.data.frame(matrix(NA, nrow = dim(widedata)[1] / 2, ncol = dim(widedata)[2] - 2))
  wdata[, 1:3] = unique(widedata[, c("Location", "target", "replicate")])
  
  names(wdata) = c("Location", "target", "replicate", names(widedata)[6:dim(widedata)[2]])
  
  for (i in 1:dim(wdata)[1]) {
    wdata[i, 4:dim(wdata)[2]] =
      as.numeric(apply(widedata[(2 * i - 1):(2 * i), 6:dim(widedata)[2]],
                       2, function(x) ifelse(is.na(x[2]), x[1], x[2])))
  }
  modeldata = as.data.table(melt(wdata, id.vars = c("Location", "target", "replicate"),
                                 variable.name = "date", value.name = "value.raw", na.rm = TRUE))
  modeldata$date = as.Date(gsub("value.raw.", "", modeldata$date))
  
  if (!is.null(covariate)){
    modeldata = merge(modeldata,
                      data[method == "solids",
                           c("Location", "date", "target", "replicate",
                             "dinflvol", "infltemp", "ph", "temp",
                             "tss", "rain_t", "rain_y")],
                      by = c("Location", "target", "replicate", "date"),
                      all =  FALSE, all.x = TRUE, all.y = FALSE)} else {
                        modeldata = merge(modeldata,
                                          data[method == "solids",
                                               c("Location", "date", "target", "replicate")],
                                          by = c("Location", "target", "replicate", "date"),
                                          all =  FALSE, all.x = TRUE, all.y = FALSE)                   
                      }
  modeldata = unique(modeldata)
  modeldata$log10.value.raw = log10(modeldata$value.raw)
  
  modeldata = modeldata[order(date)]
  modeldata = modeldata[target %in% model_target, ]
  if (lab != "all") {
    modeldata = modeldata[Location %in% lab, ]
  }
  if (city != "all") {
    sites = c()
    if ("Edmonton" %in% city) sites = c(sites, "EGB")
    if ("Halifax" %in% city) sites = c(sites, "HDA", "HHA", "HMC")
    if ("Montreal" %in% city) sites = c(sites, "MMN", "MMS")
    if ("Toronto" %in% city) sites = c(sites, "TAB", "THC", "THU", "TNT")
    if ("Vancouver" %in% city) sites = c(sites, "VAI", "VII", "VLG", "VLI", "VNL")
    modeldata = modeldata[Location %in% sites, ]
  }
  return(modeldata)
}
modeldata <- DataPrep(RawData, "N2")
data_orig <- range(modeldata$date)

# Import up-to-date PHAC 
phac_d <- read_csv("~/Wastewater/toEnglish/Data/covid19-wastewater.csv")

df = data.frame(Location = c("Vancouver Annacis Island",
                             "Vancouver Iona Island",
                             "Vancouver Lions Gate",
                             "Vancouver Lulu Island",
                             "Vancouver Northwest Langley",
                             "Edmonton Goldbar",
                             "Toronto Ashbridges Bay",
                             "Toronto Highland Creek",
                             "Toronto Humber",
                             "Toronto North Toronto",
                             "Halifax Dartmouth",
                             "Halifax Halifax",
                             "Halifax Millcove",
                             "Montreal North",
                             "Montreal South"),
                loc = c("VAI","VII","VLG","VLI","VNL","EGB","TAB","THC","THU","TNT",
                        "HDA","HHA","HMC","MMN","MMS"))

work_d <- phac_d %>% 
  filter(fractionid=="solid") %>% 
  left_join(df) %>% 
  filter(!is.na(loc)) %>% 
  mutate(Location = loc) %>% 
  dplyr::select(Date, Location, viral_load) %>% 
  rename("date"=Date) %>%
  mutate(source = "phac") %>% 
  rbind(modeldata %>% 
          group_by(Location,date) %>% 
          summarise(viral_load = mean(value.raw)) %>% 
          mutate(source = "paper")) %>% 
  group_by(Location, date) %>% 
  mutate(ll = length(viral_load)) %>%
  group_by(Location) %>% 
  mutate(max_date = max(date[which(source == "paper")])) %>% 
  mutate(viral_load = ifelse(source == "phac" & date <= max_date, NA_real_, viral_load)) %>% 
  arrange(Location, date) %>% 
  filter(!is.na(viral_load)) %>% 
  group_by(Location, date) %>% 
  mutate(ll = length(viral_load))

work_d <- work_d %>% 
  mutate(log10.value.raw = log10(viral_load)) %>% 
  mutate(log10.value.raw = ifelse(log10.value.raw < -1, NA_real_, log10.value.raw)) 

work_d <- work_d %>% 
  dplyr::select(- 'source', - 'max_date',-'ll') %>% 
  rename('value.raw'= viral_load)


modeldata <- work_d%>% 
  mutate(target = "N2") %>% 
  arrange(date, Location) %>% 
  as.data.table()



## Log full data
modeldata2 <- modeldata %>% 
  mutate(log10.value.raw = log(value.raw))


modeldata2 <- modeldata2 %>% 
  mutate(log10.value.raw = ifelse(log10.value.raw < -2, NA_real_,log10.value.raw)) %>% 
  mutate(value.raw = ifelse(is.na(log10.value.raw), NA_real_,value.raw))



# save(file="./Comparison work/Xiaotian et al/work_d_xiaotian_phac.RData", list="modeldata2")

date_range <-modeldata2$date %>% range()

modeldata2 <- data.frame(date = rep(seq(date_range[1], date_range[2],1),each=15),
                         Location = rep(unique(modeldata2$Location),length(seq(date_range[1], date_range[2],1)) )) %>% 
  left_join(data.frame(modeldata2), by = c("Location","date")) %>% 
  mutate(target = "N2") %>% 
  as.data.table(.)

dates_with_newdata <- modeldata2 %>% 
  filter(date >= ymd('2022-12-01') & date <= ymd('2023-03-31')) %>% 
  filter(!is.na(value.raw))%$% date %>% unique()


dates_totry = dates_with_newdata
# [seq(1, length(dates_with_newdata),7)]


for (i in 1:length(dates_totry)) {
  df <- modeldata2 %>% 
    filter(date <= dates_totry[i])
  
  
  trainingdata= df 
  testingdata <- modeldata2 %>% 
    filter(date > dates_totry[i] & date <= (dates_totry[i] + 7))
  
  set.seed(2)
  
  # modeldata = trainingdata
  # ID = c("Location", "target")
  # date = "date"
  # value = "log10.value.raw"
  # covariate = NULL
  # iteration = 5000
  # burnin = 2500
  # cores = 4
  # chains=4
  # nbasis = 50
  
  ## Log full data
  res <- WWmodel_meanconcentrate(modeldata = trainingdata, 
                                 ID = c("Location", "target"), 
                                 date = "date", 
                                 value = "log10.value.raw", 
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
                       value = "log10.value.raw", 
                       covariate = NULL, 
                       iteration = 5000, 
                       burnin = 2500,
                       nbasis = 50,
                       fpca.fit = res$fpca.fit)
  
  res_df <- res$df
  
  save(file=paste0("~/Wastewater/toEnglish/xiaotian_phac/nowcast_nowave_morebasis_xiaotian/model",i,".RData"), list = c("trainingdata","testingdata", "res_df","pred_df"))
}






