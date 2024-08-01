# R version 4.4.1 (2024-06-14)
# Platform: aarch64-apple-darwin20
# Running under: macOS Monterey 12.4

library(readr) # readr_2.1.5
library(dplyr) #dplyr_1.1.4
library(data.table) #data.table_1.15.4

##########################################
#### Data preprocess
##########################################

## This data is PHAC data that was provided to me by the authors of
## https://www.nature.com/articles/s41598-022-17543-y. This is the data the authors used to 
## develop their statistical framework for SARS-CoV-2 wastewater concentrations. This is also the 
## paper that we compare our framework with in Section 4.1

raw_d <- readr::read_rds(file = './data/ww-db-2021-11-18.rds')

RawData = as.data.table(raw_d)
RawData<- RawData %>% 
   dplyr::select("Location",  "date", "method", "target", "replicate", "labname", "value.raw","uqnd")

#### DataPrep is a function that the authors created to prepare the data for analysis.
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

modeldata <- DataPrep(RawData, "N2") ## prepares the N2 data

data_orig <- range(modeldata$date) # saves the range of the paper's data

# Import public up-to-date PHAC 
# This is data I downloaded from https://health-infobase.canada.ca/wastewater/
# It's no longer available in this format. 
phac_d <- read_csv("./data/covid19-wastewater.csv")

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


### This is to join the paper's PHAC data with the newer public PHAC data.

work_d <- phac_d %>% 
  filter(fractionid=="solid") %>% 
  left_join(df) %>% 
  filter(!is.na(loc)) %>% # so I am analyzing the same stations as the paper. 
  mutate(Location = loc) %>% 
  dplyr::select(Date, Location, viral_load) %>% 
  rename("date"=Date) %>%
  mutate(source = "phac") %>% 
  rbind(modeldata %>% 
              group_by(Location,date) %>% 
              summarise(viral_load = mean(value.raw)) %>% # the public phac data is the average of 2 N2 replicates, so did this to match the two datasets
              mutate(source = "paper")) %>% 
  group_by(Location, date) %>% 
  mutate(ll = length(viral_load)) %>%
  group_by(Location) %>% 
  mutate(max_date = max(date[which(source == "paper")])) 

work_d %>% group_by(source) %>% 
  summarise(m1 = min(date),
            m2 = max(date))

## If there is any overlap between the new PHAC and the paper's PHAC, want to use the paper's data. 
## It was mostly the same when there was overlap.
## Slight differences because I think the paper had access to details about censoring, or 
## unreliable values and were able to remove those. 

work_d <- work_d %>% 
  mutate(viral_load = ifelse(source == "phac" & date <= max_date, NA_real_, viral_load)) %>% 
  arrange(Location, date) %>% 
  filter(!is.na(viral_load)) %>% 
  group_by(Location, date) %>% 
  mutate(ll = length(viral_load))

work_d <- work_d %>% 
  dplyr::select(- 'source', - 'max_date',-'ll') %>% 
  rename('value'= viral_load)


work_d <- work_d%>% 
  mutate(target = "N2") %>% 
  arrange(date, Location) %>% 
  as.data.table()

# I believe this is similar to the data cleaning that done in the original paper
# In any case, this step removes only 6 observations (0.153%)
work_d <- work_d %>% 
  mutate(value = ifelse( log(value) < -2, NA,value)) 

work_d <- work_d %>% 
  rename("sample_date" = date, 
         "site_id" =  Location)

save(file="./data/work_d.RData", list="work_d")