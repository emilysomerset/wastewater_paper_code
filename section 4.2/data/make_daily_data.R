make_daily_data <- function (dat, analysis_unit, covars, LOD, syears = NULL, smonths = NULL, 
          construct_catchment_covariates) 
{

  ww.long <- ww_wide_to_long(dat, analysis_unit$space, show.warning = TRUE) # wide to long format
  site_geo <- get_sites_geo(ww.long$site_uniq) # eastings, northings, latitutde, longitude
  ww.with.geo <- merge_ww_with_site_geo(ww.long$ww.long, site_geo, 
                                        by = analysis_unit$space) # merge wastewater with site 
  
  syears  <- NULL #2021  #  which year of data to be modelled
  smonths <- NULL #6:9
  
  ww <- select_region_time(ww.with.geo, sel.year = syears, 
                           sel.month = smonths)  # can leave null 
  
  
  # ww <- add_time_site_indices(ww, analysis_unit)  # adds year and week index 
  
  ww <- ww %>% 
    mutate(Date = dmy(Date)) %>% 
    mutate(time_index = as.numeric(as.factor(Date))) %>% 
    mutate(site_index = as.numeric(as.factor(uwwName)))
  
  #  LOD threshold
  L <- 160
  LOD <- L
  # values below the detection limit are set to 80
  ww <- replace_LOD(ww, raw_gc_column = "gc", LOD.means = LOD/2, 
                    LOD.sds = NULL, dist = "Uniform", isim_LOD = NULL, LOD_holder = "tLOD")
  
  # takes the mean of the logs of the concentrates for each week
  # ww %>% group_by(time_index_yw,uwwCode) %>% summarise(mm = mean(log(gc), na.rm=TRUE)) %>% head(2) %>% colnames()
  ww.daily <- ww
  if (!is.null(covars)) {
    X_names <- covars$covariates[covars$aggregation != "none"]
    op <- covars$aggregation[covars$aggregation != "none"]
    if (construct_catchment_covariates) {
      catchments_lookup <- get_catchments_lookup(ww.weekly)
      covariates.data <- LSOA.to.catchment.covariates(X_names, 
                                                      op, catchments_lookup)
    }
    else {
      data_file <- system.file("extdata", "catchment_covariates_no_genomic.rds", 
                               package = "publicWW")
      covariates.data <- readRDS(data_file)
    }
    ww.weekly <- merge_covariates_ww_data(ww.weekly, covariates.data, 
                                          X_names)
    covars_time_dependent <- covars[covars$aggregation == 
                                      "none", ]
    if (nrow(covars_time_dependent) != 0) {
      genomics.weekly <- get_weekly_genomics_data(ww, syears, 
                                                  smonths, analysis_unit, covars_time_dependent$covariates)
      ww.weekly <- dplyr::left_join(ww.weekly, genomics.weekly, 
                                    by = c("time_index"), all.x = TRUE)
    }
  }
  
  ww.daily$region <- as.factor(ww.daily$RegionName)
  ww.daily$region_index <- add_region_index(ww.daily$RegionName)
  if (!is.null(covars)) {
    ww.weekly <- scale_covariates(ww.weekly, covars$covariates)
  }
  
  ww.daily <- ww.daily %>% 
    mutate(log_e_gc = log(gc))
  return(ww.daily)
}
