
fcst_trend = function(h.ahead, modeldata, model_res, ID,
                      date, value, covariate = NULL, iteration, burnin, nbasis,fpca.fit) {
  # --- reshape input data
  Ymat = reshape(modeldata[, names(modeldata) %in% c(ID, date, value), with = FALSE],
                 idvar = ID, timevar = date, direction = "wide")
  Ymat = Ymat[order(Location)]
  IDdt = Ymat[, names(Ymat) %in% ID, with = FALSE]
  Ymat = Ymat[, !ID, with = FALSE]
  names(Ymat) = gsub(".*[.]", "", names(Ymat))
  
  I = dim(Ymat)[1] 
  T = dim(Ymat)[2]
  

  Ymat_mu = sweep(as.matrix(Ymat), 2, apply(Ymat, 2, function(x) mean(x, na.rm = TRUE)))
  # fpca.fit = fpca.sc(as.matrix(Ymat), pve= 0.99,nbasis = nbasis, var = TRUE, 
  #                    simul = TRUE)
  
  # print(fpca.fit$npc)
  
  PHI = t(fpca.fit$efunctions)
  L0 = dim(PHI)[1]
  Lambda = sqrt(fpca.fit$evalues)
  Ymat_random = as.data.table(matrix(rep(fpca.fit$mu, I), nrow = I, byrow = TRUE))
  names(Ymat_random) = names(Ymat)
  Ymat_fix = Ymat - Ymat_random
  
  
  
  Yhat_pred = as.numeric(forecast::forecast(
    forecast::auto.arima(fpca.fit$mu, d = 2),
    h = h.ahead)$mean)  ## d is the order of first-differencing
  PHIpred = t(apply(PHI, 1, function(x) as.numeric(forecast::forecast(
    forecast::auto.arima(x, d = 2), h = h.ahead)$mean)))
  Xpred = list()
  XPHIpred = list()

  list_of_draws = rstan::extract(model_res$fit)
  after = iteration - burnin + 1
  Ypred = list()
  for (i in 1:I) {
    Ypred[[i]] = matrix(rep(Yhat_pred, after), nrow = after, byrow = TRUE)
    for (l in 1:L0) {
      Ypred[[i]] = Ypred[[i]] +
        (list_of_draws[[paste("beta_xi", l, sep = "")]])[burnin:iteration, i] %*% t(PHIpred[l, ])
    }
  }
  
  for (i in 1:I) {
    Ypred[[i]] = matrix(rep(Yhat_pred, after), nrow = after, byrow = TRUE)
    for (l in 1:L0) {
      Ypred[[i]] = Ypred[[i]] + (list_of_draws[[paste("beta_xi", 
                                                    l, sep = "")]])[burnin:iteration, i] %*% t(PHIpred[l, 
                                                    ])
    }
  }
  
  loc.names = IDdt$Location[(1:I)] 
  
  tmp = lapply(Ypred, as.data.frame) #makes every element in a list a data frame
  for(i in seq_along(Ypred)) tmp[[i]]$Location <- loc.names[i]  #adds a column for the site
  
  df = do.call('rbind', tmp) %>%
    tidyr::pivot_longer(-Location) %>%
    mutate(name = as.numeric(stringr::str_remove(name, '^V'))) %>%
    rename(time.ahead = name)
  date_max = max(modeldata$date)
  df$date = date_max + df$time.ahead
  
  ## add the extra variation due to sigma
  set.seed(2)
  df=df %>% 
    group_by(Location, time.ahead) %>% 
    mutate(vv = 1:2501) %>% 
    mutate(pred = rnorm(after, value, sqrt(list_of_draws[["sigma"]][burnin:iteration])),
           sigma2 = sqrt(list_of_draws[["sigma"]][burnin:iteration]))
  
  df_mu = data.frame(date = seq(min(df$date), max(df$date),1)) %>% 
    mutate(mu_total =Yhat_pred)
  
  df <- df %>% 
    left_join(df_mu, by = "date")
  
  return(df)
}
