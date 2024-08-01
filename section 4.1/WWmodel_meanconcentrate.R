WWmodel_meanconcentrate = function(modeldata,
                                   ID,
                                   date,
                                   value,
                                   covariate = NULL,
                                   iteration,
                                   burnin,
                                   cores = 1,
                                   chains = 4,
                                   nbasis = nbasis) {
  Ymat = reshape(modeldata[, names(modeldata) %in% c(ID, date, value), with = FALSE],
                 idvar = ID, timevar = date, direction = "wide")
  Ymat = Ymat[order(Location)]
  IDdt = Ymat[, names(Ymat) %in% ID, with = FALSE]
  Ymat = Ymat[, !ID, with = FALSE]
  names(Ymat) = gsub(".*[.]", "", names(Ymat))
  I = dim(Ymat)[1] 
  T = dim(Ymat)[2]
  
  Ymat_mu = sweep(as.matrix(Ymat), 2, apply(Ymat, 2, function(x) mean(x, na.rm = TRUE)))
  fpca.fit = fpca.sc(as.matrix(Ymat), pve= 0.99,nbasis = nbasis, var = TRUE, 
                     simul = TRUE)
  
  print(fpca.fit$npc)
  
  PHI = t(fpca.fit$efunctions)
  L0 = dim(PHI)[1]
  Lambda = sqrt(fpca.fit$evalues)
  Ymat_random = as.data.table(matrix(rep(fpca.fit$mu, I), 
                                     nrow = I, byrow = TRUE))
  names(Ymat_random) = names(Ymat)
  Ymat_fix = Ymat - Ymat_random
  Y = as.numeric(t(Ymat_fix))
  missing = which(is.na(Y))
  Y = Y[-missing]
  
  XImat = list()
  for (l in 1:L0) {
    XImat[[l]] = matrix(0, nrow = I * T, ncol = I)
    for (i in 1:I) XImat[[l]][((i-1) * T + 1):(( i) * T), i] = PHI[l, ]
    XImat[[l]] = XImat[[l]][-missing, ]
    if (I == 1) 
      XImat[[l]] = matrix(XImat[[l]], ncol = 1)
  }
  names(XImat) = paste("XI", 1:L0, sep = "")
  if (is.null(covariate)) {
    data = c(list(n = length(Y), I = I, L = L0), XImat, list(Y = Y))
    stancode = createStanModel(L0, 0)
    regression_model = stan_model(model_code = stancode, verbose = TRUE)
  }
  
  fit = rstan::sampling(regression_model, data = data, chains = chains, 
                        iter = iteration, cores = cores, refresh = 100)
  list_of_draws = rstan::extract(fit)
  after = iteration - burnin + 1
  Yhat = list()
  for (i in 1:I) {
    Yhat[[i]] = matrix(rep(unlist(Ymat_random[i, ]), 
                           after), nrow = after, byrow = TRUE)
    for (l in 1:L0) {
      Yhat[[i]] = Yhat[[i]] + (list_of_draws[[paste("beta_xi", 
                                                    l, sep = "")]])[burnin:iteration, i] %*% t(PHI[l, 
                                                    ])
    }
    if (!is.null(covariate)) {
      for (p in 1:P) {
        Yhat[[i]] = Yhat[[i]] + (list_of_draws[[paste("beta", 
                                                      p, sep = "")]])[burnin:iteration, ] %*% XPHI[[p]] * 
          matrix(rep(unlist(X[[1]][2 * i, ]), after), 
                 nrow = after, byrow = TRUE)
      }
    }
  }
  
  
  loc.names = IDdt$Location[(1:I)] 
  
  tmp = lapply(Yhat, as.data.frame) #makes every element in a list a data frame
  for(i in seq_along(Yhat)) tmp[[i]]$Location <- loc.names[i]  #adds a column for the site
  
  df = do.call('rbind', tmp) %>%
    tidyr::pivot_longer(-Location) %>%
    mutate(name = as.numeric(stringr::str_remove(name, '^V'))) %>%
    rename(date = name)
  date_min = min(modeldata$date)
  df$date = date_min + df$date - 1
  
  ## add the extra variation due to sigma
  set.seed(2)
  df=df %>% 
    group_by(Location, date) %>% 
    mutate(vv = 1:2501) %>% 
    mutate(pred = rnorm(after, value, sqrt(list_of_draws[["sigma"]][burnin:iteration])),
           sigma2 = list_of_draws[["sigma"]][burnin:iteration])
  
  df_mu = data.frame(date = seq(min(df$date), max(df$date),1)) %>% 
    mutate(mu_total =unlist(Ymat_random[i, ]))
  
  df <- df %>% 
    left_join(df_mu, by = "date")
  
  return(list(fit = fit, Yhat = Yhat, df=df,fpca.fit=fpca.fit))
}

createStanModel <- function (L0, P, Lp = NULL) 
{
  scode = "data {\n  int<lower=0> n;\n  int<lower=0> I;\n  int<lower=0> L; \n"
  if (P > 0) {
    for (p in 1:P) {
      scode = paste(scode, paste("int<lower=0> L", p, ";", 
                                 sep = ""), "\n")
    }
    for (p in 1:P) {
      scode = paste(scode, paste("matrix[n,L", p, "] X", 
                                 p, ";", sep = ""), "\n")
    }
  }
  for (l in 1:L0) {
    scode = paste(scode, paste("matrix[n,I] XI", l, ";", 
                               sep = ""), "\n")
  }
  scode = paste(scode, "vector[n] Y; \n} \n parameters { \n")
  if (P > 0) {
    for (p in 1:P) {
      scode = paste(scode, paste("vector[L", p, "] beta", 
                                 p, ";", sep = ""), "\n", paste("real<lower=0> tau", 
                                                                p, ";", sep = ""), "\n")
    }
  }
  scode = paste(scode, "real<lower=0> delta; \n")
  for (l in 1:L0) {
    scode = paste(scode, paste("vector[I] beta_xi", l, ";", 
                               sep = ""), "\n")
  }
  for (l in 1:L0) {
    scode = paste(scode, paste("real<lower=0> lambda", l, 
                               ";", sep = ""), "\n")
  }
  scode = paste(scode, "real<lower=0> sigma; \n} \n model { \n delta ~ gamma(2, 100); \n")
  if (P > 0) {
    for (p in 1:P) {
      scode = paste(scode, paste("tau", p, " ~ gamma((L", 
                                 p, "+1)/2, delta/2); \n", sep = ""))
      for (l in 1:Lp[p]) {
        scode = paste(scode, paste("beta", p, "[", l, 
                                   "] ~ normal(0, sqrt(tau", p, ")*sigma); \n", 
                                   sep = ""))
      }
    }
  }
  for (l in 1:L0) {
    scode = paste(scode, paste("beta_xi", l, " ~ normal(0, sqrt(lambda", 
                               l, ")); \n", sep = ""))
    scode = paste(scode, paste("lambda", l, " ~ inv_gamma(0.1, 0.1); \n", 
                               sep = ""))
  }
  scode = paste(scode, "sigma ~ inv_gamma(0.1, 0.1); \n")
  model = "Y ~ normal("
  if (P > 0) {
    for (p in 1:P) {
      model = paste(model, "X", p, "*beta", p, " + ", sep = "")
    }
  }
  for (l in 1:L0) {
    model = paste(model, "XI", l, "*beta_xi", l, sep = "")
    if (l < L0) 
      model = paste(model, " + ", sep = "")
  }
  model = paste(model, ", sqrt(sigma)); \n", sep = "")
  scode = paste(scode, model, "}")
  return(scode)
}
