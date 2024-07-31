process_results <- function(df_full, tmbdat, samps1, polyOrder,  id_group, id_group_name){
  
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
    X = global_poly_helper(refined_x, p = p)
    X <- as.matrix(X[, 1:(p - degree)])
    for (i in 1:ncol(X)) {
      X[, i] <- (factorial(i + degree - 1)/factorial(i - 1)) * 
        X[, i]
    }
    B = as(local_poly_helper(knots, refined_x = refined_x, p = (p - 
                                                           degree)), "dgTMatrix")
    fitted_samps_deriv <- X %*% global_samps[(1 + degree):p, 
    ] + B %*% samps
    result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_deriv)))
    result
  }
  
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
  u_samps1 <- samps1$samps[(ncol(P) + ncol(X) + ncol(daily)+1):(ncol(P) + ncol(X) + ncol(daily)+ncol(obs)),]
  u_deriv_samps1 <- samps1$samps[(ncol(P) + ncol(X) + ncol(daily)+ncol(obs)+1):(nrow(samps1$samps)),]
  
  v <- compute_post_fun(samps = coefsamps1, global_samps = global_samps1, 
                         knots = knots, 
                         refined_x = df_full$t,
                         p = polyOrder, degree = 0)
  
  vderiv <- compute_post_fun(samps = coefsamps1, global_samps = global_samps1, 
                              knots = knots, 
                              refined_x = df_full$t,
                              p = polyOrder, degree = 1)
  
  Xfstat_full <-model.matrix(~site_id, 
                             data = df_full %>% mutate(site_id = factor(site_id)), 
                             contrasts.arg = list(site_id = "contr.sum"))[,-1]
  
  norm <- Xfstat_full
  norm[,2] <- -1
  norm[,1]<- -1
  
  v_u_fixed <- v[,-1]+ u_samps1 + Xfstat_full%*%Xf_samps1
  v_u_deriv <-vderiv[,-1]+ u_deriv_samps1
  v_u_fixed_deriv <-vderiv[,-1]+ u_deriv_samps1
  u_fixed <- u_samps1 + Xfstat_full%*%Xf_samps1
  v_fixed <- v[,-1] + Xfstat_full%*%Xf_samps1
  
  ## Nominal AR2 + ospline+ fixed effects
  df_full$exp_v_u_fixed <- as.numeric(apply(exp(v_u_fixed), MARGIN=1,median))
  df_full$exp_v_u_fixed_lwr <- as.numeric(apply(exp(v_u_fixed), MARGIN=1,quantile, 0.025))
  df_full$exp_v_u_fixed_upr <- as.numeric(apply(exp(v_u_fixed), MARGIN=1,quantile, 0.975))
  
  df_full$exp_v_u_fixed_deriv <- as.numeric(apply(exp(v_u_fixed)*(v_u_deriv), MARGIN=1,median))
  df_full$exp_v_u_fixed_deriv_lwr <- as.numeric(apply(exp(v_u_fixed)*(v_u_deriv), MARGIN=1,quantile, 0.025))
  df_full$exp_v_u_fixed_deriv_upr <- as.numeric(apply(exp(v_u_fixed)*(v_u_deriv), MARGIN=1,quantile, 0.975))
  
  
  ## Nominal AR2 + fixed effects
  df_full$exp_u_fixed <- as.numeric(apply(exp(u_fixed), MARGIN=1,median))
  df_full$exp_u_fixed_lwr <- as.numeric(apply(exp(u_fixed), MARGIN=1,quantile, 0.025))
  df_full$exp_u_fixed_upr <- as.numeric(apply(exp(u_fixed), MARGIN=1,quantile, 0.975))
  
  df_full$exp_u_fixed_deriv <- as.numeric(apply(exp(u_fixed)*(u_deriv_samps1), MARGIN=1,median))
  df_full$exp_u_fixed_deriv_lwr <- as.numeric(apply(exp(u_fixed)*(u_deriv_samps1), MARGIN=1,quantile, 0.025))
  df_full$exp_u_fixed_deriv_upr <- as.numeric(apply(exp(u_fixed)*(u_deriv_samps1), MARGIN=1,quantile, 0.975))
  
  ## Log AR2 + ospline+ fixed effects
  df_full$v_u_fixed <- as.numeric(apply(v_u_fixed, MARGIN=1,median))
  df_full$v_u_fixed_lwr <- as.numeric(apply(v_u_fixed, MARGIN=1,quantile, 0.025))
  df_full$v_u_fixed_upr <- as.numeric(apply(v_u_fixed, MARGIN=1,quantile, 0.975))
  
  df_full$v_u_fixed_deriv <- as.numeric(apply(v_u_deriv, MARGIN=1,median))
  df_full$v_u_fixed_deriv_lwr <- as.numeric(apply(v_u_deriv, MARGIN=1,quantile, 0.025))
  df_full$v_u_fixed_deriv_upr <- as.numeric(apply(v_u_deriv, MARGIN=1,quantile, 0.975))
  
  ## Log AR2 + fixed effects
  df_full$u_fixed <- as.numeric(apply(u_fixed, MARGIN=1,median))
  df_full$u_fixed_lwr <- as.numeric(apply(u_fixed, MARGIN=1,quantile, 0.025))
  df_full$u_fixed_upr <- as.numeric(apply(u_fixed, MARGIN=1,quantile, 0.975))
  
  df_full$u_fixed_deriv <- as.numeric(apply(u_deriv_samps1, MARGIN=1,median))
  df_full$u_fixed_deriv_lwr <- as.numeric(apply(u_deriv_samps1, MARGIN=1,quantile, 0.025))
  df_full$u_fixed_deriv_upr <- as.numeric(apply(u_deriv_samps1, MARGIN=1,quantile, 0.975))
  
  ## Log Ospline + fixed effects
  df_full$v_fixed <- as.numeric(apply(v_fixed, MARGIN=1,median))
  df_full$v_fixed_upr <- as.numeric(apply(v_fixed, MARGIN=1,quantile,p=0.975))
  df_full$v_fixed_lwr <- as.numeric(apply(v_fixed, MARGIN=1,quantile,p=0.025))
  
  df_full$exp_v_fixed_deriv <- as.numeric(apply(exp(v_fixed)*vderiv[,-1], MARGIN=1,median))
  df_full$exp_v_fixed_deriv_upr <- as.numeric(apply(exp(v_fixed)*vderiv[,-1], MARGIN=1,quantile, p=0.975))
  df_full$exp_v_fixed_deriv_lwr<- as.numeric(apply(exp(v_fixed)*vderiv[,-1], MARGIN=1,quantile, p=0.025))
  
  ## Log Ospline 
  df_full$v <- as.numeric(apply(v[,-1], MARGIN=1,median))
  df_full$v_upr <- as.numeric(apply(v[,-1], MARGIN=1,quantile,p=0.975))
  df_full$v_lwr <- as.numeric(apply(v[,-1], MARGIN=1,quantile,p=0.025))
  
  df_full$v_deriv <- as.numeric(apply(vderiv[,-1], MARGIN=1,median))
  df_full$v_deriv_upr <- as.numeric(apply(vderiv[,-1], MARGIN=1,quantile, p=0.975))
  df_full$v_deriv_lwr<- as.numeric(apply(vderiv[,-1], MARGIN=1,quantile, p=0.025))
  
  ## Ospline
  df_full$exp_v <- as.numeric(apply(exp(v[,-1]), MARGIN=1,median))
  df_full$exp_v_upr <- as.numeric(apply(exp(v[,-1]), MARGIN=1,quantile, p = 0.975))
  df_full$exp_v_lwr <- as.numeric(apply(exp(v[,-1]), MARGIN=1,quantile, p = 0.025))
  
  df_full$exp_v_deriv <- as.numeric(apply((exp(v[,-1]) * vderiv[,-1]), MARGIN=1,median))
  df_full$exp_v_deriv_upr <- as.numeric(apply((exp(v[,-1]) * vderiv[,-1]), MARGIN=1,quantile, 0.975))
  df_full$exp_v_deriv_lwr<- as.numeric(apply((exp(v[,-1]) * vderiv[,-1]), MARGIN=1,quantile, 0.025))
  
if ( id_group == 1){
  post_samps_df <- df_full %>% 
    dplyr::select('sample_date','site_id') %>% 
    cbind(as.data.frame(v_u_fixed)) %>% 
    melt(id.vars = 1:2)
  
  post_samps_df_deriv <- df_full %>% 
    dplyr::select('sample_date','site_id') %>% 
    cbind(as.data.frame(v_u_fixed_deriv)) %>% 
    melt(id.vars = 1:2)
  
  output_postsamps = post_samps_df %>% 
    rename("v_u_fixed" = value) %>% 
    cbind(post_samps_df_deriv%>% dplyr::select(value) %>% rename("v_u_fixed_deriv"= value) ) 
  
  tmp<- output_postsamps  %>% 
    group_by(variable, sample_date) %>% 
    summarise(ave_exps = mean(exp(v_u_fixed)),
              ave_exps_deriv = mean(v_u_fixed_deriv*exp(v_u_fixed))) %>% 
    ungroup() %>% 
    group_by(sample_date) %>% 
    summarise(ave_exp_v_u_fixed = median(ave_exps),
              ave_exp_v_u_fixed_upr = quantile(ave_exps, 0.975),
              ave_exp_v_u_fixed_lwr = quantile(ave_exps, 0.025),
              ave_exp_v_u_fixed_deriv = median(ave_exps_deriv),
              ave_exp_v_u_fixed_deriv_upr = quantile(ave_exps_deriv, 0.975),
              ave_exp_v_u_fixed_deriv_lwr = quantile(ave_exps_deriv, 0.025),
              post_prob_ave_exp_v_u_fixed_deriv = length(which(ave_exps_deriv>0))/length(ave_exps_deriv))}
  
  return(list(df_full=df_full, station_ave_df = tmp))
}