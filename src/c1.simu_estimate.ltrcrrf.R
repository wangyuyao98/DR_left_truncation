# Simulate data and then estimate using the cf, IPW.Q, IPW.T1, IPW.T2 with  LTRCRRF:forests for estimating the nuisance parameters
# alpha: (1-alpha)*100% CI is computed 
c1.simu_estimate.ltrcrrf <- function(seed, nu, covariates.X, covariates.Q, 
                             K, mtry, ntree, trim, trim.C){
  
  # simulate data
  set.seed(seed)
  gen_result = gen(n, multi, Tmodel, Qmodel, T.min, Q.min, Q.max, tau,
                   beta0.T, beta0.Q, beta.T, beta.Q, beta.T2, beta.Q2,
                   shape.T, shape.Q, epsT.sd, epsQ.sd, C.min, shape.C, scale.C)
  dat = gen_result$dat    # the observed data
  dat.full = gen_result$dat.full    # full data that corresponds to the observed sample with sample size n 
  
  
  ### cf estimator with LTRCforests::ltrcrrf to fit the two parameters
  # cross-fitting is applied
  est.ltrcrrf = c1.est_cf.ltrcrrf(dat, nu, covariates.X, covariates.Q, K,
                                  mtry, ntree, trim, trim.C)
  estimate.ltrcrrf = c(cf = est.ltrcrrf$est_DR, 
                        IPW.Q = est.ltrcrrf$est_IPQW, 
                        Reg.T1 = est.ltrcrrf$est_IPSW, 
                        Reg.T2 = est.ltrcrrf$est_IPSW2)
  se.ltrcrrf = c(se_cf = est.ltrcrrf$se_DR,
                 se_IPW.Q = est.ltrcrrf$se_IPQW,
                 se_Reg.T1 = est.ltrcrrf$se_IPSW,
                 se_Reg.T2 = est.ltrcrrf$se_IPSW2)
  
  beta = c(beta_hat.Q = est.ltrcrrf$beta, beta_hat.T = est.ltrcrrf$beta2)
  
  se_boot.ltrcrrf = rep(NA, 4)
  names(se_boot.ltrcrrf) = c("bootse.cf", "bootse.IPW.Q", "bootse.Reg.T1", "bootse.Reg.T2")

  tab.ltrcrrf = c(estimate.ltrcrrf, se.ltrcrrf, se_boot.ltrcrrf)
  
  
  return(list(tab.ltrcrrf = tab.ltrcrrf,
              beta = beta, 
              seed = seed))
}


