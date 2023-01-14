# Simulate data under c2 type pf censoring (noninformative censoring on the residual scale)
# and then estimate using the DR, IPW.Q, IPW.T1, IPW.T2 with coxph for estimating the nuisance parameters, 
# plus the naive estimator and the full data estimator
c.simu_estimate.cox_naive_full <- function(seed, nu, covariates.T, covariates.Q, 
                                            trim.C, n.boot,
                                            n, multi, Tmodel, Qmodel, T.min, Q.min, Q.max, tau,
                                            beta0.T, beta0.Q, beta.T, beta.Q, beta.T2, beta.Q2,
                                            shape.T, shape.Q, epsT.sd, epsQ.sd, 
                                            C.min, shape.C, scale.C, type.C){
  
  # simulate data
  set.seed(seed)
  gen_result = gen(n, multi, Tmodel, Qmodel, T.min, Q.min, Q.max, tau,
                   beta0.T, beta0.Q, beta.T, beta.Q, beta.T2, beta.Q2,
                   shape.T, shape.Q, epsT.sd, epsQ.sd, C.min, shape.C, scale.C, type.C)
  dat = gen_result$dat    # the observed data
  dat.full = gen_result$dat.full    # full data that corresponds to the observed sample with sample size n 
  
  # Add higher order and interaction terms to dat if necessary
  dat$Z1sq = (dat$Z1)^2
  dat$Z1Z2 = dat$Z1 * dat$Z2
  
  
  ### DR estimator with coxph to fit the two parameters
  # DR estimator based on Cox models without cross fitting
  est.cox = c.est_DR.cox(dat, nu, covariates.T, covariates.Q, trim.C, type.C)
  
  estimate.cox = c(DR = est.cox$est_DR, 
                   IPQW = est.cox$est_IPQW, 
                   IPSW = est.cox$est_IPSW, 
                   IPSW2 = est.cox$est_IPSW2)
  se.cox = c(se_DR = est.cox$se_DR,
             se_IPQW = est.cox$se_IPQW,
             se_IPSW = est.cox$se_IPSW,
             se_IPSW2 = est.cox$se_IPSW2)
  beta = c(beta = est.cox$beta, beta2 = est.cox$beta2)
  coxcoef.T = est.cox$beta.T
  coxcoef.Q2 = est.cox$beta.Q2
  
  
  # SE and CI constructed using boostrap
  est_boot.cox = boot(dat, c.estimate.cox_boot, R = n.boot, nu = nu, 
                      covariates.T = covariates.T, covariates.Q = covariates.Q, 
                      trim.C = trim.C, type.C = type.C)
  se_boot.cox = apply(est_boot.cox$t, MARGIN = 2, sd)
  beta_se_boot.cox = se_boot.cox[c(5,6)]
  se_boot.cox = se_boot.cox[-c(5,6)]
  names(se_boot.cox) = paste("se_boot_", names(estimate.cox), sep = "")
  
  tab.cox = c(estimate.cox, se.cox, se_boot.cox) 
  
  
  ### naive estimator - ignore left truncation, IPCW to deal with censoring
  n = nrow(dat)
  if(sum(1-dat$delta) == 0){
    Sd = function(x){return(1)}
  }else{
    kmfit.D = survfit(Surv(X-Q, 1-delta)~1, data = dat, type = "kaplan-meier")
    Sd = stepfun(kmfit.D$time,  pmax(c(1, kmfit.D$surv), trim.C) ) # Bound Sc away from 0 so we do not get extreme weights
  }
  est.naive = sum(dat$delta*nu(dat$X)/Sd(dat$X))/sum(dat$delta/Sd(dat$X))
  se.naive = NA
  names(est.naive) = "est.naive"
  names(se.naive) = "se.naive"
  
  # SE and CI constructed using boostrap
  est_boot.naive = boot(dat, c.estimator.naive_boot, R = n.boot, nu = nu, trim.C = trim.C)
  se_boot.naive = apply(est_boot.naive$t, MARGIN = 2, sd)
  names(se_boot.naive) = "se_boot_naive"
  
  tab.naive = c(est.naive, se.naive, se_boot.naive)  
  
  
  ### full data estimator
  n.full = nrow(dat.full)
  est.full = mean(nu(dat.full$TT))
  se.full = sd(nu(dat.full$TT))/sqrt(n.full)
  names(est.full) = "est.full"
  names(se.full) = "se.full"
  
  # SE and CI constructed using boostrap
  est_boot.full = boot(dat.full, estimator.full_boot, R = n.boot, nu = nu)
  se_boot.full = apply(est_boot.full$t, MARGIN = 2, sd)
  names(se_boot.full) = "se_boot_full"
  
  tab.full = c(est.full, se.full, se_boot.full) 
  
  
  
  return(list(tab.cox = tab.cox, 
              tab.naive = tab.naive,
              tab.full = tab.full,
              beta = beta, 
              seed = seed,
              coxcoef.T = coxcoef.T,
              coxcoef.Q2 = coxcoef.Q2))
}


c.estimate.cox_boot <- function(dat, id, nu, covariates.T, covariates.Q, trim.C, type.C){
  dat = dat[id,]
  result.cox = c.est_DR.cox(dat, nu, covariates.T, covariates.Q, trim.C, type.C)
  est.cox = c(DR = result.cox$est_DR, 
              IPQW = result.cox$est_IPQW, 
              IPSW = result.cox$est_IPSW, 
              IPSW2 = result.cox$est_IPSW2,
              beta = result.cox$beta,
              beta2 = result.cox$beta2)
  
  return(est.cox)
}


c.estimator.naive_boot <- function(dat, id, nu, trim.C){
  dat = dat[id,]
  if(sum(1-dat$delta) == 0){
    Sd = function(x){return(1)}
  }else{
    kmfit.D = survfit(Surv(X-Q, 1-delta)~1, data = dat, type = "kaplan-meier")
    Sd = stepfun(kmfit.D$time,  pmax(c(1, kmfit.D$surv), trim.C) ) # Bound Sc away from 0 so we do not get extreme weights
  }
  est.naive =  sum(dat$delta*nu(dat$X)/Sd(dat$X))/sum(dat$delta/Sd(dat$X))
  
  return(est.naive)
}

estimator.full_boot <- function(dat.full, id, nu){
  dat.full = dat.full[id,]
  est.full = mean(nu(dat.full$TT))
  
  return(est.full)
}




