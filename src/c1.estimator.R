# function for computing the dr and cf estimator
# with IPCW weights to handle independent censoring


# dr-Cox-Cox ------------------------------------------------------------------
# dr, IPW.Q, Reg.T1, Reg.T2 estimators that use Cox models for fitting T~Z and (tau-Q)~Z models 
# The Cox models are fitted on the entire sample
# The argument `dat' should contains X := min(T,C), Q, delta := I(T<C), and the covariates
c1.est_DR.cox <- function(dat, nu, covariates.X, covariates.Q, trim.C = 1e-7){
    n = nrow(dat)
    
    # compute the survival curve Sc for C from KM estimator
    # If no censoring, Sc is the function that is constant 1
    if(sum(1-dat$delta) == 0){
        Sc = function(x){return(1)}
    }else{
        kmfit.C = survfit(Surv(Q, X, 1-delta)~1, data = dat, type = "kaplan-meier")
        Sc = stepfun(kmfit.C$time,  pmax(c(1, kmfit.C$surv), trim.C) )  # Bound Sc away from 0 so we do not get extreme weights
    }

    tau = max(dat$X)+1
    dat$delta.1 = rep(1, nrow(dat))
    formula.X = formula(paste("Surv(Q, X, delta.1) ~ ", 
                              paste(covariates.X, collapse = "+"), collapse = ""))
    formula.Q2 = formula(paste("Surv(tau-X, tau-Q, delta.1) ~ ", 
                               paste(covariates.Q, collapse = "+"), collapse = ""))
    Z.X = as.matrix(dat[,covariates.X])
    Z.Q = as.matrix(dat[,covariates.Q])
    
    # remove the subjects with T-Q < 10^{-7} when fitting Cox models
    # such subjects will trigger error when fitting a Cox model
    id.cox = ((dat$X - dat$Q) >= 10^{-7})
    dat.cox = dat[id.cox, ]
    
    # fit Cox-PH model for T|Z
    fit.X = coxph(formula.X, data = dat.cox)
    # fit Cox-PH model for (tau-Q)|Z
    fit.Q2 = coxph(formula.Q2, data = dat.cox)
    basehaz.X = basehaz(fit.X, centered = FALSE)
    basehaz.Q2 = basehaz(fit.Q2, centered = FALSE)
    beta.X = coef(fit.X)
    beta.Q2 = coef(fit.Q2)
    X = dat$X
    Q = dat$Q
    Q2 = tau - Q
    delta = dat$delta
    
    Gtz = (baseS(tau-X, basehaz.Q2))^exp(Z.Q %*% beta.Q2)
    Gqz = (baseS(tau-Q, basehaz.Q2))^exp(Z.Q %*% beta.Q2)
    Fqz = 1-(baseS(Q, basehaz.X))^exp(Z.X %*% beta.X)
    mqz = diag(integral_F(Q, nu, delta, Sc, Z.X, basehaz.X, beta.X))
    
    # Define a constant 1 function nu1
    nu1 <- function(x){
      return(rep(1, length(x)))
    }
    m1qz = diag(integral_F(Q, nu1, delta, Sc, Z.X, basehaz.X, beta.X))
    
    # DR estimator
    # compute the denominator of the estimator
    DDen1 = delta/(Sc(X)*Gtz)
    DDen2 = m1qz/(Gqz*(1-Fqz)) 
    DDen3 = integral_JG_m(Q, X, delta, Sc, Z.Q, Z.X, 
                          basehaz.Q2, basehaz.X, beta.Q2, beta.X, tau, nu1)
    
    den1 = sum(DDen1)
    den2 = sum(DDen2)
    den3 = sum(DDen3)
    
    # compute the numerator of the estimator
    NNum1 = delta*nu(X)/(Sc(X)*Gtz)
    NNum2 = mqz/(Gqz*(1-Fqz))
    NNum3 = integral_JG_m(Q, X, delta, Sc, Z.Q, Z.X, 
                          basehaz.Q2, basehaz.X, beta.Q2, beta.X, tau, nu)
    
    num1 = sum(NNum1)
    num2 = sum(NNum2)
    num3 = sum(NNum3)
    
    # The estimator based on EIF
    est_DR = (num1+num2-num3)/(den1+den2-den3)
    beta = 1/(den1/n)   # truncation rate estimator based on Q model
    IF_DR = beta*((NNum1+NNum2-NNum3) - est_DR*(DDen1+DDen2-DDen3))
    se_DR = sqrt(mean(IF_DR^2))/sqrt(n)
    
    
    # The IPQW estimator
    est_IPQW = num1/den1
    SF_IPQW = NNum1 - est_IPQW*DDen1
    se_IPQW = sqrt(mean(SF_IPQW^2)/(den1/n)^2)/sqrt(n)
    
    
    # The IPSW estimator
    Enutz = as.vector(integral_F(tau, nu, delta, Sc, Z.X, basehaz.X, beta.X))
    Enu1tz = as.vector(integral_F(tau, nu1, delta, Sc, Z.X, basehaz.X, beta.X))
    
    DDen4_1 = delta/Sc(X) + m1qz/(1-Fqz)
    DDen4 = Enu1tz/(1-Fqz)
    
    NNum4 = delta*nu(X)/Sc(X) + mqz/(1-Fqz)
    NNum5 = Enutz/(1-Fqz)
    
    den4_1 = sum(DDen4_1)
    den4 = sum(DDen4)
    num4 = sum(NNum4)
    num5 = sum(NNum5)
    
    est_IPSW = num4/den4_1
    est_IPSW2 = num5/den4
    SF_IPSW = NNum4 - est_IPSW*DDen4
    se_IPSW = sqrt(mean(SF_IPSW^2)/(den4/n)^2)/sqrt(n)
    SF_IPSW2 = NNum5 - est_IPSW2*DDen4
    se_IPSW2 = sqrt(mean(SF_IPSW2^2)/(den4/n)^2)/sqrt(n)
    
    
    # estimate truncation rate
    beta2 = 1/(den4/n)   # estimator based on T model
    
    return(list(est_DR = est_DR, 
                est_IPQW = est_IPQW,
                est_IPSW = est_IPSW, est_IPSW2 = est_IPSW2,
                beta.X = beta.X, beta.Q2 = beta.Q2,
                se_DR = se_DR,
                se_IPQW = se_IPQW,
                se_IPSW = se_IPSW,
                se_IPSW2 = se_IPSW2,
                beta = beta,
                beta2 = beta2))
    
}



# function that take input t0 (posibaly a vector) and estimate the survival probability(s) at t0
# with dr-Cox-Cox, IPW.Q-Cox, Reg.T1, Reg.T2 estimators
c1.estSurv_DR.cox <- function(t0, dat, covariates.X, covariates.Q, trim.C = 1e-7){
  t0.list = t0
  m = length(t0.list)
  n = nrow(dat)
  
  tau = max(dat$X)+1
  dat$delta.1 = rep(1, nrow(dat))
  
  # Define a constant 1 function nu1
  nu1 <- function(x){
    return(rep(1, length(x)))
  }
  
  # compute the survival curve Sc for C from KM estimator
  # If no censoring, Sc is the function that is constant 1
  if(sum(1-dat$delta) == 0){
    Sc = function(x){return(1)}
  }else{
    kmfit.C = survfit(Surv(Q, X, 1-delta)~1, data = dat, type = "kaplan-meier")
    Sc = stepfun(kmfit.C$time,  pmax(c(1, kmfit.C$surv), trim.C) )  # Bound Sc away from 0 so we do not get extreme weights
  }
  
  
  formula.X = formula(paste("Surv(Q, X, delta.1) ~ ", 
                            paste(covariates.X, collapse = "+"), collapse = ""))
  formula.Q2 = formula(paste("Surv(tau-X, tau-Q, delta.1) ~ ", 
                             paste(covariates.Q, collapse = "+"), collapse = ""))
  Z.X = as.matrix(dat[,covariates.X])
  Z.Q = as.matrix(dat[,covariates.Q])
  
  # remove the subjects with T-Q < 10^{-7} when fitting Cox models
  # such subjects will trigger error when fitting a Cox model
  id.cox = ((dat$X - dat$Q) >= 10^{-7})
  dat.cox = dat[id.cox, ]
  
  # fit Cox-PH model for T|Z
  fit.X = coxph(formula.X, data = dat.cox)
  # fit Cox-PH model for (tau-Q)|Z
  fit.Q2 = coxph(formula.Q2, data = dat.cox)
  basehaz.X = basehaz(fit.X, centered = FALSE)
  basehaz.Q2 = basehaz(fit.Q2, centered = FALSE)
  beta.X = coef(fit.X)
  beta.Q2 = coef(fit.Q2)
  X = dat$X
  Q = dat$Q
  Q2 = tau - Q
  delta = dat$delta
  
  m1qz = diag(integral_F(Q, nu1, delta, Sc, Z.X, basehaz.X, beta.X))
  
  Gtz = as.vector((baseS(tau-X, basehaz.Q2))^exp(Z.Q %*% beta.Q2))
  Gqz = as.vector((baseS(tau-Q, basehaz.Q2))^exp(Z.Q %*% beta.Q2))
  Fqz = as.vector(1-(baseS(Q, basehaz.X))^exp(Z.X %*% beta.X))
  
  
  # DR estimator
  # compute the denominator of the estimator
  DDen1 = delta/(Sc(X)*Gtz)
  DDen2 = m1qz/(Gqz*(1-Fqz)) 
  DDen3 = integral_JG_m(Q, X, delta, Sc, Z.Q, Z.X, 
                        basehaz.Q2, basehaz.X, beta.Q2, beta.X, tau, nu1)
  
  Enu1tz = as.vector(integral_F(tau, nu1, delta, Sc, Z.X, basehaz.X, beta.X))
  
  DDen4_1 = delta/Sc(X) + m1qz/(1-Fqz)
  DDen4 = Enu1tz/(1-Fqz)
  
  
  # compute the numerator of the estimator
  NNum1 = matrix(nrow = n, ncol = m)
  NNum2 = matrix(nrow = n, ncol = m)
  NNum3 = matrix(nrow = n, ncol = m)
  NNum4 = matrix(nrow = n, ncol = m)
  NNum5 = matrix(nrow = n, ncol = m)
  
  for(i in 1:m){
    t0 = t0.list[i]
    nu <- function(t,t00=t0){
      # indicator function
      result = as.numeric(t>t00)
      
      return(result)
    }
    
    mqz = diag(integral_F(Q, nu, delta, Sc, Z.X, basehaz.X, beta.X))
    
    NNum1[,i] = delta*nu(X)/(Sc(X)*Gtz)
    NNum2[,i] = mqz/(Gqz*(1-Fqz))
    NNum3[,i] = integral_JG_m(Q, X, delta, Sc, Z.Q, Z.X, 
                              basehaz.Q2, basehaz.X, beta.Q2, beta.X, tau, nu)
    
    Enutz = as.vector(integral_F(tau, nu, delta, Sc, Z.X, basehaz.X, beta.X))
    
    NNum4[,i] = delta*nu(X)/Sc(X) + mqz/(1-Fqz)
    NNum5[,i] = Enutz/(1-Fqz)
  }
  
  den1 = sum(DDen1)
  den2 = sum(DDen2)
  den3 = sum(DDen3)
  den4_1 = sum(DDen4_1)
  den4 = sum(DDen4)
  
  num1 = colSums(NNum1)
  num2 = colSums(NNum2)
  num3 = colSums(NNum3)
  num4 = colSums(NNum4)
  num5 = colSums(NNum5)
  
  # estimate the untruncated rate
  beta2 = 1/(den4/n)  # estimator based on Q model
  beta = 1/(den1/n)   # truncation rate estimator based on T model
  
  # The estimator based on EIF
  est_DR = (num1+num2-num3)/(den1+den2-den3)
  
  IF_DR = beta*((NNum1+NNum2-NNum3) - est_DR*matrix(rep(DDen1+DDen2-DDen3, m), nrow = n, byrow = FALSE))
  se_DR = sqrt(colMeans(IF_DR^2))/sqrt(n)
  
  # The IPQW estimator
  est_IPQW = num1/den1
  SF_IPQW = NNum1 - est_IPQW*matrix(rep(DDen1, m), nrow = n, byrow = FALSE)
  se_IPQW = sqrt(colMeans(SF_IPQW^2)/(den1/n)^2)/sqrt(n)
  
  # The IPSW estimator
  est_IPSW = num4/den4_1
  est_IPSW2 = num5/den4
  SF_IPSW = NNum4 - est_IPSW*matrix(rep(DDen4, m), nrow = n, byrow = FALSE)
  se_IPSW = sqrt(colMeans(SF_IPSW^2)/(den4/n)^2)/sqrt(n)
  SF_IPSW2 = NNum5 - est_IPSW2*matrix(rep(DDen4, m), nrow = n, byrow = FALSE)
  se_IPSW2 = sqrt(colMeans(SF_IPSW2^2)/(den4/n)^2)/sqrt(n)
  
  
  return(list(est_DR = est_DR, 
              est_IPQW = est_IPQW,
              est_IPSW = est_IPSW, est_IPSW2 = est_IPSW2,
              beta.X = beta.X, beta.Q2 = beta.Q2,
              se_DR = se_DR,
              se_IPQW = se_IPQW,
              se_IPSW = se_IPSW,
              se_IPSW2 = se_IPSW2,
              beta = beta,
              beta2 = beta2))
}




# functions for boostrap 
est_DR.cox.naive.boot <- function(dat, id, nu, covariates.X, covariates.Q, trim.C){
    dat = dat[id,]
    est.cox = c1.est_DR.cox(dat, nu, covariates.X, covariates.Q, trim.C)
    
    est.cox.vector = c(est_DR = est.cox$est_DR, est_IPQW = est.cox$est_IPQW, 
                       est_IPSW = est.cox$est_IPSW, est_IPSW2 = est.cox$est_IPSW2,
                       est_naive = mean(nu(dat$AgeDeath)),
                       se_DR = est.cox$se_DR, se_IPQW = est.cox$se_IPQW, 
                       se_IPSW = est.cox$se_IPSW, se_IPSW2 = est.cox$se_IPSW2,
                       se_naive = sd(nu(dat$AgeDeath)/sqrt(nrow(dat))), 
                       beta = est.cox$beta, beta2 = est.cox$beta2)
    
    return(est.cox.vector)
}


c1.estSurv_DR.cox.boot <- function(dat, id, t0, covariates.X, covariates.Q, trim.C){
  dat = dat[id,]
  
  estSurv.cox = c1.estSurv_DR.cox(t0, dat, covariates.X, covariates.Q, trim.C)
  
  estSurv.cox.vector = c(est_DR = estSurv.cox$est_DR, 
                         est_IPQW = estSurv.cox$est_IPQW, 
                         est_IPSW = estSurv.cox$est_IPSW, 
                         est_IPSW2 = estSurv.cox$est_IPSW2,
                         beta = estSurv.cox$beta, beta2 = estSurv.cox$beta2)
}

c1.splitboot_estSurv_DR.cox <- function(seed, n.boot, dat, t0, covariates.X, covariates.Q, trim.C){
  set.seed(seed)
  bootresult.cox = boot(dat, c1.estSurv_DR.cox.boot, R = n.boot,
                          t0 = t0, covariates.X = covariates.X,
                          covariates.Q = covariates.Q, trim.C = trim.C)
  
  return(bootresult.cox)  
}







# cf-RF-RF -------------------------------------------------------------------------

## function for computing the cf estimator using LTRCforests::ltrcrrf() to fit F and G
## The * notation in the code denote how time consuming it is to run that line. 
## More *'s means more time consuming.

# K-fold cross fiting is applyed
# K: number of folds in cross fitting
c1.est_cf.ltrcrrf <- function(dat, nu, covariates.X, covariates.Q, K, 
                              mtry = Inf, ntree = 100, trim = 1e-7, trim.C = 1e-7){
  
  tau = max(dat$X)+1   # used to reverse the time scale
  tau1 = min(dat$X)   # tau1 and tau2 are used to bound the estimates of F and G
  tau2 = max(dat$Q)
  
  dat$delta.1 = rep(1, nrow(dat))   # indicator used for fitting X and Q distributions
  dat$Q2 = tau - dat$Q    
  dat$X2 = tau - dat$X
  
  formula.X = formula(paste("Surv(Q, X, delta.1) ~ ", 
                            paste(covariates.X, collapse = "+"), collapse = ""))
  formula.Q2 = formula(paste("Surv(X2, Q2, delta.1) ~ ", 
                             paste(covariates.Q, collapse = "+"), collapse = ""))
  
  n = nrow(dat)
  
  # Define a constant 1 function nu1
  nu1 <- function(x){
    return(rep(1, length(x)))
  }
  
  
  # split the data into K folds 
  folds <- cvFolds(n, K)
  DDen1 = rep(NA, n)
  DDen2 = rep(NA, n)
  DDen3 = rep(NA, n)
  DDen4_1 = rep(NA, n)
  DDen4 = rep(NA, n)
  NNum1 = rep(NA, n)
  NNum2 = rep(NA, n)
  NNum3 = rep(NA, n)
  NNum4 = rep(NA, n)
  NNum5 = rep(NA, n)
  
  ct = 0
  for(k in 1:K){
    dat.est = dat[folds$subsets[folds$which == k], ]
    dat.fit = dat[folds$subsets[folds$which != k], ]  
    dat.fit = dat.fit[dat.fit$X - dat.fit$Q > 1e-7, ]  # remove the observation that result in "an interval has effective length 0"
    
    nk = nrow(dat.est)
    
    newdat = dat.est
    newdat$Q = 0
    newdat$X2 = 0
    newdat$X = tau
    newdat$Q2 = tau
    newdat$delta.1 = 1
    
    # Use dat.fit to estimate Sc from PL estimator
    # If no censoring, Sc is the function that is constant 1
    if(sum(1-dat.fit$delta) == 0){
      Sc = function(x){return(1)}
    }else{
      kmfit.C = survfit(Surv(Q, X, 1-delta)~1, data = dat.fit, type = "kaplan-meier")
      Sc = stepfun(kmfit.C$time,  pmax(c(1, kmfit.C$surv), trim.C) )  # Bound Sc away from 0 so we do not get extreme weights
    }
    
    # use dat.fit to estimate F and G
    X.obj <- ltrcrrf(formula = formula.X, data = dat.fit,
                     mtry = mtry, ntree = ntree)
    Q2.obj <- ltrcrrf(formula = formula.Q2, data = dat.fit,
                      mtry = mtry, ntree = ntree)
    
    jumps.X = sort(dat.fit$X)
    X.pred.dF = predictProb(object = X.obj, newdata = newdat, 
                            time.eval = jumps.X)   #**
    # row in $survival.probs corresponds to a time point in jumps.T, column corresponds to a subject
    
    # Bound the estimates such that \hat{F}(\tau_2|Z_i) <= 1-trim for all i
    id.X = which(jumps.X <= tau2)
    X.pred.dF$survival.probs[id.X,] <- pmin(X.pred.dF$survival.probs[id.X,], 1-trim)
    
    
    jumps.Q = sort(dat.fit$Q)
    odr = order(tau-jumps.Q)
    Q.pred.dG = predictProb(object = Q2.obj, newdata = newdat, 
                            time.eval = (tau-jumps.Q)[odr])   #**
    # Bound the estimates such that \hat{G}(\tau_1|Z_i) >= trim for all i
    id.Q = which(jumps.Q >= tau1)
    Q.pred.dG$survival.probs[id.Q, ] <- pmax(Q.pred.dG$survival.probs[id.Q, ], trim)
    
    Gtz = diag(predict.survProb.mx(Q.pred.dG, dat.est$X2))
    Gqz = diag(predict.survProb.mx(Q.pred.dG, dat.est$Q2))
    Fqz = 1-diag(predict.survProb.mx(X.pred.dF, dat.est$Q))
    mqz = int_nu_dF.LTRCforests(X.pred.dF, dat.est$Q, nu, newdat$delta, Sc)
    
    m1qz = int_nu_dF.LTRCforests(X.pred.dF, dat.est$Q, nu1, newdat$delta, Sc)
    
    Q.est = dat.est$Q
    X.est = dat.est$X
    delta.est = dat.est$delta
    n.est = nrow(dat.est)
    atrisk = matrix(nrow = n.est, ncol = length(jumps.Q))
    for(i in 1:n.est){
      atrisk[i,] = (Q.est[i] <= jumps.Q & jumps.Q < X.est[i])    # at risk indicator for subject i
    }
    
    Fvz.mx = 1-predict.survProb.mx(X.pred.dF, jumps.Q)
    Gvz.mx = t(Q.pred.dG$survival.probs)[,odr]
    
    # compute the denominator of the estimator
    DDen1[(ct+1):(ct+nk)] = delta.est/(Sc(X.est)*pmax(Gtz,trim))
    DDen2[(ct+1):(ct+nk)] = m1qz/pmax(Gqz*(1-Fqz), trim)
    DDen3[(ct+1):(ct+nk)] = int_JG_m.LTRCforests(X.pred.dF, jumps.Q, Fvz.mx, Gvz.mx, atrisk, 
                                                 tau, nu1, trim, delta.est, Sc)  #*
    
    # compute the numerator of the estimator
    NNum1[(ct+1):(ct+nk)] = delta.est*nu(X.est)/(Sc(X.est)*pmax(Gtz,trim))
    NNum2[(ct+1):(ct+nk)] = mqz/(pmax(Gqz, trim)*pmax(1-Fqz, trim))
    NNum3[(ct+1):(ct+nk)] = int_JG_m.LTRCforests(X.pred.dF, jumps.Q, Fvz.mx, Gvz.mx, atrisk, 
                                                 tau, nu, trim, delta.est, Sc)  #*
    
    # denominator and numerator for the IPSW estimator
    Enutz = as.vector(int_nu_dF.mx.LTRCforests(X.pred.dF, tau, nu, delta.est, Sc))
    Enu1tz = as.vector(int_nu_dF.mx.LTRCforests(X.pred.dF, tau, nu1, delta.est, Sc))
    
    DDen4_1[(ct+1):(ct+nk)] = delta.est/Sc(X.est) + m1qz/pmax(1-Fqz, trim)
    DDen4[(ct+1):(ct+nk)] = Enu1tz/pmax(1-Fqz, trim)
    
    NNum4[(ct+1):(ct+nk)] = delta.est*nu(X.est)/Sc(X.est) + mqz/pmax(1-Fqz, trim)
    NNum5[(ct+1):(ct+nk)] = Enutz/pmax(1-Fqz, trim)
    
    ct = ct + nk
  }
  
  den1 = sum(DDen1)
  den2 = sum(DDen2)
  den3 = sum(DDen3)
  den4_1 = sum(DDen4_1)
  den4 = sum(DDen4)
  
  num1 = sum(NNum1)
  num2 = sum(NNum2)
  num3 = sum(NNum3)
  num4 = sum(NNum4)
  num5 = sum(NNum5)
  
  # The DR estimator
  est_DR = (num1+num2-num3)/(den1+den2-den3)
  beta = 1/(den1/n)
  IF_DR = beta*((NNum1+NNum2-NNum3) - est_DR*(DDen1+DDen2-DDen3))
  se_DR = sqrt(mean(IF_DR^2))/sqrt(n)
  
  # The IPQW estimator
  est_IPQW = num1/den1
  SF_IPQW = NNum1 - est_IPQW*DDen1
  se_IPQW = sqrt(mean(SF_IPQW^2)/(den1/n)^2)/sqrt(n)
  
  # The IPSW estimator
  est_IPSW = num4/den4_1
  est_IPSW2 = num5/den4
  SF_IPSW = NNum4 - est_IPSW*DDen4
  se_IPSW = sqrt(mean(SF_IPSW^2)/(den4/n)^2)/sqrt(n)
  SF_IPSW2 = NNum5 - est_IPSW2*DDen4
  se_IPSW2 = sqrt(mean(SF_IPSW2^2)/(den4/n)^2)/sqrt(n)
  
  beta2 = 1/(den4/n)
  
  return(list(est_DR = est_DR,
              est_IPQW = est_IPQW,
              est_IPSW = est_IPSW,
              est_IPSW2 = est_IPSW2,
              se_DR = se_DR,
              se_IPQW = se_IPQW,
              se_IPSW = se_IPSW,
              se_IPSW2 = se_IPSW2,
              beta = beta,
              beta2 = beta2)
  )
  
}


### function that take input t0 (could be a vector) and estimate the survival probability at t0
c1.estSurv_cf.ltrcrrf <- function(t0, dat, covariates.X, covariates.Q, K, 
                              mtry = Inf, ntree = 100, trim = 1e-7, trim.C = 1e-7){
  t0.list = t0
  n = nrow(dat)
  m = length(t0.list)
  tau = max(dat$X)+1   # used to reverse the time scale
  tau1 = min(dat$X)   # tau1 and tau2 are used to bound the estimates of F and G
  tau2 = max(dat$Q)
  
  dat$delta.1 = rep(1, nrow(dat))   # indicator used for fitting X and Q distributions
  dat$Q2 = tau - dat$Q    
  dat$X2 = tau - dat$X
  
  formula.X = formula(paste("Surv(Q, X, delta.1) ~ ", 
                            paste(covariates.X, collapse = "+"), collapse = ""))
  formula.Q2 = formula(paste("Surv(X2, Q2, delta.1) ~ ", 
                             paste(covariates.Q, collapse = "+"), collapse = ""))
  
  # Define a constant 1 function nu1
  nu1 <- function(x){
    return(rep(1, length(x)))
  }
  
  # split the data into K folds 
  folds <- cvFolds(n, K)
  DDen1 = rep(NA, n)
  DDen2 = rep(NA, n)
  DDen3 = rep(NA, n)
  DDen4_1 = rep(NA, n)
  DDen4 = rep(NA, n)
  NNum1 = matrix(nrow = n, ncol = m)
  NNum2 = matrix(nrow = n, ncol = m)
  NNum3 = matrix(nrow = n, ncol = m)
  NNum4 = matrix(nrow = n, ncol = m)
  NNum5 = matrix(nrow = n, ncol = m)
  
  ct = 0
  for(k in 1:K){
    dat.est = dat[folds$subsets[folds$which == k], ]
    dat.fit = dat[folds$subsets[folds$which != k], ]  
    dat.fit = dat.fit[dat.fit$X - dat.fit$Q > 1e-7, ]  # remove the observation that result in "an interval has effective length 0"
    
    nk = nrow(dat.est)
    
    newdat = dat.est
    newdat$Q = 0
    newdat$X2 = 0
    newdat$X = tau
    newdat$Q2 = tau
    newdat$delta.1 = 1
    
    # Use dat.fit to estimate Sc from PL estimator
    # If no censoring, Sc is the function that is constant 1
    if(sum(1-dat.fit$delta) == 0){
      Sc = function(x){return(1)}
    }else{
      kmfit.C = survfit(Surv(Q, X, 1-delta)~1, data = dat.fit, type = "kaplan-meier")
      Sc = stepfun(kmfit.C$time,  pmax(c(1, kmfit.C$surv), trim.C) )  # Bound Sc away from 0 so we do not get extreme weights
    }
    
    # use dat.fit to estimate F and G
    X.obj <- ltrcrrf(formula = formula.X, data = dat.fit,
                     mtry = mtry, ntree = ntree)
    Q2.obj <- ltrcrrf(formula = formula.Q2, data = dat.fit,
                      mtry = mtry, ntree = ntree)
    
    jumps.X = sort(dat.fit$X)
    X.pred.dF = predictProb(object = X.obj, newdata = newdat, 
                            time.eval = jumps.X)   #**
    # row in $survival.probs corresponds to a time point in jumps.T, column corresponds to a subject
    
    # bound the estimates such that \hat{F}(\tau_2|Z_i) <= 1-trim for all i
    id.X = which(jumps.X <= tau2)
    X.pred.dF$survival.probs[id.X,] <- pmin(X.pred.dF$survival.probs[id.X,], 1-trim)
    
    
    jumps.Q = sort(dat.fit$Q)
    odr = order(tau-jumps.Q)
    Q.pred.dG = predictProb(object = Q2.obj, newdata = newdat, 
                            time.eval = (tau-jumps.Q)[odr])   #**
    # bound the estimates such that \hat{G}(\tau_1|Z_i) >= trim for all i
    id.Q = which(jumps.Q >= tau1)
    Q.pred.dG$survival.probs[id.Q, ] <- pmax(Q.pred.dG$survival.probs[id.Q, ], trim)
    
    Gtz = diag(predict.survProb.mx(Q.pred.dG, dat.est$X2))
    Gqz = diag(predict.survProb.mx(Q.pred.dG, dat.est$Q2))
    Fqz = 1-diag(predict.survProb.mx(X.pred.dF, dat.est$Q))
    
    m1qz = int_nu_dF.LTRCforests(X.pred.dF, dat.est$Q, nu1, newdat$delta, Sc)
    
    Q.est = dat.est$Q
    X.est = dat.est$X
    delta.est = dat.est$delta
    n.est = nrow(dat.est)
    atrisk = matrix(nrow = n.est, ncol = length(jumps.Q))
    for(i in 1:n.est){
      atrisk[i,] = (Q.est[i] <= jumps.Q & jumps.Q < X.est[i])    # at risk indicator for subject i
    }
    
    Fvz.mx = 1-predict.survProb.mx(X.pred.dF, jumps.Q)
    Gvz.mx = t(Q.pred.dG$survival.probs)[,odr]
    
    # compute the denominator of the estimator
    DDen1[(ct+1):(ct+nk)] = delta.est/(Sc(X.est)*pmax(Gtz,trim))
    DDen2[(ct+1):(ct+nk)] = m1qz/pmax(Gqz*(1-Fqz), trim)
    DDen3[(ct+1):(ct+nk)] = int_JG_m.LTRCforests(X.pred.dF, jumps.Q, Fvz.mx, Gvz.mx, atrisk, 
                                                 tau, nu1, trim, delta.est, Sc)  #*
    
    # denominator and numerator for the IPSW estimator
    Enu1tz = as.vector(int_nu_dF.mx.LTRCforests(X.pred.dF, tau, nu1, delta.est, Sc))
    
    DDen4_1[(ct+1):(ct+nk)] = delta.est/Sc(X.est) + m1qz/pmax(1-Fqz, trim)
    DDen4[(ct+1):(ct+nk)] = Enu1tz/pmax(1-Fqz, trim)
    
    
    for(i in 1:m){
      t0 = t0.list[i]
      nu <- function(t,t00=t0){
        # indicator function
        result = as.numeric(t>t00)
        
        return(result)
      }
      
      mqz = int_nu_dF.LTRCforests(X.pred.dF, dat.est$Q, nu, newdat$delta, Sc)
      
      # compute the numerators of the estimators
      NNum1[(ct+1):(ct+nk), i] = delta.est*nu(X.est)/(Sc(X.est)*pmax(Gtz,trim))
      NNum2[(ct+1):(ct+nk), i] = mqz/(pmax(Gqz, trim)*pmax(1-Fqz, trim))
      NNum3[(ct+1):(ct+nk), i] = int_JG_m.LTRCforests(X.pred.dF, jumps.Q, Fvz.mx, Gvz.mx, atrisk, 
                                                      tau, nu, trim, delta.est, Sc)  #*
      
      Enutz = as.vector(int_nu_dF.mx.LTRCforests(X.pred.dF, tau, nu, delta.est, Sc))
      
      NNum4[(ct+1):(ct+nk), i] = delta.est*nu(X.est)/Sc(X.est) + mqz/pmax(1-Fqz, trim)
      NNum5[(ct+1):(ct+nk), i] = Enutz/pmax(1-Fqz, trim)
    }
    
    ct = ct + nk
  }
  
  den1 = sum(DDen1)
  den2 = sum(DDen2)
  den3 = sum(DDen3)
  den4 = sum(DDen4)
  den4_1 = sum(DDen4_1)
  
  num1 = colSums(NNum1)
  num2 = colSums(NNum2)
  num3 = colSums(NNum3)
  num4 = colSums(NNum4)
  num5 = colSums(NNum5)
  
  # estimators for the untruncated rate
  beta = 1/(den1/n)
  beta2 = 1/(den4/n)
  
  # The DR estimator
  est_DR = (num1+num2-num3)/(den1+den2-den3)
  IF_DR = beta*((NNum1+NNum2-NNum3) - est_DR*matrix(rep(DDen1+DDen2-DDen3, m), nrow = n, byrow = FALSE))
  se_DR = sqrt(colMeans(IF_DR^2))/sqrt(n)
  
  # The IPQW estimator
  est_IPQW = num1/den1
  SF_IPQW = NNum1 - est_IPQW*matrix(rep(DDen1, m), nrow = n, byrow = FALSE)
  se_IPQW = sqrt(colMeans(SF_IPQW^2)/(den1/n)^2)/sqrt(n)
  
  # The IPSW estimator
  est_IPSW = num4/den4_1
  est_IPSW2 = num5/den4
  SF_IPSW = NNum4 - est_IPSW*matrix(rep(DDen4, m), nrow = n, byrow = FALSE)
  se_IPSW = sqrt(colMeans(SF_IPSW^2)/(den4/n)^2)/sqrt(n)
  SF_IPSW2 = NNum5 - est_IPSW2*matrix(rep(DDen4, m), nrow = n, byrow = FALSE)
  se_IPSW2 = sqrt(colMeans(SF_IPSW2^2)/(den4/n)^2)/sqrt(n)
  
  
  return(list(est_DR = est_DR,
              est_IPQW = est_IPQW,
              est_IPSW = est_IPSW,
              est_IPSW2 = est_IPSW2,
              se_DR = se_DR,
              se_IPQW = se_IPQW,
              se_IPSW = se_IPSW,
              se_IPSW2 = se_IPSW2,
              beta = beta,
              beta2 = beta2))
}



# function for bootstrap
c1.estSurv_cf.ltrcrrf.boot <- function(dat, id, t0, covariates.X, covariates.Q, K, 
                                       mtry, ntree, trim, trim.C){
  dat = dat[id,]
  
  result = c1.estSurv_cf.ltrcrrf(t0, dat, covariates.X, covariates.Q, K,
                                 mtry, ntree, trim, trim.C)
  
  est.ltrcrrf = c(est_DR = result$est_DR,
                  est_IPQW = result$est_IPQW,
                  est_IPSW = result$est_IPSW,
                  est_IPSW2 = result$est_IPSW2, 
                  beta = result$beta,
                  beta2 = result$beta2)
  
  return(est.ltrcrrf)
  
}

c1.splitboot_estSurv.ltrcrrf <- function(seed, n.boot, dat, t0, covariates.X, covariates.Q, 
                                         K, mtry, ntree, trim, trim.C){
  set.seed(seed)
  bootresult.ltrcrrf = boot(dat, c1.estSurv_cf.ltrcrrf.boot, R = n.boot,
                        t0 = t0, covariates.X = covariates.X, covariates.Q = covariates.Q,
                        K = K, mtry = mtry, ntree = ntree, trim = trim, trim.C = trim.C)
  
  return(bootresult.ltrcrrf)  
}









# naive ----------------------------------------------------------------------------
c1.estSurv_naive <- function(t0, dat){
  n = nrow(dat)
  # fit Sc by product-limit estimator
  if(sum(1-dat$delta) == 0){
    Sc = function(x){return(1)}
  }else{
    kmfit.C = survfit(Surv(Q, X, 1-delta)~1, data = dat, type = "kaplan-meier")
    Sc = stepfun(kmfit.C$time,  pmax(c(1, kmfit.C$surv), trim.C) ) # Bound Sc away from 0 so we do not get extreme weights
  }
  
  nu <- function(t,t00=t0){
    # indicator function
    result = as.numeric(t>t00)
    
    return(result)
  }
  
  est.naive = sum(dat$delta*nu(dat$X)/Sc(dat$X))/sum(dat$delta/Sc(dat$X))
  names(est.naive) = "est.naive"
  # se.naive = NA
  # names(se.naive) = "se.naive"
  
  return(c(est.naive = est.naive))

}

# function for bootstrap
estSurv_KM.boot <- function(dat, id, t0){
  dat = dat[id,]

  kmfit.T = survfit(Surv(X, delta)~1, data = dat)
  surv.km = stepfun(kmfit.T$time, c(1, kmfit.T$surv))
  
  return(surv.km(t0))
}

estSurv_pl.boot <- function(dat, id, t0){
  dat = dat[id,]
  
  plfit.T = survfit(Surv(Q, X, delta)~1, data = dat)
  surv.pl = stepfun(plfit.T$time, c(1, plfit.T$surv))
  
  return(surv.pl(t0))
}





