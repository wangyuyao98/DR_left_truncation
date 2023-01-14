# function for computing the dr and cf estimator
# with IPCW weights to handle c2 type non-informative censoring: Q<C a.s.


# dr-Cox-Cox ------------------------------------------------------------------
# dr, IPW.Q, Reg.T1, Reg.T2 estimators that use Cox models for fitting T~Z and (tau-Q)~Z models 
# The Cox models are fitted on the entire sample
# The argument `dat' should contains X := min(T,C), Q, delta := I(T<C), and the covariates
c2.est_DR.cox <- function(dat, nu, covariates.T, covariates.Q, trim.C = 1e-7){
  
  n = nrow(dat)
  tau = max(dat$X) + 1
  delta = dat$delta
  X = dat$X
  Q = dat$Q
  Q2 = tau-Q
  
  dat$delta.1 = rep(1, nrow(dat))
  
  # compute the survival curve Sd for the residual censoring time D:=C-Q from KM estimator
  # If no censoring, Sd is the function that is constant 1
  if(sum(1-dat$delta) == 0){
    Sd = function(x){return(1)}
  }else{
    kmfit.D = survfit(Surv(X-Q, 1-delta)~1, data = dat, type = "kaplan-meier")
    Sd = stepfun(kmfit.D$time,  pmax(c(1, kmfit.D$surv), trim.C) )  # Bound Sc away from 0 so we do not get extreme weights
  }
  # IPCW weights
  wd = rep(0, n)
  id1 = which(delta == 1)
  wd[id1] = 1/pmax(Sd(X-Q)[id1], trim.C)
  dat$wd = wd
  
  formula.T = formula(paste("Surv(Q, X, delta.1) ~ ", 
                            paste(covariates.T, collapse = "+"), collapse = ""))
  formula.Q2 = formula(paste("Surv(tau-X, tau-Q, delta.1) ~ ", 
                             paste(covariates.Q, collapse = "+"), collapse = ""))
  Z.T = as.matrix(dat[,covariates.T])
  Z.Q = as.matrix(dat[,covariates.Q])
  
  # remove the subjects with T-Q < 10^{-7} when fitting Cox models - such subjects will trigger error when fitting a Cox model
  # IPCW weighted uncensored subjects, i.e., those with delta = 1
  id.cox = (dat$delta == 1) & ((dat$X - dat$Q) >= 10^{-7}) 
  dat.cox = dat[id.cox, ]
 
  # fit Cox-PH model for T|Z
  fit.T = coxph(formula.T, data = dat.cox, weights = dat.cox$wd)
  # fit Cox-PH model for (tau-Q)|Z
  fit.Q2 = coxph(formula.Q2, data = dat.cox, weights = dat.cox$wd)
  basehaz.T = basehaz(fit.T, centered = FALSE)
  basehaz.Q2 = basehaz(fit.Q2, centered = FALSE)
  beta.T = coef(fit.T)
  beta.Q2 = coef(fit.Q2)
  
  Gtz = (baseS(tau-X, basehaz.Q2))^exp(Z.Q %*% beta.Q2)
  Gqz = (baseS(tau-Q, basehaz.Q2))^exp(Z.Q %*% beta.Q2)
  Fqz = 1-(baseS(Q, basehaz.T))^exp(Z.T %*% beta.T)
  mqz = diag(integral_F(Q, nu, Z.T, basehaz.T, beta.T))
  
  # compute the denominator of the estimator
  DDen1 = wd*1/Gtz
  DDen2 = wd*Fqz/(Gqz*(1-Fqz))
  DDen3 = wd*integral_JG_F(Q, X, Z.Q, Z.T, basehaz.Q2, basehaz.T, beta.Q2, beta.T, tau)
  
  den1 = sum(DDen1)
  den2 = sum(DDen2)
  den3 = sum(DDen3)
  
  # compute the numerator of the estimator
  NNum1 = wd*nu(X)/Gtz
  NNum2 = wd*mqz/(Gqz*(1-Fqz))
  NNum3 = wd*integral_JG_m(Q, X, Z.Q, Z.T, basehaz.Q2, basehaz.T, beta.Q2, beta.T, tau, nu)
  
  num1 = sum(NNum1)
  num2 = sum(NNum2)
  num3 = sum(NNum3)
  
  # The estimator based on EIF
  est_DR = (num1+num2-num3)/(den1+den2-den3)
  beta = 1/(den1/n)
  IF_DR = beta*((NNum1+NNum2-NNum3) - est_DR*(DDen1+DDen2-DDen3))
  se_DR = sqrt(mean(IF_DR^2))/sqrt(n)
  
  # The IPQW estimator
  est_IPQW = num1/den1
  SF_IPQW = NNum1 - est_IPQW*DDen1
  se_IPQW = sqrt(mean(SF_IPQW^2)/(den1/n)^2)/sqrt(n)
  
  # The IPSW estimator
  Enutz = as.vector(integral_F(tau, nu, Z.T, basehaz.T, beta.T))
  
  DDen4 = wd*1/(1-Fqz)
  NNum4 = wd*(nu(X) + mqz/(1-Fqz))
  NNum5 = wd*Enutz/(1-Fqz)
  
  den4 = sum(DDen4)
  num4 = sum(NNum4)
  num5 = sum(NNum5)
  
  est_IPSW = num4/den4
  est_IPSW2 = num5/den4
  SF_IPSW = NNum4 - est_IPSW*DDen4
  se_IPSW = sqrt(mean(SF_IPSW^2)/(den4/n)^2)/sqrt(n)
  SF_IPSW2 = NNum5 - est_IPSW2*DDen4
  se_IPSW2 = sqrt(mean(SF_IPSW2^2)/(den4/n)^2)/sqrt(n)
  
  beta2 = 1/(den4/n)
  
  return(list(est_DR = est_DR, 
              est_IPQW = est_IPQW,
              est_IPSW = est_IPSW, est_IPSW2 = est_IPSW2,
              beta.T = beta.T, beta.Q2 = beta.Q2,
              se_DR = se_DR,
              se_IPQW = se_IPQW,
              se_IPSW = se_IPSW,
              se_IPSW2 = se_IPSW2,
              beta = beta,
              beta2 = beta2))
  
}


c.est_DR.cox <- function(dat, nu, covariates.T, covariates.Q, trim.C, type.C){
  if(type.C == "c1"){
    est.cox = c1.est_DR.cox(dat, nu, covariates.T, covariates.Q, trim.C)
  }else if(type.C == "c2"){
    est.cox = c2.est_DR.cox(dat, nu, covariates.T, covariates.Q, trim.C)
  }else{
    stop("The censoring type should be either c1 or c2!")
  }
  
  return(est.cox)
}



## function that take input t0 (posibaly a vector) and estimate the survival probability(s) at t0
## with dr-Cox-Cox, IPW.Q-Cox, Reg.T1, Reg.T2 estimators

# c2 type censoring 
c2.estSurv_DR.cox <- function(t0, dat, covariates.T, covariates.Q, trim.C){
  t0.list = t0
  m = length(t0.list)
  n = nrow(dat)
  
  tau = max(dat$X) + 1
  dat$delta.1 = rep(1, nrow(dat))
  
  delta = dat$delta
  X = dat$X
  Q = dat$Q
  Q2 = tau-Q
  
  # compute the survival curve Sd for the residual censoring time D:=C-Q from KM estimator
  # If no censoring, Sd is the function that is constant 1
  if(sum(1-dat$delta) == 0){
    Sd = function(x){return(1)}
  }else{
    kmfit.D = survfit(Surv(X-Q, 1-delta)~1, data = dat, type = "kaplan-meier")
    Sd = stepfun(kmfit.D$time,  pmax(c(1, kmfit.D$surv), trim.C) )  # Bound Sc away from 0 so we do not get extreme weights
  }
  # IPCW weights
  wd = rep(0, nrow(dat))
  id1 = which(dat$delta == 1)
  wd[id1] = 1/Sd(dat$X-dat$Q)[id1]
  dat$wd = wd
  
  formula.T = formula(paste("Surv(Q, X, delta.1) ~ ", 
                            paste(covariates.T, collapse = "+"), collapse = ""))
  formula.Q2 = formula(paste("Surv(tau-X, tau-Q, delta.1) ~ ", 
                             paste(covariates.Q, collapse = "+"), collapse = ""))
  Z.T = as.matrix(dat[,covariates.T])
  Z.Q = as.matrix(dat[,covariates.Q])
  
  # remove the subjects with T-Q < 10^{-7} when fitting Cox models - such subjects will trigger error when fitting a Cox model
  # IPCW weighted uncensored subjects, i.e., those with delta = 1
  id.cox = (dat$delta == 1) & ((dat$X - dat$Q) >= 10^{-7}) 
  dat.cox = dat[id.cox, ]
  
  # fit Cox-PH model for T|Z
  fit.T = coxph(formula.T, data = dat.cox, weights = dat.cox$wd)
  # fit Cox-PH model for (tau-Q)|Z
  fit.Q2 = coxph(formula.Q2, data = dat.cox, weights = dat.cox$wd)
  basehaz.T = basehaz(fit.T, centered = FALSE)
  basehaz.Q2 = basehaz(fit.Q2, centered = FALSE)
  beta.T = coef(fit.T)
  beta.Q2 = coef(fit.Q2)
  
  Gtz = (baseS(tau-X, basehaz.Q2))^exp(Z.Q %*% beta.Q2)
  Gqz = (baseS(tau-Q, basehaz.Q2))^exp(Z.Q %*% beta.Q2)
  Fqz = 1-(baseS(Q, basehaz.T))^exp(Z.T %*% beta.T)
  
  # compute the denominator of the estimator
  DDen1 = wd*1/Gtz
  DDen2 = wd*Fqz/(Gqz*(1-Fqz))
  DDen3 = wd*integral_JG_F(Q, X, Z.Q, Z.T, basehaz.Q2, basehaz.T, beta.Q2, beta.T, tau)
  DDen4 = wd*1/(1-Fqz)
  
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
    
    mqz = diag(integral_F(Q, nu, Z.T, basehaz.T, beta.T))
    
    NNum1[,i] = wd*nu(X)/Gtz
    NNum2[,i] = wd*mqz/(Gqz*(1-Fqz))
    NNum3[,i] = wd*integral_JG_m(Q, X, Z.Q, Z.T, basehaz.Q2, basehaz.T, beta.Q2, beta.T, tau, nu)
    
    Enutz = as.vector(integral_F(tau, nu, Z.T, basehaz.T, beta.T))
    
    NNum4[,i] = wd*(nu(X) + mqz/(1-Fqz))
    NNum5[,i] = wd*Enutz/(1-Fqz)
    
  }
  
  den1 = sum(DDen1)
  den2 = sum(DDen2)
  den3 = sum(DDen3)
  den4 = sum(DDen4)
  
  num1 = colSums(NNum1)
  num2 = colSums(NNum2)
  num3 = colSums(NNum3)
  num4 = colSums(NNum4)
  num5 = colSums(NNum5)
  
  # estimate the untruncated rate
  beta = 1/(den1/n)   # estimator based on Q model
  beta2 = 1/(den4/n)  # truncation rate estimator based on T model
  
  # The estimator based on EIF
  est_DR = (num1+num2-num3)/(den1+den2-den3)
  
  IF_DR = beta*((NNum1+NNum2-NNum3) - est_DR*matrix(rep(DDen1+DDen2-DDen3, m), nrow = n, byrow = FALSE))
  se_DR = sqrt(colMeans(IF_DR^2))/sqrt(n)
  
  # The IPQW estimator
  est_IPQW = num1/den1
  SF_IPQW = NNum1 - est_IPQW*matrix(rep(DDen1, m), nrow = n, byrow = FALSE)
  se_IPQW = sqrt(colMeans(SF_IPQW^2)/(den1/n)^2)/sqrt(n)
  
  # The IPSW estimator
  est_IPSW = num4/den4
  est_IPSW2 = num5/den4
  SF_IPSW = NNum4 - est_IPSW*matrix(rep(DDen4, m), nrow = n, byrow = FALSE)
  se_IPSW = sqrt(colMeans(SF_IPSW^2)/(den4/n)^2)/sqrt(n)
  SF_IPSW2 = NNum5 - est_IPSW2*matrix(rep(DDen4, m), nrow = n, byrow = FALSE)
  se_IPSW2 = sqrt(colMeans(SF_IPSW2^2)/(den4/n)^2)/sqrt(n)
  
  
  
  return(list(est_DR = est_DR, 
              est_IPQW = est_IPQW,
              est_IPSW = est_IPSW, est_IPSW2 = est_IPSW2,
              beta.T = beta.T, beta.Q2 = beta.Q2,
              se_DR = se_DR,
              se_IPQW = se_IPQW,
              se_IPSW = se_IPSW,
              se_IPSW2 = se_IPSW2,
              beta = beta,
              beta2 = beta2))
  
}


# function that can handle both `c1' and `c2' types of non-informative censoring
c.estSurv_DR.cox <- function(t0, dat, covariates.T, covariates.Q, trim.C, type.C){
  if(type.C == "c1"){
    result = c1.estSurv_DR.cox(t0, dat, covariates.T, covariates.Q, trim.C)
  }else if(type.C == "c2"){
    result = c2.estSurv_DR.cox(t0, dat, covariates.T, covariates.Q, trim.C)
  }else{
    stop("The censoring type should be either c1 or c2!")
  }
  
  return(result)
}


c.estSurv_DR.cox.boot <- function(dat, id, t0, covariates.T, covariates.Q, trim.C, type.C){
  dat = dat[id,]
  
  result = c.estSurv_DR.cox(t0, dat, covariates.T, covariates.Q, trim.C, type.C)
  
  est = c(DR = result$est_DR, 
          IPQW = result$est_IPQW, 
          IPSW = result$est_IPSW, 
          IPSW2 = result$est_IPSW2)
  
  return(est)
}

c.estSurv_bootse_DR.cox <- function(t0, dat, covariates.T, covariates.Q, trim.C, type.C, n.boot){
  bootresult.cox = boot(dat, c.estSurv_DR.cox.boot, R = n.boot,
                        t0 = t0, covariates.T = covariates.T, 
                        covariates.Q = covariates.Q, trim.C = trim.C, type.C = type.C)

  return(bootresult.cox)
}
  
  
  
  
  
  

# naive estimator -----------------------------------------------------------------
c2.estSurv_naive <- function(t0, dat){
  n = nrow(dat)
  # fit Sc by product-limit estimator
  if(sum(1-dat$delta) == 0){
    Sd = function(x){return(1)}
  }else{
    kmfit.D = survfit(Surv(X-Q, 1-delta)~1, data = dat, type = "kaplan-meier")
    Sd = stepfun(kmfit.D$time,  c(1, kmfit.D$surv) )  # Bound Sc away from 0 so we do not get extreme weights
  }
  
  nu <- function(t,t00=t0){
    # indicator function
    result = as.numeric(t>t00)
    
    return(result)
  }
  
  wd = rep(0, nrow(dat))
  id1 = which(dat$delta == 1)
  wd[id1] = 1/Sd(dat$X-dat$Q)[id1]
  dat$wd = wd
  
  est.naive = sum(wd*nu(dat$X))/sum(wd)
  names(est.naive) = "est.naive"
  # se.naive = NA
  # names(se.naive) = "se.naive"
  
  return(c(est.naive = est.naive))
  
}



c.estSurv_naive <- function(t0, dat, type.C){
  if(type.C == "c1"){
    result = c1.estSurv_naive(t0, dat)
  }else if(type.C == "c2"){
    result = c2.estSurv_naive(t0, dat)
  }else{
    stop("The censoring type should be either c1 or c2!")
  }
  
  return(result)
}


