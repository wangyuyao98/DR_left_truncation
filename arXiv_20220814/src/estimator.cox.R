# estimators that use Cox models for fitting 

## Does not use cross-fitting
est_DR.cox<- function(dat, tau, nu){
    n = nrow(dat)
    
    # remove the subjects with T-Q < 10^{-7} when fitting Cox models
    # such subjects will trigger error when fitting a Cox model
    id.cox = ((dat$time - dat$Q) >= 10^{-7})
    dat.cox = dat[id.cox, ]
    
    # fit Cox-PH model for T|Z
    fit.T = coxph(Surv(Q, time, event = rep(1,nrow(dat.cox))) ~ Z1+Z2, data = dat.cox)
    # fit Cox-PH model for (tau-Q)|Z
    fit.Q2 = coxph(Surv(tau-time, tau-Q, event = rep(1,nrow(dat.cox))) ~ Z1+Z2, data = dat.cox)
    basehaz.T = basehaz(fit.T, centered = FALSE)
    basehaz.Q2 = basehaz(fit.Q2, centered = FALSE)
    beta.T = coef(fit.T)
    beta.Q2 = coef(fit.Q2)
    time = dat$time
    Q = dat$Q
    Q2 = tau - Q
    Z = as.matrix(dat[,-(1:4)])
    
    Gtz = (baseS(tau-time, basehaz.Q2))^exp(Z %*% beta.Q2)
    Gqz = (baseS(tau-Q, basehaz.Q2))^exp(Z %*% beta.Q2)
    Fqz = 1-(baseS(Q, basehaz.T))^exp(Z %*% beta.T)
    mqz = diag(integral_F(Q, nu, Z, basehaz.T, beta.T))
    
    # compute the denominator of the estimator
    DDen1 = 1/Gtz
    DDen2 = Fqz/(Gqz*(1-Fqz))
    DDen3 = integral_JG_F(Q, time, Z, basehaz.Q2, basehaz.T, beta.Q2, beta.T, tau)
    
    den1 = sum(DDen1)
    den2 = sum(DDen2)
    den3 = sum(DDen3)

    
    # compute the numerator of the estimator
    NNum1 = nu(time)/Gtz
    NNum2 = mqz/(Gqz*(1-Fqz))
    NNum3 = integral_JG_m(Q, time, Z, basehaz.Q2, basehaz.T, beta.Q2, beta.T, tau, nu)
    
    num1 = sum(nu(time)/Gtz)
    num2 = sum(mqz/(Gqz*(1-Fqz)))
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
    Enutz = as.vector(integral_F(tau, nu, Z, basehaz.T, beta.T))
    
    DDen4 = 1/(1-Fqz)
    NNum4 = nu(time) + mqz/(1-Fqz)
    NNum5 = Enutz/(1-Fqz)
    
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
                beta2 = beta2
                )
           )
    
}
