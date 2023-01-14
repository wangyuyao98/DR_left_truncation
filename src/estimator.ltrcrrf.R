## function for computing the DR estimator using LTRCforests package
## The * notation in the code denote how time consuming it is to run that line. 
## More *'s means more time consuming.

# K-fold cross fiting is applyed
# K: number of folds in cross fitting
# formula.T: furmula for fitting model of T|Z
# formula.Q2: formula for fitting model for \tau-Q|Z
est_DR.ltrcrrf <- function(dat, formula.T, formula.Q2, nu, K, 
                           mtry = Inf, ntree = 100, trim = 1e-7){
    tau = max(dat$time)+1
    dat$Q2 = tau - dat$Q
    dat$T2 = tau - dat$time
    n = nrow(dat)
    dat$delta.T = rep(1, n)
    dat$delta.Q = rep(1, n)
    
    tau1 = min(dat$time)
    tau2 = max(dat$Q)
    
    # split the data into K folds 
    folds <- cvFolds(n, K)
    DDen1 = rep(NA, n)
    DDen2 = rep(NA, n)
    DDen3 = rep(NA, n)
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
        nk = nrow(dat.est)
        
        newdat = dat.est
        newdat$Q = 0
        newdat$T2 = 0
        newdat$time = tau
        newdat$Q2 = tau
        newdat$delta.T = 1
        newdat$delta.Q = 1
        
        # use dat.fit to estimate F and G
        T.obj <- ltrcrrf(formula = formula.T, data = dat.fit,
                         mtry = mtry, ntree = ntree)
        Q2.obj <- ltrcrrf(formula = formula.Q2, data = dat.fit,
                         mtry = mtry, ntree = ntree)
        
        jumps.T = sort(dat.fit$time)
        T.pred.dF = predictProb(object = T.obj, newdata = newdat,   
                                time.eval = jumps.T)   #**
        # row in $survival.probs corresponds to a time point in jumps.T, column corresponds to a subject
        
        # bound the estimates such that \hat{F}(\tau_2|Z_i) <= 1-trim for all i
        id.T = which(jumps.T <= tau2)
        T.pred.dF$survival.probs[id.T,] <- pmin(T.pred.dF$survival.probs[id.T,], 1-trim)
        
        
        jumps.Q = sort(dat.fit$Q)
        odr = order(tau-jumps.Q)
        Q.pred.dG = predictProb(object = Q2.obj, newdata = newdat, 
                                time.eval = (tau-jumps.Q)[odr])   #**
        # bound the estimates such that \hat{G}(\tau_1|Z_i) >= trim for all i
        id.Q = which(jumps.Q >= tau1)
        Q.pred.dG$survival.probs[id.Q, ] <- pmax(Q.pred.dG$survival.probs[id.Q, ], trim)
        
        Gtz = diag(predict.survProb.mx(Q.pred.dG, dat.est$T2))
        Gqz = diag(predict.survProb.mx(Q.pred.dG, dat.est$Q2))
        Fqz = 1-diag(predict.survProb.mx(T.pred.dF, dat.est$Q))
        mqz = int_nu_dF.LTRCforests(T.pred.dF, dat.est$Q, nu)
        
        Q.est = dat.est$Q
        time.est = dat.est$time
        n.est = nrow(dat.est)
        atrisk = matrix(nrow = n.est, ncol = length(jumps.Q))
        for(i in 1:n.est){
            atrisk[i,] = (Q.est[i] <= jumps.Q & jumps.Q < time.est[i])    # at risk indicator for subject i
        }
        
        
        Fvz.mx = 1-predict.survProb.mx(T.pred.dF, jumps.Q)
        Gvz.mx = t(Q.pred.dG$survival.probs)[,odr]
        
        intJG_F = int_FJ_dG.LTRCforests(Fvz.mx, Gvz.mx, atrisk, tau, nu, trim)
        intJG_m = int_JG_m.LTRCforests(T.pred.dF,jumps.Q, Fvz.mx, Gvz.mx, atrisk, tau, nu, trim)  #*

        # compute the denominator of the estimator
        DDen1[(ct+1):(ct+nk)] = 1/pmax(Gtz,trim)
        DDen2[(ct+1):(ct+nk)] = Fqz/pmax(Gqz*(1-Fqz), trim)
        DDen3[(ct+1):(ct+nk)] = int_FJ_dG.LTRCforests(Fvz.mx, Gvz.mx, atrisk, tau, nu, trim)
        
        # compute the numerator of the estimator
        NNum1[(ct+1):(ct+nk)] = nu(dat.est$time)/pmax(Gtz,trim)
        NNum2[(ct+1):(ct+nk)] = mqz/(pmax(Gqz, trim)*pmax(1-Fqz, trim))
        NNum3[(ct+1):(ct+nk)] = int_JG_m.LTRCforests(T.pred.dF,jumps.Q, Fvz.mx, Gvz.mx, atrisk, tau, nu, trim)
           #*
        
        # denominator and numerator for the IPSW estimator
        Enutz = as.vector(int_nu_dF.mx.LTRCforests(T.pred.dF, tau, nu))
        DDen4[(ct+1):(ct+nk)] = 1/pmax(1-Fqz, trim)
        NNum4[(ct+1):(ct+nk)] = nu(dat.est$time) + mqz/pmax(1-Fqz, trim)
        NNum5[(ct+1):(ct+nk)] = Enutz/pmax(1-Fqz, trim)
        
        ct = ct + nk
    }
    
    den1 = sum(DDen1)
    den2 = sum(DDen2)
    den3 = sum(DDen3)
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
    est_IPSW = num4/den4
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




