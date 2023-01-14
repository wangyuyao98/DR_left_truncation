## functions used to compute the PL estimate

c.simu_estimate.pl <- function(seed, t0, n.boot,
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
    
    # estimator and bootstrap SE
    result_pl = boot(dat, est_pl.boot, R = n.boot, t0 = t0)
    est_pl = result_pl$t0
    se_boot_pl = sd(result_pl$t)
    
    tab_pl = c(est_pl = est_pl, se_boot_pl = se_boot_pl)
    
    return(list(seed = seed, 
                tab.pl= tab_pl))
    
}

est_pl.boot <- function(dat, id, t0){
    dat = dat[id,]
    
    fit.KM.T = survfit(Surv(Q, X, delta)~1, data = dat)
    sf = stepfun(fit.KM.T$time, c(1,fit.KM.T$surv))
    
    return(sf(t0))
    
}


est_pl <- function(t0, dat, alpha, n.boot){
    # estimator and bootstrap SE
    result_pl = boot(dat, est_pl.boot, R = n.boot, t0 = t0)
    est_pl = result_pl$t0
    se_boot_pl = sd(result_pl$t)
    
    return(list(est_pl = est_pl, se_boot_pl = se_boot_pl))
}







