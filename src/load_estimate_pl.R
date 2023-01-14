# product limit estimator for survival probabilities. 

load_estimate_pl <- function(iteration_id, folder, t0, alpha, n.boot){
    
    load(paste('datasets/',folder,'/data',iteration_id,'.rda', sep = ""))
    qz = qnorm(alpha/2, lower = F)
    n = nrow(dat)
    
    # estimator and bootstrap SE
    result_pl = boot(dat, est_pl.boot, R = n.boot, t0 = t0)
    est_pl = result_pl$t0
    se_boot_pl = sd(result_pl$t)
    cover_boot_pl = (est_pl-qz*se_boot_pl < psi_true) & (psi_true < est_pl+qz*se_boot_pl)
    
    tab_pl = c(est_pl = est_pl, se_boot_pl = se_boot_pl, cover_boot_pl = cover_boot_pl)
    
    return(list(data_id = iteration_id, 
                tab_pl= tab_pl))
    
}



est_pl.boot <- function(dat, id, t0){
    
    dat = dat[id,]
    
    fit.KM.T = survfit(Surv(Q, time, rep(1,nrow(dat)))~1, data = dat)
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




