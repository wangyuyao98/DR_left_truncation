# load data and then estimate using the DR, IPW.Q, IPW.T1, IPW.T2 with 
# coxph for fitting the nuisance parameters, 
# and naive and full estimators

# alpha: (1-alpha)*100% CI is computed 
load_estimate <- function(iteration_id, folder, nu, tau, alpha, n.boot){
  
  load(paste('datasets/',folder,'/data',iteration_id,'.rda', sep = ""))
  
  qz = qnorm(alpha/2, lower = F)
  
  
  ### DR estimator with coxph to fit the two parameters
  # DR estimator based on Cox models without cross fitting
  est.cox = est_DR.cox(dat, tau, nu)
  estimate.cox = c(DR = est.cox$est_DR, 
                   IPQW = est.cox$est_IPQW, 
                   IPSW = est.cox$est_IPSW, 
                   IPSW2 = est.cox$est_IPSW2)
  se.cox = c(se_DR = est.cox$se_DR,
             se_IPQW = est.cox$se_IPQW,
             se_IPSW = est.cox$se_IPSW,
             se_IPSW2 = est.cox$se_IPSW2)
  cover.cox = (estimate.cox-qz*se.cox < psi_true) & (psi_true < estimate.cox+qz*se.cox)
  names(cover.cox) = paste("cover_", names(estimate.cox), sep = "")

  # SE and CI constructed using boostrap
  est_boot.cox = boot(dat, estimate.cox_boot, R = n.boot, nu = nu, tau = tau)
  se_boot.cox = apply(est_boot.cox$t, MARGIN = 2, sd)
  cover_boot.cox = (estimate.cox-qz*se_boot.cox < psi_true) & (psi_true < estimate.cox+qz*se_boot.cox)
  names(se_boot.cox) = paste("se_boot_", names(estimate.cox), sep = "")
  names(cover_boot.cox) = paste("cover_boot_", names(estimate.cox), sep = "")
  
  tab.cox = c(estimate.cox, se.cox, se_boot.cox, cover.cox, cover_boot.cox)
  
  
  
  ### naive estimator
  n = nrow(dat)
  est.naive = mean(nu(dat$time))
  se.naive = sd(nu(dat$time))/sqrt(n)
  cover.naive = (est.naive-qz*se.naive < psi_true) & (psi_true < est.naive+qz*se.naive)
  names(est.naive) = "est.naive"
  names(se.naive) = "se.naive"
  names(cover.naive) = "cover.naive"
  
  # SE and CI constructed using boostrap
  est_boot.naive = boot(dat, estimator.naive_boot, R = n.boot, nu = nu)
  se_boot.naive = apply(est_boot.naive$t, MARGIN = 2, sd)
  cover_boot.naive = (est.naive-qz*se_boot.naive < psi_true) & (psi_true < est.naive+qz*se_boot.naive)
  names(se_boot.naive) = "se_boot_naive"
  names(cover_boot.naive) = "cover_boot_naive"
  
  tab.naive = c(est.naive, se.naive, se_boot.naive, cover.naive, cover_boot.naive)
  
  
  
  ### full data estimator
  n.full = nrow(dat.full)
  est.full = mean(nu(dat.full$time))
  se.full = sd(nu(dat.full$time))/sqrt(n.full)
  cover.full = (est.full-qz*se.full < psi_true) & (psi_true < est.full+qz*se.full)
  names(est.full) = "est.full"
  names(se.full) = "se.full"
  names(cover.full) = "cover.full"
  
  # SE and CI constructed using boostrap
  est_boot.full = boot(dat.full, estimator.full_boot, R = n.boot, nu = nu)
  se_boot.full = apply(est_boot.full$t, MARGIN = 2, sd)
  cover_boot.full = (est.full-qz*se_boot.full < psi_true) & (psi_true < est.full+qz*se_boot.full)
  names(se_boot.full) = "se_boot_full"
  names(cover_boot.full) = "cover_boot_full"
  
  tab.full = c(est.full, se.full, se_boot.full, cover.full, cover_boot.full)
  
  
  
  return(list(tab.cox = tab.cox, 
              tab.naive = tab.naive,
              tab.full = tab.full,
              data.id = iteration_id))
}


estimate.cox_boot <- function(dat, id, nu, tau){
  dat = dat[id,]
  result.cox = est_DR.cox(dat, tau, nu)
  est.cox = c(DR = result.cox$est_DR, 
                   IPQW = result.cox$est_IPQW, 
                   IPSW = result.cox$est_IPSW, 
                   IPSW2 = result.cox$est_IPSW2)
  
  return(est.cox)
}


estimator.naive_boot <- function(dat, id, nu){
  dat = dat[id,]
  est.naive = mean(nu(dat$time))
  
  return(est.naive)
}

estimator.full_boot <- function(dat.full, id, nu){
  dat.full = dat.full[id,]
  est.full = mean(nu(dat.full$time))
  
  return(est.full)
}



