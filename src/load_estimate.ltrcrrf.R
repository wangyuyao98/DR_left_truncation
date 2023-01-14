# load data and then estimate using the DR, IPW.Q, IPW.T1, IPW.T2 with 
# LTRCforests::ltrcrrf for fitting the nuisance parameters

# alpha: (1-alpha)*100% CI is computed 
load_estimate.ltrcrrf <- function(iteration_id, folder, nu, alpha, n.boot,
                                  K, mtry, ntree, trim){
  
  load(paste('datasets/',folder,'/data',iteration_id,'.rda', sep = ""))
  
  qz = qnorm(alpha/2, lower = F)
  
  ### DR estimator with LTRCforests::ltrcrrf to fit the two parameters
  est.ltrcrrf = est_DR.ltrcrrf(dat, formula.T, formula.Q2, K, mtry, ntree, trim)
  estimate.ltrcrrf = c(DR = est.ltrcrrf$est_DR, 
                   IPQW = est.ltrcrrf$est_IPQW, 
                   IPSW = est.ltrcrrf$est_IPSW, 
                   IPSW2 = est.ltrcrrf$est_IPSW2)
  se.ltrcrrf = c(se_DR = est.ltrcrrf$se_DR,
                 se_IPQW = est.ltrcrrf$se_IPQW,
                 se_IPSW = est.ltrcrrf$se_IPSW,
                 se_IPSW2 = est.ltrcrrf$se_IPSW2)
  cover.ltrcrrf = (estimate.ltrcrrf-qz*se.ltrcrrf < psi_true) & (psi_true < estimate.ltrcrrf+qz*se.ltrcrrf)
  names(cover.ltrcrrf) = paste("cover_", names(estimate.ltrcrrf), sep = "")

  
  
  # SE and CI constructed using boostrap
  se_boot.ltrcrrf = rep(NA, 4)
  cover_boot.ltrcrrf = rep(NA, 4)
  # est_boot.ltrcrrf = boot(dat, estimate.ltrcrrf_boot, R = n.boot, 
  #                         nu = nu, K = K, mtry = mtry, ntree = ntree, trim = trim)
  # se_boot.ltrcrrf = apply(est_boot.ltrcrrf$t, MARGIN = 2, sd)
  # cover_boot.ltrcrrf = (estimate.ltrcrrf-qz*se_boot.ltrcrrf < psi_true) & (psi_true < estimate.ltrcrrf+qz*se_boot.ltrcrrf)
  names(se_boot.ltrcrrf) = paste("se_boot_", names(estimate.ltrcrrf), sep = "")
  names(cover_boot.ltrcrrf) = paste("cover_boot_", names(estimate.ltrcrrf), sep = "")
  
  
  tab.ltrcrrf = c(estimate.ltrcrrf, se.ltrcrrf, se_boot.ltrcrrf, cover.ltrcrrf, 
                  cover_boot.ltrcrrf)
  
  return(list(data.id = iteration_id,
              tab.ltrcrrf = tab.ltrcrrf))
}


estimate.ltrcrrf_boot <- function(dat, id, nu,
                                  K, mtry, ntree, trim){
  dat = dat[id,]
  est.ltrcrrf = est_DR.ltrcrrf(dat, formula.T, formula.Q2, K, mtry, ntree, trim)
  estimate.ltrcrrf = c(DR = est.ltrcrrf$est_DR, 
                       IPQW = est.ltrcrrf$est_IPQW, 
                       IPSW = est.ltrcrrf$est_IPSW, 
                       IPSW2 = est.ltrcrrf$est_IPSW2)
  
  return(estimate.ltrcrrf)
}


