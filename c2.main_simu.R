# main function for generating data and compute the DR, IPW.Q, IPW.T1, IPW.T2 with 
# coxph for fitting the nuisance parameters, as well as the naive and full estimators
# Data generated with noninformative censoring on the residual time scale (c2 type)

# simulate and save the datasets
rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(boot)
library(survival)
library(parallel)   # package for parallel computing
library(LTRCforests)
library(cvTools)  # used for creating folds for cross-fitting
library(xtable)
source("src/c.gen.R")
source("src/functions.cox.R")
source("src/c2.estimator.R")
source("src/c.simu_estimate.cox_naive_full.R")

dir.create("c2.results")
dir.create("c2.datasets")

## parameters used to generate data 
load("seeds.rda")   #  load the seeds for generating data sets
Tmodel = "weibull1"
Qmodel = "Cox1"

# create the folder to store the info for the data sets
DGM = paste("X", Tmodel, "_Q", Qmodel, sep = "")
folder = paste("c2.datasets/", DGM, sep = "")
dir.create(folder)

n = 1000
beta.T = c(0.3,0.5)
beta.Q = c(0.3,0.5)
beta.T2 = c(-0.1,-0.1,0.3,0.3)
beta.Q2 = c(-0.1,-0.1,0.3,0.3)
T.min = 1
Q.min = 0
Q.max = 4.5
tau = Q.max
multi = 20
beta0.T = -2
beta0.Q = -2
shape.T = 2
shape.Q = 1.5
epsT.sd = 1
epsQ.sd = 1
C.min = 0
shape.C = 2
scale.C = 4
type.C = "c2"

alpha = 0.05
n.boot = 2
itern = 70
trim.C = 0.05   # bound the Sc(t) from below 



# The given transformation of the survival time in defining the parameter of interest
t0 = 3
nu <- function(t,t0=3){
    # # identity function
    # result = t
    
    # indicator function
    result = as.numeric(t>t0)
    
    return(result)
}

# compute the true parameter of interest and store it
nn = 5*10^5
set.seed(123)
gen_result = gen(nn, multi, Tmodel, Qmodel, T.min, Q.min, Q.max, tau,
                 beta0.T, beta0.Q, beta.T, beta.Q, beta.T2, beta.Q2,
                 shape.T, shape.Q, epsT.sd, epsQ.sd, C.min, shape.C, scale.C, type.C)
dat = gen_result$dat
dat.full = gen_result$dat.full
psi_true = mean(nu(gen_result$dat.full$TT))

TR = 1- nrow(dat)/nrow(dat.full)  # truncation rate
CR.full = mean(dat.full$delta == 0)  # censoring rate in full data
CR = mean(dat$delta ==0)  # censoring rate in observed data
TR
CR.full
CR
psi_true

save(psi_true, TR, CR.full, CR, file = paste(folder, "/psi_true.rda", sep = ""))

## dr-Cox-Cox, IPW.Q-Cox, Reg.T1-Cox, Reg.T2-Cox -----------------------------------------------
covariates.T = c("Z1","Z2")
# covariates.T = c("Z1sq","Z1Z2")
covariates.Q = c("Z1","Z2")
# covariates.Q = c("Z1sq","Z1Z2")

cn = 7
mcn = cn*10   # number of tasks in mclapply
rn = ceiling(itern/mcn)   # number of repeating the mclapply

result = vector(mode = "list", length = rn)

start_time <- Sys.time()
system.time({
    for(i in 1:rn){
        print(paste("Iteration",i))
        
        index = (mcn*(i-1)+1):(mcn*(i-1) + min(mcn, itern-mcn*(i-1)))
        result[[i]] = mclapply(seeds[index], c.simu_estimate.cox_naive_full, mc.cores = getOption("mc.cores", cn),
                               nu = nu, covariates.T = covariates.T, covariates.Q = covariates.Q,
                               trim.C = trim.C,  n.boot = n.boot,
                               n = n, multi = multi, Tmodel = Tmodel, Qmodel = Qmodel, 
                               T.min = T.min, Q.min = Q.min, Q.max = Q.max, tau = tau,
                               beta0.T = beta0.T, beta0.Q = beta0.Q, beta.T = beta.T, 
                               beta.Q = beta.Q, beta.T2 = beta.T2, beta.Q2 = beta.Q2,
                               shape.T = shape.T, shape.Q = shape.Q, epsT.sd = epsT.sd, epsQ.sd = epsQ.sd, 
                               C.min = C.min, shape.C = shape.C, scale.C = scale.C, type.C = type.C)
        save(result, nu, alpha, n.boot, trim.C, Tmodel, Qmodel, i,
             file = paste("c2.results/T", Tmodel, "_Q", Qmodel,
                          "_cox-T", paste(covariates.T, collapse = ""), "-Q", paste(covariates.Q, collapse = ""),
                          "_naive_full_R", n.boot, ".rda", sep = ""))
        
        now_time <- Sys.time()
        print(paste("Time elapsed", now_time - start_time))
    }
    
})
end_time <- Sys.time()
end_time - start_time

# Organizing the results 
result.mx.cox = matrix(nrow = itern, ncol = length(result[[1]][[1]]$tab.cox)+1)
colnames(result.mx.cox) = c("seed", names(result[[1]][[1]]$tab.cox))

result.mx.naive = matrix(nrow = itern, ncol = length(result[[1]][[1]]$tab.naive)+1)
colnames(result.mx.naive) = c("seed", names(result[[1]][[1]]$tab.naive))

result.mx.full = matrix(nrow = itern, ncol = length(result[[1]][[1]]$tab.full)+1)
colnames(result.mx.full) = c("seed", names(result[[1]][[1]]$tab.full))

rn = length(result)
mcn = length(result[[1]])
for(i in 1:rn){
    index = (mcn*(i-1)+1):(mcn*(i-1) + min(mcn, itern-mcn*(i-1)))
    count_i  = length(index)
    for(j in 1:count_i){
        result.mx.cox[mcn*(i-1)+j,] = c(seed = result[[i]][[j]]$seed, 
                                        result[[i]][[j]]$tab.cox)
        result.mx.naive[mcn*(i-1)+j,] = c(seed = result[[i]][[j]]$seed, 
                                          result[[i]][[j]]$tab.naive)
        result.mx.full[mcn*(i-1)+j,] = c(seed = result[[i]][[j]]$seed, 
                                         result[[i]][[j]]$tab.full)
    }
}

# check if the results are complete
if(length(unique(result.mx.cox[,1])) == itern){
    print("The results are complete.")
}else{
    print("The results are imcomplete.")
}


result.mx.cox2 = result.mx.cox[,-1]
result.mx.naive2 = result.mx.naive[,-1]
result.mx.full2 = result.mx.full[,-1]

qz = qnorm(alpha/2, lower = F)
est.mx = result.mx.cox2[,1:4]
se.mx = result.mx.cox2[,5:8]
bootse.mx = result.mx.cox2[,9:12]
cover.mx = (est.mx - qz*se.mx < psi_true & psi_true < est.mx + qz*se.mx)
bootcover.mx = (est.mx - qz*bootse.mx < psi_true & psi_true < est.mx + qz*bootse.mx)
tab.cox = cbind(bias = apply(est.mx, MARGIN = 2, mean) - psi_true,
                percent_bias = (apply(est.mx, MARGIN = 2, mean) - psi_true)/psi_true*100,
                sd = apply(est.mx, MARGIN = 2, sd),
                se = apply(se.mx, MARGIN = 2, mean),
                se.boot = apply(bootse.mx, MARGIN = 2, mean),
                cover = apply(cover.mx, MARGIN = 2, mean),
                cover.boot = apply(bootcover.mx, MARGIN = 2, mean))

est.naive = result.mx.naive2[,1]
se.naive = result.mx.naive2[,2]
bootse.naive = result.mx.naive2[,3]
cover.naive = est.naive - qz*se.naive < psi_true & psi_true < est.naive + qz*se.naive
bootcover.naive = est.naive - qz*bootse.naive < psi_true & psi_true < est.naive + qz*bootse.naive
tab.naive = cbind(bias = mean(est.naive) - psi_true,
                  percent_bias = (mean(est.naive) - psi_true)/psi_true*100,
                  sd = sd(est.naive),
                  se = mean(se.naive ),
                  se.boot = mean(bootse.naive),
                  cover = mean(cover.naive),
                  cover.boot = mean(bootcover.naive))

est.full = result.mx.full2[,1]
se.full = result.mx.full2[,2]
bootse.full = result.mx.full2[,3]
cover.full = est.full - qz*se.full < psi_true & psi_true < est.full + qz*se.full
bootcover.full = est.full - qz*bootse.full < psi_true & psi_true < est.full + qz*bootse.full
tab.full = cbind(bias = mean(est.full) - psi_true,
                 percent_bias = (mean(est.full) - psi_true)/psi_true*100,
                 sd = sd(est.full),
                 se = mean(se.full),
                 se.boot = mean(bootse.full),
                 cover = mean(cover.full),
                 cover.boot = mean(bootcover.full))

rownames(tab.naive) = "naive"
rownames(tab.full) = "full"
tab = rbind(tab.cox, tab.naive,  tab.full)

tab2 = cbind(bias = sprintf("%.4f", tab[,1]),
             perc.bias =  sprintf("%.1f", tab[,2]),
             sd = sprintf("%.3f", tab[,3]),
             se_bootse = sprintf("%.3f/%.3f", tab[,4], tab[,5]),
             CP_bootCP = sprintf("%.3f/%.3f", tab[,6], tab[,7]))


save(result, tab, tab2,
     nu, alpha, n.boot, Tmodel, Qmodel, i, covariates.T, covariates.Q,
     file = paste("c2.results/T", Tmodel, "_Q", Qmodel,
                  "_cox-T", paste(covariates.T, collapse = ""), "-Q", paste(covariates.Q, collapse = ""),
                  "_naive_full_R", n.boot, ".rda", sep = ""))



method.names = c("dr-Cox-Cox", "IPW.Q-Cox", "Reg.T1-Cox", "Reg.T2-Cox", "naive", "full")
rownames(tab2) = method.names
# print(
xtable(tab2,
       type = "latex", 
       digits = c(rep(3,6)),
       align = rep("c",6),
       file = "tab.tex")
#   , include.rownames = FALSE)
print(paste("Tmodel:", Tmodel, ";   Qmodel:", Qmodel))
print(paste("covariates.T:", paste(covariates.T, collapse = ", ")))
print(paste("covariates.Q:", paste(covariates.Q, collapse = ", ")))






# PL ------------------------------------------------------------------------------
source("c.simu_estimate_pl.R")

itern = 500
n.boot = 100

cn = 7
mcn = cn*10   # number of tasks in mclapply
rn = ceiling(itern/mcn)   # number of repeating the mclapply
result = vector(mode = "list", length = rn)

start_time <- Sys.time()
system.time({
    for(i in 1:rn){
        print(paste("Iteration",i))
        
        index = (mcn*(i-1)+1):(mcn*(i-1) + min(mcn, itern-mcn*(i-1)))
        result[[i]] = mclapply(seeds[index], c.simu_estimate.pl, mc.cores = getOption("mc.cores", cn),
                               t0 = t0, n.boot = n.boot,
                               n = n, multi = multi, Tmodel = Tmodel, Qmodel = Qmodel, 
                               T.min = T.min, Q.min = Q.min, Q.max = Q.max, tau = tau,
                               beta0.T = beta0.T, beta0.Q = beta0.Q, beta.T = beta.T, 
                               beta.Q = beta.Q, beta.T2 = beta.T2, beta.Q2 = beta.Q2,
                               shape.T = shape.T, shape.Q = shape.Q, epsT.sd = epsT.sd, epsQ.sd = epsQ.sd, 
                               C.min = C.min, shape.C = shape.C, scale.C = scale.C, type.C = type.C)
        save(result, t0, alpha, n.boot, Tmodel, Qmodel, i,
             file = paste("c2.results/T", Tmodel, "_Q", Qmodel,"_pl_R", n.boot, ".rda", sep = ""))
        
        now_time <- Sys.time()
        print(paste("Time elapsed", now_time - start_time))
    }
    
})
end_time <- Sys.time()
end_time - start_time


# Organizing the results 
result.mx.pl = matrix(nrow = itern, ncol = length(result[[1]][[1]]$tab.pl)+1)
colnames(result.mx.pl) = c("seed", names(result[[1]][[1]]$tab.pl))

rn = length(result)
mcn = length(result[[1]])
for(i in 1:rn){
    index = (mcn*(i-1)+1):(mcn*(i-1) + min(mcn, itern-mcn*(i-1)))
    count_i  = length(index)
    for(j in 1:count_i){
        result.mx.pl[mcn*(i-1)+j,] = c(seed = result[[i]][[j]]$seed, 
                                       result[[i]][[j]]$tab.pl)
    }
}

# check if the results are complete
if(length(unique(result.mx.pl[,1])) == length(unique(seeds))){
  print("The results are complete.")
}else{
  print("The results are incomplete.")
}


result.mx.pl2 = result.mx.pl[,-1]

qz = qnorm(alpha/2, lower = F)
est.mx = result.mx.pl2[,1]
bootse.mx = result.mx.pl2[,2]
bootcover.mx = (est.mx - qz*bootse.mx < psi_true & psi_true < est.mx + qz*bootse.mx)
tab= cbind(bias = mean(est.mx) - psi_true,
           percent_bias = (mean(est.mx) - psi_true)/psi_true*100,
           sd = sd(est.mx),
           se = NA,
           se.boot = mean(bootse.mx),
           cover = NA,
           cover.boot = mean(bootcover.mx))

tab2 = cbind(bias = sprintf("%.4f", tab[,1]),
             perc.bias =  sprintf("%.1f", tab[,2]),
             sd = sprintf("%.3f", tab[,3]),
             se_bootse = sprintf("- /%.3f", tab[,5]),
             CP_bootCP = sprintf("- /%.3f", tab[,7]))

save(result, tab, tab2,
     nu, alpha, n.boot, Tmodel, Qmodel, 
     file = paste("c2.results/T", Tmodel, "_Q", Qmodel,"_pl_R", n.boot, ".rda", sep = ""))


method.names = c("pl")
rownames(tab2) = method.names
# print(
xtable(tab2,
       type = "latex", 
       digits = c(rep(3,6)),
       align = rep("c",6),
       file = "tab.tex")
#   , include.rownames = FALSE)
print(paste("Tmodel:", Tmodel, ";   Qmodel:", Qmodel))





# KM ------------------------------------------------------------------------------
source("c.simu_estimate_KM.R")

itern = 500
n.boot = 100

cn = 7
mcn = cn*10   # number of tasks in mclapply
rn = ceiling(itern/mcn)   # number of repeating the mclapply
result = vector(mode = "list", length = rn)

set.seed(123)
start_time <- Sys.time()
system.time({
  for(i in 1:rn){
    print(paste("Iteration",i))
    
    index = (mcn*(i-1)+1):(mcn*(i-1) + min(mcn, itern-mcn*(i-1)))
    result[[i]] = mclapply(seeds[index], c.simu_estimate.KM, mc.cores = getOption("mc.cores", cn),
                           t0 = t0, n.boot = n.boot,
                           n = n, multi = multi, Tmodel = Tmodel, Qmodel = Qmodel, 
                           T.min = T.min, Q.min = Q.min, Q.max = Q.max, tau = tau,
                           beta0.T = beta0.T, beta0.Q = beta0.Q, beta.T = beta.T, 
                           beta.Q = beta.Q, beta.T2 = beta.T2, beta.Q2 = beta.Q2,
                           shape.T = shape.T, shape.Q = shape.Q, epsT.sd = epsT.sd, epsQ.sd = epsQ.sd, 
                           C.min = C.min, shape.C = shape.C, scale.C = scale.C, type.C = type.C)
    save(result, t0, alpha, n.boot, Tmodel, Qmodel, i,
         file = paste("c2.results/n", n,"_", DGM,"_KM_R", n.boot, ".rda", sep = ""))
    
    now_time <- Sys.time()
    print(paste("Time elapsed", now_time - start_time))
  }
  
})
end_time <- Sys.time()
end_time - start_time


# Organizing the results 
result.mx.KM = matrix(nrow = itern, ncol = length(result[[1]][[1]]$tab.KM)+1)
colnames(result.mx.KM) = c("seed", names(result[[1]][[1]]$tab.KM))

rn = length(result)
mcn = length(result[[1]])
for(i in 1:rn){
  index = (mcn*(i-1)+1):(mcn*(i-1) + min(mcn, itern-mcn*(i-1)))
  count_i  = length(index)
  for(j in 1:count_i){
    result.mx.KM[mcn*(i-1)+j,] = c(seed = result[[i]][[j]]$seed, 
                                   result[[i]][[j]]$tab.KM)
  }
}

result.mx.KM2 = result.mx.KM[,-1]

qz = qnorm(alpha/2, lower = F)
est.mx = result.mx.KM2[,1]
bootse.mx = result.mx.KM2[,2]
bootcover.mx = (est.mx - qz*bootse.mx < psi_true & psi_true < est.mx + qz*bootse.mx)
tab= cbind(bias = mean(est.mx) - psi_true,
           percent_bias = (mean(est.mx) - psi_true)/psi_true*100,
           sd = sd(est.mx),
           se = NA,
           se.boot = mean(bootse.mx),
           cover = NA,
           cover.boot = mean(bootcover.mx))

tab2 = cbind(bias = sprintf("%.4f", tab[,1]),
             perc.bias =  sprintf("%.1f", tab[,2]),
             sd = sprintf("%.3f", tab[,3]),
             se_bootse = sprintf("- /%.3f", tab[,5]),
             CP_bootCP = sprintf("- /%.3f", tab[,7]))

save(result, tab, tab2,
     nu, alpha, n.boot, Tmodel, Qmodel,
     file = paste("c2.results/n", n, "surv",t0, "_", DGM,"_KM_R", n.boot, ".rda", sep = ""))


method.names = c("KM")
rownames(tab2) = method.names
# print(
xtable(tab2,
       type = "latex", 
       digits = c(rep(3,6)),
       align = rep("c",6),
       file = "tab.tex")
#   , include.rownames = FALSE)
print(paste("Tmodel:", Tmodel, ";   Qmodel:", Qmodel))



