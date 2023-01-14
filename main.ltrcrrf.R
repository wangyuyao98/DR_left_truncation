# The main code for loading data and compute the 'cf' estimator withLTRCforests::ltrcrrf() for estimating the nuisance parameters

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(survival)
library(LTRCforests)
library(cvTools)    # used when split the data into K folds
library(boot)
library(parallel)   # package for parallel computing
library(xtable)  # package for output tables in LaTeX format

source("src/functions.ltrcrrf.R")
source("src/estimator.ltrcrrf.R")
source("src/load_estimate.ltrcrrf.R")

dir.create("results.ltrcrrf")

# The given transformation of the survival time in defining the parameter of interest
t0 = 3





nu <- function(t,t0=3){
  # # identity function
  # result = t
  
  # indicator function
  result = as.numeric(t>t0)
  
  return(result)
}


alpha = 0.05

# parameters for LTRCforests::ltrcrrf
formula.T = Surv(Q, time, event = delta.T) ~ Z1+Z2
formula.Q2 = Surv(T2, Q2, event = delta.Q) ~ Z1+Z2
mtry = 2
ntree = 200

# trimming parameter
trim = 0.05

# number of folds for cross fitting
K = 10

# Tmodel.list = c("weibull1", "weibull1", "weibull2", "weibull2", "weibull1", "AFT2_weibull2", "AFT2_weibull2")
# Qmodel.list = c("weibull1", "weibull2", "weibull1", "weibull2", "weibull2_AFT2", "weibull1", "weibull2_AFT2")
Tmodel.list = c("weibull1")
Qmodel.list = c("Cox1")

itern = 500
n.boot = 100
cn = 7
mcn = cn*10   # number of tasks in mclapply
rn = ceiling(itern/mcn)   # number of repeating the mclapply

for(mi in 1:length(Tmodel.list)){
  set.seed(123)      # There is randomness in bootstrap and random forest. Set a seed for reproducibility.
  Tmodel = Tmodel.list[mi]
  Qmodel = Qmodel.list[mi]
  print(paste("Tmodel:", Tmodel, ";  Qmodel:", Qmodel))
  folder  = paste("T", Tmodel, "_Q", Qmodel, sep = "")  # the folder where the data sets are
  
  result = vector(mode = "list", length = rn)
  
  start_time <- Sys.time()
  system.time({
    for(i in 1:rn){
      print(paste("Iteration",i))
      
      index = (mcn*(i-1)+1):(mcn*(i-1) + min(mcn, itern-mcn*(i-1)))
      result[[i]] = mclapply(index, load_estimate.ltrcrrf, mc.cores = getOption("mc.cores", cn),
                             folder = folder, nu = nu, alpha = alpha, 
                             n.boot = n.boot, K = K, mtry = mtry, ntree = ntree, trim = trim)
      save(result, nu, alpha, n.boot, Tmodel, Qmodel, i,
           file = paste("results.ltrcrrf/",folder, "_ltrcrrf_R", n.boot, ".rda", sep = ""))
      
      now_time <- Sys.time()
      print(paste("Time elapsed", now_time - start_time))
    }
    
  })
  
  end_time <- Sys.time()
  end_time - start_time
  
  
  
  
  #################### code for organizing the results ###################
  result.mx.ltrcrrf = matrix(nrow = itern, ncol = length(result[[1]][[1]]$tab.ltrcrrf)+1)
  colnames(result.mx.ltrcrrf) = c("data_id", names(result[[1]][[1]]$tab.ltrcrrf))
  
  rn = length(result)
  mcn = length(result[[1]])
  for(i in 1:rn){
    index = (mcn*(i-1)+1):(mcn*(i-1) + min(mcn, itern-mcn*(i-1)))
    count_i  = length(index)
    for(j in 1:count_i){
      result.mx.ltrcrrf[mcn*(i-1)+j,] = c(data_id = result[[i]][[j]]$data.id, 
                                          result[[i]][[j]]$tab.ltrcrrf)
    }
  }
  
  # check if the results are complete
  if(length(unique(result.mx.ltrcrrf[,1])) == itern){
    print("The results are complete.")
  }else{
    print("The results are imcomplete.")
  }
  
  # load psi_true
  load(paste('datasets/',folder,'/data',1,'.rda', sep = ""))
  
  result.mx.ltrcrrf2 = result.mx.ltrcrrf[,-1]
  
  tab.ltrcrrf = cbind(bias = apply(result.mx.ltrcrrf2[,1:4], MARGIN = 2, mean) - psi_true,
                  percent_bias = (apply(result.mx.ltrcrrf2[,1:4], MARGIN = 2, mean) - psi_true)/psi_true*100,
                  sd = apply(result.mx.ltrcrrf2[,1:4], MARGIN = 2, sd),
                  se = apply(result.mx.ltrcrrf2[,5:8], MARGIN = 2, mean),
                  se.boot = apply(result.mx.ltrcrrf2[,9:12], MARGIN = 2, mean),
                  cover = apply(result.mx.ltrcrrf2[,13:16], MARGIN = 2, mean),
                  cover.boot = apply(result.mx.ltrcrrf2[,17:20], MARGIN = 2, mean))

  tab = tab.ltrcrrf
  
  tab2 = cbind(bias = sprintf("%.4f", tab[,1]),
               perc.bias =  sprintf("%.1f", tab[,2]),
               sd = sprintf("%.3f", tab[,3]),
               se_bootse = sprintf("%.3f/%.3f", tab[,4], tab[,5]),
               CP_bootCP = sprintf("%.3f/%.3f", tab[,6], tab[,7]))
  
  save(result, tab, tab2,
       nu, alpha, n.boot, Tmodel, Qmodel, i,
       file = paste("results.ltrcrrf/",folder, "_ltrcrrf_R", n.boot, ".rda", sep = ""))
  
  method.names = c("DR", "IPW.Q", "Reg.T1", "Reg.T2")
  xtable(cbind(NA, method.names, tab2),
         type = "latex", 
         digits = c(rep(3,8)),
         align = rep("c",8),
         file = "tab.tex")
  print(paste("Tmodel:", Tmodel, ";   Qmodel:", Qmodel))
  
  
}




# table that only shows boot SE and boot CP
tab3 = cbind(bias = sprintf("%.4f", tab[,1]),
             perc.bias =  sprintf("%.2f", tab[,2]),
             sd = sprintf("%.4f", tab[,3]),
             se_bootse = sprintf("%.4f", tab[,5]),
             CP_bootCP = sprintf("%.3f",  tab[,7]))
method.names = c("DR", "IPW.Q", "IPW.T1", "IPW.T2", "naive", "full")
print(xtable(cbind(NA, NA, NA, method.names, tab3),
             type = "latex", 
             digits = c(rep(3,10)),
             align = rep("c",10),
             file = "tab.tex"), 
      include.rownames = FALSE)
print(paste("Tmodel:", Tmodel, ";   Qmodel:", Qmodel))







