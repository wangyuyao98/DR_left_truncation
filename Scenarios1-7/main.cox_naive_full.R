# main function for loading data and compute the DR, IPW.Q, IPW.T1, IPW.T2 with 
# coxph for fitting the nuisance parameters, 
# and naive and full estimators

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(boot)
library(survival)
library(parallel)   # package for parallel computing
library(xtable)  # package for output tables in LaTeX format

source("src/load_estimate.cox_naive_full.R")
source("src/estimator.cox.R")
source("src/functions.cox.R")

dir.create("results")

# The given transformation of the survival time in defining the parameter of interest
t0 = 7
nu <- function(t,t0=7){
    # # identity function
    # result = t
    
    # indicator function
    result = as.numeric(t>t0)
    
    return(result)
}


tau = 20
alpha = 0.05

Tmodel.list = c("weibull1", "weibull1", "weibull2", "weibull2", "weibull1", "AFT2_weibull2", "AFT2_weibull2")
Qmodel.list = c("weibull1", "weibull2", "weibull1", "weibull2", "weibull2_AFT2", "weibull1", "weibull2_AFT2")

itern = 500
n.boot = 100
cn = 7
mcn = cn*10   # number of tasks in mclapply
rn = ceiling(itern/mcn)   # number of repeating the mclapply

for(mi in 1:length(Tmodel.list)){
  set.seed(123) # There is randomness in bootstrap. Set a seed for reproducibility.
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
            result[[i]] = mclapply(index, load_estimate, mc.cores = getOption("mc.cores", cn),
                                   folder = folder, nu = nu, tau = tau, alpha = alpha, 
                                   n.boot = n.boot)
            save(result, nu, tau, alpha, n.boot, Tmodel, Qmodel, i,
                 file = paste("results/",folder, "_cox_naive_full_R", n.boot, ".rda", sep = ""))
            
            now_time <- Sys.time()
            print(paste("Time elapsed", now_time - start_time))
        }
        
    })
    
    end_time <- Sys.time()
    end_time - start_time
    
    
    
    
    #################### code for organizing the results ###################
    
    result.mx.cox = matrix(nrow = itern, ncol = length(result[[1]][[1]]$tab.cox)+1)
    colnames(result.mx.cox) = c("data_id", names(result[[1]][[1]]$tab.cox))
    
    result.mx.naive = matrix(nrow = itern, ncol = length(result[[1]][[1]]$tab.naive)+1)
    colnames(result.mx.naive) = c("data_id", names(result[[1]][[1]]$tab.naive))
    
    result.mx.full = matrix(nrow = itern, ncol = length(result[[1]][[1]]$tab.full)+1)
    colnames(result.mx.full) = c("data_id", names(result[[1]][[1]]$tab.full))
    
    rn = length(result)
    mcn = length(result[[1]])
    for(i in 1:rn){
        index = (mcn*(i-1)+1):(mcn*(i-1) + min(mcn, itern-mcn*(i-1)))
        count_i  = length(index)
        for(j in 1:count_i){
            result.mx.cox[mcn*(i-1)+j,] = c(data_id = result[[i]][[j]]$data.id, 
                                            result[[i]][[j]]$tab.cox)
            result.mx.naive[mcn*(i-1)+j,] = c(data_id = result[[i]][[j]]$data.id, 
                                              result[[i]][[j]]$tab.naive)
            result.mx.full[mcn*(i-1)+j,] = c(data_id = result[[i]][[j]]$data.id, 
                                             result[[i]][[j]]$tab.full)
        }
    }
    
    # check if the results are complete
    if(length(unique(result.mx.cox[,1])) == itern){
        print("The results are complete.")
    }else{
        print("The results are imcomplete.")
    }
    
    # load psi_true
    load(paste('datasets/',folder,'/data',1,'.rda', sep = ""))
    
    result.mx.cox2 = result.mx.cox[,-1]
    result.mx.naive2 = result.mx.naive[,-1]
    result.mx.full2 = result.mx.full[,-1]
    
    tab.cox = cbind(bias = apply(result.mx.cox2[,1:4], MARGIN = 2, mean) - psi_true,
                    percent_bias = (apply(result.mx.cox2[,1:4], MARGIN = 2, mean) - psi_true)/psi_true*100,
                    sd = apply(result.mx.cox2[,1:4], MARGIN = 2, sd),
                    se = apply(result.mx.cox2[,5:8], MARGIN = 2, mean),
                    se.boot = apply(result.mx.cox2[,9:12], MARGIN = 2, mean),
                    cover = apply(result.mx.cox2[,13:16], MARGIN = 2, mean),
                    cover.boot = apply(result.mx.cox2[,17:20], MARGIN = 2, mean))
    
    tab.naive = cbind(bias = mean(result.mx.naive2[,1]) - psi_true,
                      percent_bias = (mean(result.mx.naive2[,1]) - psi_true)/psi_true*100,
                      sd = sd(result.mx.naive2[,1]),
                      se = mean(result.mx.naive2[,2]),
                      se.boot = mean(result.mx.naive2[,3]),
                      cover = mean(result.mx.naive2[,4]),
                      cover.boot = mean(result.mx.naive2[,5]))
    
    tab.full = cbind(bias = mean(result.mx.full2[,1]) - psi_true,
                     percent_bias = (mean(result.mx.full2[,1]) - psi_true)/psi_true*100,
                     sd = sd(result.mx.full2[,1]),
                     se = mean(result.mx.full2[,2]),
                     se.boot = mean(result.mx.full2[,3]),
                     cover = mean(result.mx.full2[,4]),
                     cover.boot = mean(result.mx.full2[,5]))
    
    rownames(tab.naive) = "naive"
    rownames(tab.full) = "full"
    tab = rbind(tab.cox, tab.naive,  tab.full)
    
    tab2 = cbind(bias = sprintf("%.4f", tab[,1]),
                 perc.bias =  sprintf("%.1f", tab[,2]),
                 sd = sprintf("%.3f", tab[,3]),
                 se_bootse = sprintf("%.3f/%.3f", tab[,4], tab[,5]),
                 CP_bootCP = sprintf("%.3f/%.3f", tab[,6], tab[,7]))
    
    save(result, tab, tab2,
         nu, tau, alpha, n.boot, Tmodel, Qmodel, i,
         file = paste("results/",folder, "_cox_naive_full_R", n.boot, ".rda", sep = ""))
    
    
    method.names = c("DR", "IPW.Q", "Reg.T1", "Reg.T2", "naive", "full")
    print(xtable(cbind(NA, method.names, tab2),
                 type = "latex", 
                 digits = c(rep(3,8)),
                 align = rep("c",8),
                 file = "tab.tex"), 
          include.rownames = FALSE)
    print(paste("Tmodel:", Tmodel, ";   Qmodel:", Qmodel))
    
    
}




# table that only shows boot SE and boot CP
tab3 = cbind(bias = sprintf("%.4f", tab[,1]),
             perc.bias =  sprintf("%.1f", tab[,2]),
             sd = sprintf("%.3f", tab[,3]),
             se_bootse = sprintf("%.3f", tab[,5]),
             CP_bootCP = sprintf("%.3f",  tab[,7]))
method.names = c("DR", "IPW.Q", "Reg.T1", "Reg.T2", "naive", "full")
print(xtable(cbind(NA, NA, NA, method.names, tab3),
             type = "latex", 
             digits = c(rep(3,10)),
             align = rep("c",10),
             file = "tab.tex"), 
      include.rownames = FALSE)
print(paste("Tmodel:", Tmodel, ";   Qmodel:", Qmodel))

