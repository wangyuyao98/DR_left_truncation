# The main code for loading data and computing the naive estimator and the full data estimator.

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
t0 = 3
nu <- function(t,t0=3){
    # # identity function
    # result = t
    
    # indicator function
    result = as.numeric(t>t0)
    
    return(result)
}

alpha = 0.05

# Tmodel.list = c("weibull1", "weibull1", "weibull2", "weibull2")
# Qmodel.list = c("Cox1", "Cox2", "Cox1", "Cox2")
Tmodel.list = c("weibull1")
Qmodel.list = c("Cox1")
# Tmodel.list = c("AFT2_weibull2")
# Qmodel.list = c("Cox2_step2")

covariates.T = c("Z1","Z2")
covariates.Q = c("Z1","Z2")


itern = 100
n.boot = 2
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
            result[[i]] = mclapply(index, load_estimate.naive_full, 
                                   mc.cores = getOption("mc.cores", cn),
                                   folder = folder, nu = nu, 
                                   alpha = alpha, n.boot = n.boot)
            save(result, nu, alpha, n.boot, Tmodel, Qmodel, i,
                 file = paste("results/",folder, "_naive_full_R", n.boot, ".rda", sep = ""))
            
            now_time <- Sys.time()
            print(paste("Time elapsed", now_time - start_time))
        }
        
    })
    
    end_time <- Sys.time()
    end_time - start_time
    
    
    #################### code for organizing the results ###################
    
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
            result.mx.naive[mcn*(i-1)+j,] = c(data_id = result[[i]][[j]]$data.id, 
                                              result[[i]][[j]]$tab.naive)
            result.mx.full[mcn*(i-1)+j,] = c(data_id = result[[i]][[j]]$data.id, 
                                             result[[i]][[j]]$tab.full)
        }
    }
    
    # load psi_true
    load(paste('datasets/',folder,'/data',1,'.rda', sep = ""))
    
    result.mx.naive2 = result.mx.naive[,-1]
    result.mx.full2 = result.mx.full[,-1]
    
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
    tab = rbind(tab.naive,  tab.full)
    
    tab2 = cbind(bias = sprintf("%.4f", tab[,1]),
                 perc.bias =  sprintf("%.1f", tab[,2]),
                 sd = sprintf("%.3f", tab[,3]),
                 se_bootse = sprintf("%.3f/%.3f", tab[,4], tab[,5]),
                 CP_bootCP = sprintf("%.3f/%.3f", tab[,6], tab[,7]))
    
    save(result, tab, tab2,
         nu, alpha, n.boot, Tmodel, Qmodel, i,
         file = paste("results/",folder, "_naive_full_R", n.boot, ".rda", sep = ""))
    
    
    method.names = c("naive", "full")
    # print(
    xtable(cbind(NA, method.names, tab2),
           type = "latex", 
           digits = c(rep(3,8)),
           align = rep("c",8),
           file = "tab.tex")
    #   , include.rownames = FALSE)
    print(paste("Tmodel:", Tmodel, ";   Qmodel:", Qmodel))
    
    
}

tab2
