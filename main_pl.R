# main function for computing the product-limit estimator for survival probabilities. 

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(boot)
library(survival)
library(parallel)   # package for parallel computing
library(xtable)  # package for output tables in LaTeX format

source("src/load_estimate_pl.R")

dir.create("results")

# The given transformation of the survival time in defining the parameter of interest
t0 = 3
alpha = 0.05

Tmodel.list = c("weibull1")
Qmodel.list = c("Cox1")
# Tmodel.list = c("weibull1", "weibull1", "weibull2", "weibull2")
# Qmodel.list = c("Cox1", "Cox2", "Cox1", "Cox2")

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
            result[[i]] = mclapply(index, load_estimate_pl, mc.cores = getOption("mc.cores", cn),
                                   folder = folder, t0 = t0,
                                   alpha = alpha, n.boot = n.boot)
            save(result, t0, alpha, n.boot, Tmodel, Qmodel, i,
                 file = paste("results/",folder, "_pl_R", n.boot, ".rda", sep = ""))
            
            now_time <- Sys.time()
            print(paste("Time elapsed", now_time - start_time))
        }
        
    })
    
    end_time <- Sys.time()
    end_time - start_time
    
    
    #################### code for organizing the results ###################
    
    result.mx_pl = matrix(nrow = itern, ncol = length(result[[1]][[1]]$tab_pl)+1)
    colnames(result.mx_pl) = c("data_id", names(result[[1]][[1]]$tab_pl))
    
    rn = length(result)
    mcn = length(result[[1]])
    for(i in 1:rn){
        index = (mcn*(i-1)+1):(mcn*(i-1) + min(mcn, itern-mcn*(i-1)))
        count_i  = length(index)
        for(j in 1:count_i){
            result.mx_pl[mcn*(i-1)+j,] = c(data_id = result[[i]][[j]]$data_id, 
                                            result[[i]][[j]]$tab_pl)
        }
    }
    
    # check if the results are complete
    if(length(unique(result.mx_pl[,1])) == itern){
        print("The results are complete.")
    }else{
        print("The results are imcomplete.")
    }
    
    # load psi_true
    load(paste('datasets/',folder,'/data',1,'.rda', sep = ""))
    
    result.mx_pl2 = result.mx_pl[,-1]
    tab = cbind(bias = mean(result.mx_pl2[,"est_pl"]) - psi_true,
                perc.bias = (mean(result.mx_pl2[,"est_pl"]) - psi_true)/psi_true*100,
                sd = sd(result.mx_pl2[,"est_pl"]),
                se = NA,
                bootse = mean(result.mx_pl2[,"se_boot_pl"]),
                CP = NA,
                bootCP = mean(result.mx_pl2[,"cover_boot_pl"]))
    
    tab2 = cbind(bias = sprintf("%.4f", mean(result.mx_pl2[,"est_pl"]) - psi_true),
                 perc.bias =  sprintf("%.1f", (mean(result.mx_pl2[,"est_pl"]) - psi_true)/psi_true*100),
                 sd = sprintf("%.3f", sd(result.mx_pl2[,"est_pl"])),
                 se_bootse = sprintf("- /%.3f", mean(result.mx_pl2[,"se_boot_pl"])),
                 CP_bootCP = sprintf("- /%.3f", mean(result.mx_pl2[,"cover_boot_pl"]))) 
    
    
    save(result, result.mx_pl2, tab, tab2, t0, alpha, n.boot, Tmodel, Qmodel, i,
         file = paste("results/",folder, "_pl_R", n.boot, ".rda", sep = ""))
    
    method.names = c("pl")
    print(xtable(cbind(NA, method.names, tab2),
           type = "latex", 
           digits = c(rep(3,8)),
           align = rep("c",8),
           file = "tab.tex"))
    #   , include.rownames = FALSE)
    print(paste("Tmodel:", Tmodel, ";   Qmodel:", Qmodel))
    
    
}


