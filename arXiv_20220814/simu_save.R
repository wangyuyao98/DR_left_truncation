# simulate and save the datasets

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(boot)
library(survival)
source("src/gen.R")

dir.create("datasets")
setwd("./datasets")

Tmodel.list = c("weibull1", "weibull1", "weibull2", "weibull2", "weibull1", "AFT2_weibull2", "AFT2_weibull2")
Qmodel.list = c("weibull1", "weibull2", "weibull1", "weibull2", "weibull2_AFT2", "weibull1", "weibull2_AFT2")

for(mi in 1:length(Tmodel.list)){
    
    Tmodel = Tmodel.list[mi]
    Qmodel = Qmodel.list[mi]
    print(paste("Tmodel:", Tmodel, ";  Qmodel:", Qmodel))

# # T model and Q model used for simulating data
# # Tmodel = 'weibull1'
# Tmodel = 'AFT2_weibull2'
# # Qmodel = 'weibull1'
# Qmodel = 'weibull2_AFT2'

# the folder to store the datasets
folder = paste("T", Tmodel, "_Q", Qmodel, sep = "")
dir.create(folder)

# parameters for datasets2
n = 1000
itern = 500
beta.T = c(0.3,0.5)
beta.Q = c(0.3,0.5)
beta.T2 = c(0.3,0.5,0.6,0.5)
beta.Q2 = c(0.3,0.5,0.6,0.5)
tau = 20
T.min = 5
Q.max = 8
multi = 20
beta0.T = -1
beta0.Q = -1
shape.T = 2
shape.Q = 2
epsT.sd = 1
epsQ.sd = 1



# The given transformation of the survival time in defining the parameter of interest
t0 = 7
nu <- function(t,t0=7){
    # # identity function
    # result = t
    
    # indicator function
    result = as.numeric(t>t0)
    
    return(result)
}



###  compute the true parameter of interest
nn = 5*10^5
set.seed(123)
gen_result = gen.mixture(nn, multi, Tmodel, Qmodel, T.min, Q.max, tau,
                         beta0.T, beta0.Q, beta.T, beta.Q, beta.T2, beta.Q2,
                         shape.T, shape.Q, epsT.sd)
dat = gen_result$dat
dat.full = gen_result$dat.full
psi_true = mean(nu(gen_result$dat.all$time))

# truncation rate
TR = 1- nrow(dat)/nrow(dat.full)
# censoring rate in observed data
CR.T = mean((dat.full$time > tau)[dat.full$time>dat.full$Q])
CR.Q = mean((dat.full$Q < 0)[dat.full$time>dat.full$Q])

print(paste("average of", length(gen_result$dat.all$time), 
            "observations in full data;   psi_true = ", psi_true))
print(paste('Truncation rate in the whole population:', round(TR,3)))
print(paste("CR.T = ", CR.T))
print(paste("CR.Q = ", CR.Q))


# simulate the datasets
set.seed(123)
for(i in 1:itern){
    gen_result = gen.mixture(n, multi, Tmodel, Qmodel, T.min, Q.max, tau,
                             beta0.T, beta0.Q, beta.T, beta.Q, beta.T2, beta.Q2,
                             shape.T, shape.Q, epsT.sd, epsQ.sd)
    
    dat = gen_result$dat
    dat.full = gen_result$dat.full    # full data that corresponds to the observed sample with sample size n 
    
    save(dat, dat.full, psi_true, TR, CR.T, CR.Q,
         file = paste(folder,'/data',i,'.rda', sep = ""))
}

}

