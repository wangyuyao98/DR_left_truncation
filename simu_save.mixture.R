# simulate and save the data sets

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(boot)
library(survival)
source("src/gen.R")

dir.create("datasets")
setwd("./datasets")

# Tmodel.list = c("weibull1", "weibull1", "weibull2", "weibull2", "weibull1", "AFT2_weibull2", "AFT2_weibull2")
# Qmodel.list = c("weibull1", "weibull2", "weibull1", "weibull2", "weibull2_AFT2", "weibull1", "weibull2_AFT2")
Tmodel.list = c("AFT2_weibull2")
Qmodel.list = c("Cox2_step2")

n = 1000
itern = 500
beta0.T_1 = 0.5     # For weibull2 -2.5; For AFT2 0.5
beta0.T_2 = -2.5    # For weibull2 -2.5; For AFT2 0.5
beta0.Q_1 = -0.5    # for Cox2 -0.5; for AFT2 0.5
beta0.Q_2 = -0.5     # for Cox2 -0.5; for AFT2 0.5
beta.T_1 = c(0.5,0.5)
beta.T_2 = c(0.5,0.5)
beta.Q_1 = c(0.5,0.5)
beta.Q_2 = c(0.5,0.5)
epsQ.min = -1
epsQ.max = 1
shape.T = 2

T.min = 1
Q.min = 0
Q.max = 4.5
tau = Q.max
multi = 20

Q2.min = tau - Q.max
Q2.max = tau - Q.min

# check if the simulated Q from AFT model is less than Q2.max
z1 = seq(-1,1, by = 0.01)
z1_sq = cbind(z1, z1^2)
Q2.min + exp(beta0.Q_1 + max(z1_sq%*%beta.Q_1) + epsQ.max)
Q2.min + exp(beta0.Q_1 + max(z1_sq%*%beta.Q_1) + epsQ.max) <= Q2.max
Q2.min + exp(beta0.Q_2 + max(z1_sq%*%beta.Q_2) + epsQ.max)
Q2.min + exp(beta0.Q_2 + max(z1_sq%*%beta.Q_2) + epsQ.max) <= Q2.max


### check the positivity assumption
### for F (T)
# weibull2_
min(pweibull(Q.max-T.min, shape = shape.T, 
              scale = exp(-1/shape.T * (beta0.T_1 + z1_sq%*%beta.T_1)),
              lower = F))
# _weibull2
min(pweibull(Q.max-T.min, shape = shape.T, 
             scale = exp(-1/shape.T * (beta0.T_2 + z1_sq%*%beta.T_2)),
             lower = F))
# AFT2_
max(log(Q.max - T.min) - beta0.T_1 - z1_sq%*%beta.T_1)
punif(max(log(Q.max - T.min) - beta0.T_1 - z1_sq%*%beta.T_1), min = -1, max = 1, lower = F)
# _AFT2
max(log(Q.max - T.min) - beta0.T_2 - z1_sq%*%beta.T_2)
punif(max(log(Q.max - T.min) - beta0.T_2 - z1_sq%*%beta.T_2), min = -1, max = 1, lower = F)


### for G (Q)
# Cox2_
(T.min/(Q.max-Q.min))^exp(max(beta0.Q_1+z1_sq%*%beta.Q_1))
# _Cox2
(T.min/(Q.max-Q.min))^exp(max(beta0.Q_1+z1_sq%*%beta.Q_2))
# AFT2_ 
max(log(Q.max - T.min) - beta0.Q_1 - z1_sq%*%beta.Q_1)
punif(max(log(tau - T.min) - beta0.Q_1 - z1_sq%*%beta.Q_1), 
      min = epsQ.min, max = epsQ.max, lower = F)
# _AFT2
max(log(tau - T.min) - beta0.Q_2 - z1_sq%*%beta.Q_2)
punif(max(log(tau - T.min) - beta0.Q_2 - z1_sq%*%beta.Q_2), 
      min = epsQ.min, max = epsQ.max, lower = F)
# _step2
pbeta((T.min-Q.min)/(Q.max-Q.min), 2,1, lower = T)




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


# The given transformation of the survival time in defining the parameter of interest
t0 = 3
nu <- function(t,t0=3){
    # # identity function
    # result = t
    
    # indicator function
    result = as.numeric(t>t0)
    
    return(result)
}



###  compute the true parameter of interest
nn = 5*10^5
set.seed(123)
gen_result = gen.mixture(nn, multi = 100, Tmodel, Qmodel, T.min, Q.min, Q.max, tau,
                         beta0.T_1, beta0.T_2, beta0.Q_1, beta0.Q_2,
                         beta.T_1, beta.T_2, beta.Q_1, beta.Q_2, 
                         shape.T, shape.Q, epsQ.min, epsQ.max)
dat = gen_result$dat
dat.full = gen_result$dat.full
psi_true = mean(nu(gen_result$dat.all$time))

# truncation rate
TR = 1- nrow(dat)/nrow(dat.full)
# censoring rate in observed data
CR.T = mean((dat.full$time > tau)[dat.full$time>dat.full$Q])
CR.Q = mean((dat.full$Q < 0)[dat.full$time>dat.full$Q])

hist(dat.full$time)

print(paste("average of", length(gen_result$dat.all$time), 
            "observations in full data;   psi_true = ", psi_true))
print(paste('Truncation rate in the whole population:', round(TR,3)))



# simulate the datasets
set.seed(123)
for(i in 1:itern){
    gen_result = gen.mixture(n, multi = 100, Tmodel, Qmodel, T.min, Q.min, Q.max, tau,
                             beta0.T_1, beta0.T_2, beta0.Q_1, beta0.Q_2,
                             beta.T_1, beta.T_2, beta.Q_1, beta.Q_2, 
                             shape.T, shape.Q, epsQ.min, epsQ.max)
    
    dat = gen_result$dat
    dat.full = gen_result$dat.full    # full data that corresponds to the observed sample with sample size n 
    
    save(dat, dat.full, psi_true, TR, CR.T, CR.Q,
         file = paste(folder,'/data',i,'.rda', sep = ""))
}

}

