# simulate and save the datasets
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(boot)
library(survival)
source("src/gen.R")

dir.create("datasets")
setwd("./datasets")

Tmodel.list = c("weibull1")        # "weibull1", "weibull2",
Qmodel.list = c("Cox1")       # "Cox1", "Cox2"

n = 1000
itern = 500
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

### check the positivity assumption
Z11 = seq(-1,1, by = 0.01)
Z1.extreme = rep(Z11, 2)
Z2.extreme = rep(c(0.5, -0.5), each = length(Z11))
Z.extreme = cbind(Z1 = Z1.extreme, Z2 = Z2.extreme)
Z.extreme.sq = cbind(Z1 = Z1.extreme, Z2 = Z2.extreme,
                     Z1_sq = Z1.extreme^2-1/3, Z1Z2 = Z1.extreme*Z2.extreme)
# Z.extreme.sq = cbind(Z1 = Z1.extreme, Z2 = Z2.extreme, 
#                      Z1_sq = sqrt(abs(Z1.extreme))-2/3, Z1Z2 = Z1.extreme*Z2.extreme)
## for F
# Weibull1
min(pweibull(Q.max-T.min, shape = shape.T, 
              scale = exp(-1/shape.T * (beta0.T + Z.extreme%*%beta.T)),
              lower = F))
# weibull2
min(pweibull(Q.max-T.min, shape = shape.T, 
         scale = exp(-1/shape.T * (beta0.T + Z.extreme.sq%*%beta.T2)),
         lower = F))
## for G
# Cox1
min((T.min/(Q.max-Q.min))^exp(Z.extreme%*%beta.Q))
# Cox2
min((T.min/(Q.max-Q.min))^exp(Z.extreme.sq%*%beta.Q2))






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
gen_result = gen(nn, multi, Tmodel, Qmodel, T.min, Q.min, Q.max, tau,
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

hist(dat.full$time)

print(paste("average of", length(gen_result$dat.all$time), 
            "observations in full data;   psi_true = ", psi_true))
print(paste('Truncation rate in the whole population:', round(TR,3)))



# simulate the datasets
set.seed(123)
for(i in 1:itern){
    gen_result = gen(n, multi, Tmodel, Qmodel, T.min, Q.min, Q.max, tau,
                     beta0.T, beta0.Q, beta.T, beta.Q, beta.T2, beta.Q2,
                     shape.T, shape.Q, epsT.sd)
    
    dat = gen_result$dat
    dat.full = gen_result$dat.full    # full data that corresponds to the observed sample with sample size n 
    
    save(dat, dat.full, psi_true, TR, CR.T, CR.Q,
         file = paste(folder,'/data',i,'.rda', sep = ""))
}

}

