# main function for 

# main function for generating data and compute the DR, IPW.Q, IPW.T1, IPW.T2 with 
# coxph for fitting the nuisance parameters, as well as the naive and full estimators
# Data generated with noninformative censoring on the original scale (c1 type)

# simulate and save the datasets
rm(list = ls())

# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(boot)
library(survival)
library(LTRCforests)
library(cvTools)  # used for creating folds for cross-fitting
source("c1.gen_Xcox.R")
source("c1.estimator.R")
source("c1.functions.ltrcrrf.R")

load("seeds_input.rda")

## parameters used to generate data 
Tmodel = "weibull1"
Qmodel = "Cox1"

n = 1000
beta.T = c(0.3,0.5)
beta.Q = c(0.3,0.5)
beta.T2 = c(-0.1,-0.1,0.3,0.3)
beta.Q2 = c(-0.1,-0.1,0.3,0.3)
T.min = 1
Q.min = 0
Q.max = 4.5
C.min = 1
tau = Q.max
multi = 20
beta0.T = -2
beta0.Q = -2
shape.T = 2
shape.Q = 1.5
epsT.sd = 1
epsQ.sd = 1
shape.C = 2
scale.C = 4

alpha = 0.05
nb.each = 20
trim.C = 0.05   # bound the Sc(t) from below 
type.C = "c1"

# parameters used for fitting models from ltrcrrf
K = 10  # number of folds for cross fitting
mtry = 2
ntree = 100
trim = 0.05   # bound for trimming the estimated probabilities


# The given transformation of the survival time in defining the parameter of interest
tt0 = 3
t0 = tt0
nu <- function(t,t0=tt0){
    # # identity function
    # result = t
    
    # indicator function
    result = as.numeric(t>t0)
    
    return(result)
}


## simulate data
set.seed(seed)
gen_result = gen(n, multi, Tmodel, Qmodel, T.min, Q.min, Q.max, tau,
                 beta0.T, beta0.Q, beta.T, beta.Q, beta.T2, beta.Q2,
                 shape.T, shape.Q, epsT.sd, epsQ.sd, C.min, shape.C, scale.C)
dat = gen_result$dat    # the observed data



## cf-RF-RF, IPW.Q-RF, Reg.T1-RF, Reg.T2-RF 
covariates.X = c("Z1","Z2")
covariates.Q = c("Z1","Z2")

bootresult.ltrcrrf = c1.splitboot_estSurv.ltrcrrf(seed.b,
                             n.boot = nb.each, dat = dat, t0 = t0,
                             covariates.X = covariates.X, covariates.Q = covariates.Q,
                             K = K, mtry = mtry, ntree = ntree, trim = trim, trim.C = trim.C)

save(bootresult.ltrcrrf, seed, seed.b, t0, mtry, ntree, trim, trim.C, K,
     file = "result.rda")

