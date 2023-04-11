## The main code for analyzing the CNS lymphoma data ("CNS_data.txt") from a study of methotrexate-based chemotherapy (Wang et al., 2015).
## The data set is from the supplementary material of Vakulenko-Lagun et al.(2021).

rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(boot)
library(survival)
library(parallel)   # package for parallel computing
library(tranSurv)    # conditional Kendall's tau test
library(CoxR2)    # compute the R2 measures
library(LTRCforests)
library(cvTools)  # used for creating folds for cross-fitting

source("src/c1.estimator.R")
source("src/c1.functions.cox.R")
source("src/c1.functions.ltrcrrf.R")

dir.create("CNS_data/c1.results")
dir.create("CNS_data/c1.plots")

# load and process data ------------------------------------------------------
# load data
d <- read.table(file="CNS_data/CNS_data.txt", sep="\t", header=TRUE)
d$chemo.i <- ifelse(d$chemo==1, 1, 0)

dat0 = d  # original data set
dat0$X = dat0$OS
dat0$delta = dat0$death

s <- d[d$relapse==1,] # select only those who progressed to relapse
dim(s)  # 98*12

# Check dependency: conditional Kendall's tau test
cKendall(s$PFS, s$OS, s$death)    # tau = -0.0384, p-value = 0.5684

obs10 <- min(s$OS[s$chemo.i==1])
s <- s[s$OS!=obs10,] # As in the analysis in Vakulenko-Lagun et al.(2021), omitting an observation that zeroes out the risk set in the beginning of follow-up
dim(s)

# Note the censoring of this data set is c1 type
dat <- s
dat$X = dat$OS
dat$Q = dat$PFS
dat$delta = dat$death

tau = max(dat$X) + 1
CR = mean(1-dat$delta)
CR   # censoring rate

tau
min(dat$X[dat$delta == 1])


covariates = c("chemo.i", "radiation")
covariates.X = covariates
covariates.Q = covariates

dat$delta.1 = rep(1, nrow(dat))
formula.X = formula(paste("Surv(Q, X, delta.1) ~ ", 
                          paste(covariates.X, collapse = "+"), collapse = ""))
formula.Q2 = formula(paste("Surv(tau-X, tau-Q, delta.1) ~ ", 
                           paste(covariates.Q, collapse = "+"), collapse = ""))


# fit Cox-PH model for T|Z
fit.X = coxph(formula.X, data = dat)
# fit Cox-PH model for (tau-Q)|Z
fit.Q2 = coxph(formula.Q2, data = dat)

coxr2(fit.X)  # R2 = 0.0346
coxr2(fit.Q2)   # R2 = 0.344

fit.X
fit.Q2


### compute the estimators ------------------------------------------------------------=
alpha = 0.05
qz = qnorm(alpha/2, lower = FALSE)
trim.C = 1e-7    # bound the Sc(t) from below 
n.boot = 100

t0.list = seq(min(dat$X[dat$delta==1])-1, max(dat$X)+1, by = 1)

# PL estimator --------------------------------------------------------------------
plfit.T = survfit(Surv(Q, X, delta)~1, data = dat)
plsf.T = stepfun(plfit.T$time, c(1, plfit.T$surv))
est.pl = plsf.T(t0.list)

# KM estimator that ignores left truncation  -----------------------------------------
kmfit.T = survfit(Surv(X, delta)~1, data = dat)
kmsf.T = stepfun(kmfit.T$time, c(1, kmfit.T$surv))
est.km = kmsf.T(t0.list)

# naive estimator -----------------------------------------------------------------
est_naive = sapply(t0.list, c1.estSurv_naive, dat = dat)
names(est_naive) = t0.list

# dr-Cox-Cox, IPW.Q-Cox, Reg.T1-Cox, Reg.T2-Cox -------------------------------------
set.seed(123)
start_time <- Sys.time()
system.time({
  result = c1.estSurv_DR.cox(t0.list, dat, covariates.X, covariates.Q, trim.C)
})
now_time <- Sys.time()
print(paste("Time elapsed", now_time - start_time))
save(result, t0.list, covariates.X, covariates.Q, trim.C, n.boot,
     file = paste("CNS_data/c1.results/CNS_n",nrow(dat),".cox_trimC", trim.C, ".rda", sep = ""))


## cf-RF-RF ------------------------------------------------------------------------
# parameters used for ltrcrrf
K = 10  # number of folds for cross fitting
mtry = 2
ntree = 2000
trim = 0.05   # bound for trimming the estimated probabilities

alpha = 0.05
qz = qnorm(alpha/2, lower = FALSE)

set.seed(123)
start_time <- Sys.time()
system.time({
  result.rf = c1.estSurv_cf.ltrcrrf(t0.list, dat, covariates.X, covariates.Q, K,
                                 mtry, ntree, trim, trim.C)
})
now_time <- Sys.time()
print(paste("Time elapsed", now_time - start_time))
save(result.rf, t0.list, covariates.X, covariates.Q, trim.C,
     file = paste("CNS_data/c1.results/CNS_n",nrow(dat),".ltrcrrf_ntree", ntree, ".rda", sep = ""))


est_cf.RF = result.rf$est_DR
se_cf.RF = result.rf$se_DR

  


### percentile method for getting CI's --------------------------------------------------
cn = 8
nb.total = 2000
ns = cn
set.seed(123)
seeds.boot = sample(10000, size = ns)
nb = ceiling(nb.total/ns)

# load(paste("CNS_data/c1.results/CNS_n",nrow(dat),".cox_trimC", trim.C, ".rda", sep = ""))
# load(paste("CNS_data/c1.results/CNS_n",nrow(dat),".ltrcrrf_ntree", ntree, ".rda", sep = ""))

est_cf_RF = result.rf$est_DR
# finite sample correction so that the survival curve is decreasing
for(i in 1:length(est_cf_RF)){
  est_cf_RF[i] = min(est_cf_RF[1:i])
}

est = cbind(cf_RF = est_cf_RF, 
            DR = result$est_DR, 
            IPQW = result$est_IPQW, 
            IPSW = result$est_IPSW, 
            IPSW2 = result$est_IPSW2,
            naive = est_naive)
rownames(est) = t0.list


# dr-Cox-Cox, IPW.Q-Cox, Reg.T1-Cox, Reg.T2-Cox 
start_time <- Sys.time()
system.time({
  bootresult.cox = mclapply(seeds.boot, c1.splitboot_estSurv_DR.cox, mc.cores = getOption("mc.cores", cn),
                            n.boot = nb, dat = dat, t0 = t0.list,
                            covariates.X = covariates.X, covariates.Q = covariates.Q,
                            trim.C = trim.C)
})
now_time <- Sys.time()
print(paste("Time elapsed", now_time - start_time))
save(bootresult.cox, t0.list, covariates.X, covariates.Q, trim.C, n.boot,
     file = paste("CNS_data/c1.results/CNS_n",nrow(dat),".sbootcox",nb.total, "_trimC", trim.C, ".rda", sep = ""))

# load(paste("CNS_data/c1.results/CNS_n",nrow(dat),".sbootcox",nb.total, "_trimC", trim.C, ".rda", sep = ""))

# cf-RF-RF
start_time <- Sys.time()
system.time({
  result.rf = c1.estSurv_cf.ltrcrrf(t0.list, dat, covariates.X, covariates.Q, K,
                                    mtry, ntree, trim, trim.C)
  bootresult.ltrcrrf = mclapply(seeds.boot, c1.splitboot_estSurv.ltrcrrf, mc.cores = getOption("mc.cores", cn),
                            n.boot = nb, dat = dat, t0 = t0.list,
                            covariates.X = covariates.X, covariates.Q = covariates.Q,
                            K = K, mtry = mtry, ntree = ntree, trim = trim, trim.C = trim.C)
})
now_time <- Sys.time()
print(paste("Time elapsed", now_time - start_time))
save(bootresult.ltrcrrf, t0.list, covariates.X, covariates.Q, K, mtry, ntree, trim, trim.C, n.boot,
     file = paste("CNS_data/c1.results/CNS_n",nrow(dat),".sbootRF",nb.total, "_trimC", trim.C, ".rda", sep = ""))

# load(paste("CNS_data/c1.results/CNS_n",nrow(dat),".sbootRF",nb.total, "_trimC", trim.C, ".rda", sep = ""))


m = length(t0.list)
est_DR.mx = NULL
est_IPQW.mx = NULL
est_IPSW.mx = NULL
est_IPSW2.mx = NULL
est_cfRF.mx = NULL
for(i in 1:ns){
  br = bootresult.cox[[i]]$t
  br.cf_RF = bootresult.ltrcrrf[[i]]$t
  est_DR.mx = rbind(est_DR.mx, br[,1:m])
  est_IPQW.mx = rbind(est_IPQW.mx, br[,m+(1:m)])
  est_IPSW.mx = rbind(est_IPSW.mx, br[,2*m+(1:m)])
  est_IPSW2.mx = rbind(est_IPSW2.mx, br[,3*m+(1:m)])
  est_cfRF.mx = rbind(est_cfRF.mx, br.cf_RF[,1:m])
}
est_DR.mx = est_DR.mx[1:nb.total,]
est_IPQW.mx = est_IPQW.mx[1:nb.total,]
est_IPSW.mx = est_IPSW.mx[1:nb.total,]
est_IPSW2.mx = est_IPSW2.mx[1:nb.total,]
est_cfRF.mx = est_cfRF.mx[1:nb.total,]

pbootCI_DR = apply(est_DR.mx, 2, quantile, probs = c(alpha/2, 1-alpha/2))
pbootCI_IPQW = apply(est_IPQW.mx, 2, quantile, probs = c(alpha/2, 1-alpha/2))
pbootCI_IPSW = apply(est_IPSW.mx, 2, quantile, probs = c(alpha/2, 1-alpha/2))
pbootCI_IPSW2 = apply(est_IPSW2.mx, 2, quantile, probs = c(alpha/2, 1-alpha/2))
pbootCI_cfRF = apply(est_cfRF.mx, 2, quantile, probs = c(alpha/2, 1-alpha/2))

pbootCI.l = cbind(cf_RF = pbootCI_cfRF[1,],
                  DR = pbootCI_DR[1,],
                  IPQW = pbootCI_IPQW[1,],
                  IPSW = pbootCI_IPSW[1,],
                  IPSW2 = pbootCI_IPSW2[1,], 
                  naive = NA)
pbootCI.u = cbind(cf_RF = pbootCI_cfRF[2,],
                  DR = pbootCI_DR[2,],
                  IPQW = pbootCI_IPQW[2,],
                  IPSW = pbootCI_IPSW[2,],
                  IPSW2 = pbootCI_IPSW2[2,], 
                  naive = NA)


### CW estimates in Vakulenko et. al. (2021) ----------------------------------------
source("CNS_data/func_C_before_T_severalZ.R") 
# Z={chemo, RT} 
X.s <- s$OS
T.s <- s$PFS
delta <- s$death
o <- order(X.s)
Z <- cbind(s$chemo.i[o], s$radiation[o])
X.s.o <- X.s[o]
T.s.o <- T.s[o]
d.o <- delta[o]

# estimation 
CW.est <- CW.deptr(X=X.s.o, T=T.s.o, delta=d.o, Z=Z, bs=TRUE, nbs.rep=500) # 
tvw <- surv.dep.trun(surv=X.s.o, trun=T.s.o, cens=d.o, zz=Z, bs=TRUE, nbs.rep=500)
# TVW alpha= 1.994624 0.7490019 

f.CW = stepfun(CW.est$time, c(1, CW.est$surv))
fCI_l.CW = stepfun(CW.est$time, c(1, CW.est$CI.L))
fCI_u.CW = stepfun(CW.est$time, c(1, CW.est$CI.U))




### plots -----------------------------------------------------------------------------------
names <- c("cf-RF-RF", "dr-Cox-Cox", "IPW.Q-Cox", "Reg.T1-Cox", "Reg.T2-Cox", "naive")  # naive = naive KM
xlim =  c(0,max(dat$X)+1)
ylim = c(-0.5,2.5)
color.pl = 8  #"grey"
lty.pl = 6
lwd.pl = 2
color = c(2,1,"purple",4,3,"brown")
lty = c(1,1,4,5,2,3)
lwd = c(2, 2, 2, 2, 2, 3)
color.cw = "orange" # 5, 
lty.cw = 4
lwd.cw = 2
cex = 2


## bound the CI's between 0 and 1
ylim = c(0,1)
pdf(paste("./CNS_data/c1.plots/c1CNS_surv_pbootCIb01_n",nrow(dat),"_cf_ntree_", ntree, ".pdf", sep = ""), 
    width = 5.5, height = 5.5)
# plot(plfit.T, conf.int = FALSE, col = color.pl,
#      lty = lty.pl, xlim = xlim, ylim = ylim, lwd = lwd.pl,
#      xlab = "Time (months)", ylab = "Survival probability")
plot(t0.list, est.pl, type = "l", col = color.pl,
     lty = lty.pl, xlim = xlim, ylim = ylim, lwd = lwd.pl,
     xlab = "Time (months)", ylab = "Survival probability")
for(j in 5:1){
  # shade the confidence intervals
  polygon(x = c(t0.list, rev(t0.list)),
          y = pmax(pmin(c(pbootCI.l[,j], rev(pbootCI.u[,j])),1),0),
          col =  adjustcolor(color[j], alpha.f = 0.2), border = NA)
}
#CW CI's
polygon(x = c(t0.list, rev(t0.list)),
        y = pmax(pmin(c(fCI_l.CW(t0.list), rev(fCI_u.CW(t0.list))),1),0),
        col =  adjustcolor(color.cw, alpha.f = 0.2), border = NA)
j=ncol(est)
# lines(kmfit.T, conf.int = FALSE, xlim = xlim,
#       col = color[j], lty = lty[j], lwd = lwd[j])
lines(t0.list, est.km, xlim = xlim, col = color[j], lty = lty[j], lwd = lwd[j])
lines(t0.list, f.CW(t0.list), col = color.cw, lty = lty.cw, lwd = lwd.cw)
for(j in 1:(ncol(est)-1)){
  lines(t0.list, est[,j], col = color[j], lty = lty[j], lwd = lwd[j])
}
legend("topright", legend = c(names[1:(ncol(est)-1)],"CW", "PL", names[ncol(est)]), 
       lty = c(lty[1:(ncol(est)-1)], lty.cw, lty.pl, lty[ncol(est)]), 
       col = c(color[1:(ncol(est)-1)], color.cw, color.pl, color[ncol(est)]), 
       lwd = c(lwd[1:(ncol(est)-1)], lwd.cw, lwd.pl, lwd[ncol(est)]),
       bty = "n", cex = 0.9)
dev.off()



## plot without the CI's
pdf(paste("./CNS_data/c1.plots/c1CNS_surv_n",nrow(dat),"_cf_ntree_", ntree, ".pdf", sep = ""), 
    width = 5, height = 5)
plot(t0.list, est.pl, type = "l", col = color.pl,
     lty = lty.pl, xlim = xlim, ylim = ylim, lwd = lwd.pl,
     xlab = "Time (months)", ylab = "Survival probability")
j=ncol(est)
lines(t0.list, est.km, xlim = xlim, col = color[j], lty = lty[j], lwd = lwd[j])
lines(t0.list, f.CW(t0.list), col = color.cw, lty = lty.cw, lwd = lwd.cw)
for(j in 1:(ncol(est)-1)){
  lines(t0.list, est[,j], col = color[j], lty = lty[j], lwd = lwd[j])
}
legend("topright", legend = c(names[1:(ncol(est)-1)],"CW", "PL", names[ncol(est)]), 
       lty = c(lty[1:(ncol(est)-1)], lty.cw, lty.pl, lty[ncol(est)]), 
       col = c(color[1:(ncol(est)-1)], color.cw, color.pl, color[ncol(est)]), 
       lwd = c(lwd[1:(ncol(est)-1)], lwd.cw, lwd.pl, lwd[ncol(est)]),
       bty = "n", cex = 0.9)
dev.off()










### estimates for the survival at certain time points, summarized in a table -----------------------------------
# estimates for the overall survival
t0.list = c(36,60,120,180)
cn = 8
nb.total = 2000
ns = cn
set.seed(123)
seeds.boot = sample(10000, size = ns)
nb = ceiling(nb.total/ns)

# dr-Cox-Cox, IPW.Q-Cox, Reg.T1-Cox, Reg.T2-Cox 
result.cox = c1.estSurv_DR.cox(t0.list, dat, covariates.X, covariates.Q, trim.C)

start_time <- Sys.time()
system.time({
  bootresult.cox = mclapply(seeds.boot, c1.splitboot_estSurv_DR.cox, mc.cores = getOption("mc.cores", cn),
                            n.boot = nb, dat = dat, t0 = t0.list,
                            covariates.X = covariates.X, covariates.Q = covariates.Q,
                            trim.C = trim.C)
})
now_time <- Sys.time()
print(paste("Time elapsed", now_time - start_time))
save(result.cox, bootresult.cox, t0.list, covariates.X, covariates.Q, trim.C, n.boot,
     file = paste("CNS_data/c1.results/CNS_table_n",nrow(dat),".sbootcox",nb.total, "_trimC", trim.C, ".rda", sep = ""))

# load(paste("CNS_data/c1.results/CNS_table_n",nrow(dat),".sbootcox",nb.total, "_trimC", trim.C, ".rda", sep = ""))

# cf-RF-RF
start_time <- Sys.time()
system.time({
  result.ltrcrrf = c1.estSurv_cf.ltrcrrf(t0.list, dat, covariates.X, covariates.Q, K,
                                 mtry, ntree, trim, trim.C)
})
now_time <- Sys.time()
print(paste("Time elapsed", now_time - start_time))
save(result.ltrcrrf, t0.list, covariates.X, covariates.Q, trim.C,
     file = paste("CNS_data/c1.results/CNS_table_n",nrow(dat),".ltrcrrf_ntree", ntree, ".rda", sep = ""))

start_time <- Sys.time()
system.time({
  bootresult.ltrcrrf = mclapply(seeds.boot, c1.splitboot_estSurv.ltrcrrf, mc.cores = getOption("mc.cores", cn),
                            n.boot = nb, dat = dat, t0 = t0.list,
                            covariates.X = covariates.X, covariates.Q = covariates.Q,
                            K = K, mtry = mtry, ntree = ntree, trim = trim, trim.C = trim.C)
})
now_time <- Sys.time()
print(paste("Time elapsed", now_time - start_time))
save(result.ltrcrrf, bootresult.ltrcrrf, t0.list, covariates.X, covariates.Q, K, mtry, ntree, trim, trim.C, n.boot,
     file = paste("CNS_data/c1.results/CNS_table_n",nrow(dat),".sbootRF",nb.total, "_trimC", trim.C, ".rda", sep = ""))

# load(paste("CNS_data/c1.results/CNS_table_n",nrow(dat),".sbootRF",nb.total, "_trimC", trim.C, ".rda", sep = ""))


### PL and naive KM estimator
surv.pl = stepfun(plfit.T$time, c(1, plfit.T$surv))
surv.km = stepfun(kmfit.T$time, c(1, kmfit.T$surv))

est_pl = surv.pl(t0.list)
est_km = surv.km(t0.list)

bootresult_pl = boot(dat, estSurv_pl.boot, R = nb.total, t0 = t0.list)
bootresult_km = boot(dat, estSurv_KM.boot, R = nb.total, t0 = t0.list)

pbootCI_pl = apply(bootresult_pl$t, 2, quantile, probs = c(alpha/2, 1-alpha/2))
pbootCI_km = apply(bootresult_km$t, 2, quantile, probs = c(alpha/2, 1-alpha/2))




# organize the results
m = length(t0.list)
est_DR.mx = NULL
est_IPQW.mx = NULL
est_IPSW.mx = NULL
est_IPSW2.mx = NULL
est_cfRF.mx = NULL
for(i in 1:ns){
  br = bootresult.cox[[i]]$t
  br.cf_RF = bootresult.ltrcrrf[[i]]$t
  est_DR.mx = rbind(est_DR.mx, br[,1:m])
  est_IPQW.mx = rbind(est_IPQW.mx, br[,m+(1:m)])
  est_IPSW.mx = rbind(est_IPSW.mx, br[,2*m+(1:m)])
  est_IPSW2.mx = rbind(est_IPSW2.mx, br[,3*m+(1:m)])
  est_cfRF.mx = rbind(est_cfRF.mx, br.cf_RF[,1:m])
}
est_DR.mx = est_DR.mx[1:nb.total,]
est_IPQW.mx = est_IPQW.mx[1:nb.total,]
est_IPSW.mx = est_IPSW.mx[1:nb.total,]
est_IPSW2.mx = est_IPSW2.mx[1:nb.total,]
est_cfRF.mx = est_cfRF.mx[1:nb.total,]

pbootCI_DR = apply(est_DR.mx, 2, quantile, probs = c(alpha/2, 1-alpha/2))
pbootCI_IPQW = apply(est_IPQW.mx, 2, quantile, probs = c(alpha/2, 1-alpha/2))
pbootCI_IPSW = apply(est_IPSW.mx, 2, quantile, probs = c(alpha/2, 1-alpha/2))
pbootCI_IPSW2 = apply(est_IPSW2.mx, 2, quantile, probs = c(alpha/2, 1-alpha/2))
pbootCI_cfRF = apply(est_cfRF.mx, 2, quantile, probs = c(alpha/2, 1-alpha/2))

pbootCI.l = rbind(cf_RF = pbootCI_cfRF[1,],
                  DR = pbootCI_DR[1,],
                  IPQW = pbootCI_IPQW[1,],
                  IPSW = pbootCI_IPSW[1,],
                  IPSW2 = pbootCI_IPSW2[1,],
                  CW = fCI_l.CW(t0.list),
                  PL = pbootCI_pl[1,],
                  KM = pbootCI_km[1,])
pbootCI.u = rbind(cf_RF = pbootCI_cfRF[2,],
                  DR = pbootCI_DR[2,],
                  IPQW = pbootCI_IPQW[2,],
                  IPSW = pbootCI_IPSW[2,],
                  IPSW2 = pbootCI_IPSW2[2,], 
                  CW = fCI_u.CW(t0.list),
                  PL = pbootCI_pl[2,],
                  KM = pbootCI_km[2,])
est = rbind(est_cf.RF = result.ltrcrrf$est_DR,
            est_DR = result.cox$est_DR,
            est_IPW.Q = result.cox$est_IPQW,
            est_Reg.T1 = result.cox$est_IPSW,
            est_Reg.T2 = result.cox$est_IPSW2,
            CW = f.CW(t0.list),
            PL = est_pl,
            KM = est_km)
method.names = c("cf-RF-RF","dr-Cox-Cox","IPW.Q", "Reg.T1","Reg.T2", "CW", "PL", "KM")
rownames(est) = method.names
colnames(est) = t0.list
rownames(pbootCI.l) = method.names
colnames(pbootCI.l) = t0.list
rownames(pbootCI.u) = method.names
colnames(pbootCI.u) = t0.list

tab = NULL
for(i in 1:length(t0.list)){
    tabi = sprintf("%.3f (%.3f, %.3f)", est[,i], pbootCI.l[,i], pmin(pbootCI.u[,i],1))
    tab = cbind(tab, tabi)
}
rownames(tab) = method.names
colnames(tab) = t0.list
tab


library(xtable)
xtable(tab,
       type = "latex", 
       digits = c(rep(3,5)),
       align = rep("c",5),
       file = "tab.tex")




### check the PH assumption for Cox models ------------------------------------------------
library(timereg)
dat$chemo = dat$chemo.i
set.seed(123)
fit.T = cox.aalen(Surv(Q, X, delta.1) ~prop(chemo)+prop(radiation), data = dat)
summary(fit.T)  # test of proportionality: prop(chemo.i) 0.964; prop(radiation): 0.334
pdf("./CNS_data/c1.plots/CNS_T_Mresidual.pdf",width = 8, height = 4)
par(mfrow=c(1,2))
plot(fit.T, score = 1, xlim = c(min(dat$X), max(dat$X)), xlab = "Time")
dev.off()


dat$Q2 = tau - dat$Q
dat$X2 = tau - dat$X
set.seed(123)
fit.Q2 = cox.aalen(Surv(X2, Q2, delta.1)~prop(chemo)+prop(radiation), data = dat)
summary(fit.Q2)  # test of proportionality: prop(chemo.i) 0.652; prop(radiation): 0.260
pdf("./CNS_data/c1.plots/CNS_Q_Mresidual.pdf",width = 8, height = 4)
par(mfrow=c(1,2))
plot(fit.Q2, score = 1, xlim = c(min(dat$Q2), max(dat$Q2)), xlab = "Reverse time")
dev.off()



