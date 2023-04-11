### The main code for analyzing the HAAS data set with event time being time to cognitive impairment or death
### The censoring is always after left truncation (c2 type).

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(boot)
library(survival)
library(parallel)   # package for parallel computing
library(tranSurv)    # conditional Kendall's tau test
library(CoxR2)    # compute the R2 measures
library(LTRCforests)
library(cvTools)  # used for creating folds for cross-fitting

source("src/functions.cox.R")
source("src/c1.estimator.R")
source("src/c2.estimator.R")
source("src/c.simu_estimate.cox_naive_full.R")

dir.create("HAAS")
dir.create("HAAS/c2.results")
dir.create("HAAS/c2.plots")


### load and process data -----------------------------------------------------------
# load data
dat =  read.csv("~/Documents/research/R03 data/datasets/AgeX4_AgeDeathCog_n2560.csv")  # T = Age to death or cogImpair, Q = AgeX4
dat$Q = dat$AgeX4

tau = max(dat$X)+1
tau     # 102.6
dat$X2 = tau - dat$X
dat$Q2 = tau - dat$Q


# conditional Kendall's tau test
cKendall(dat$Q, dat$X, dat$delta)  # tau = 0.0426, p-val = 0.0014

# dichotomize alcohol consumption into heavy (1) and light (0) based on > 30.1
dat$AlcoholHL = as.numeric(dat$AlcoholX1 > 1.2*365.25/12)
# dat$smokeYN = as.numeric(dat$PackYearsX3 > 0)

covariates = c("Education", "E4Positive", "AlcoholHL", "PackYearsX3")
covariates.T = covariates
covariates.Q = covariates

t0.list = seq(min(dat$X[dat$delta==1])-1, max(dat$X)+1, by = 1)




# KM to estimate the censoring weights
# compute the survival curve Sd for the residual censoring time D:=C-Q from KM estimator
# If no censoring, Sd is the function that is constant 1
if(sum(1-dat$delta) == 0){
    Sd = function(x){return(1)}
}else{
    kmfit.D = survfit(Surv(X-Q, 1-delta)~1, data = dat, type = "kaplan-meier")
    Sd = stepfun(kmfit.D$time,  c(1, kmfit.D$surv) )  # Bound Sc away from 0 so we do not get extreme weights
}

# check the positivity assumption for Sc
summary(Sd(dat$X - dat$Q)[dat$delta==1])



### compute the estimators ------------------------------------------------------------
trim.C = 0
type.C = 'c2'  # The censoring type is `c2'

# PL estimator --------------------------------------------------------------------
plfit.T = survfit(Surv(Q, X, delta)~1, data = dat)
plsf.T = stepfun(plfit.T$time, c(1, plfit.T$surv))
est.pl = plsf.T(t0.list)


# KM estimator (Ignore left truncation)----------------------------------------------
kmfit.T = survfit(Surv(X, delta)~1, data = dat)
kmsf.T = stepfun(kmfit.T$time, c(1, kmfit.T$surv))
est.km = kmsf.T(t0.list)


# naive estimator (IPCW ignore lef truncation) --------------------------------------
est_naive = sapply(t0.list, c.estSurv_naive, dat = dat, type.C = type.C)



# dr-Cox-Cox, IPW.Q-Cox, Reg.T1-Cox, Reg.T2-Cox -------------------------------------
start_time <- Sys.time()
system.time({
    result.cox = c.estSurv_DR.cox(t0.list, dat, covariates.T, covariates.Q, trim.C, type.C)
    # bootresult.cox = boot(dat, c.estSurv_DR.cox, R = n.boot,
    #                       t0 = t0.list, covariates.T = covariates.T,
    #                       covariates.Q = covariates.Q, trim.C = trim.C, type.C = type.C)
})
now_time <- Sys.time()
print(paste("Time elapsed", now_time - start_time))
save(result.cox,  t0.list, covariates.T, covariates.Q, trim.C, type.C,
     file = paste("HAAS/c2.results/HAAS_n",nrow(dat),".cox_trim.C", trim.C, ".rda", sep = ""))

# load(paste("HAAS/c2.results/HAAS_n",nrow(dat),".cox_trim.C", trim.C, ".rda", sep = ""))


# Compute the bootstrap SE when using Cox to estimate F and G
n.boot = 100
cn = 6
start_time <- Sys.time()
system.time({
    bootresult.cox = mclapply(t0.list, c.estSurv_bootse_DR.cox, mc.cores = getOption("mc.cores", cn),
                      dat = dat, covariates.T = covariates.T, covariates.Q = covariates.Q,
                      trim.C = trim.C, type.C = type.C, n.boot = n.boot)
    save(bootresult.cox, t0.list, covariates.T, covariates.Q, trim.C, type.C,
         n.boot,
         file = paste("HAAS/c2.results/HAAS_n",nrow(dat),".coxboot", n.boot,
                      "_trim.C", trim.C, ".rda", sep = ""))
})
now_time <- Sys.time()
print(paste("Time elapsed", now_time - start_time))

# load(paste("HAAS/c2.results/HAAS_n",nrow(dat),".coxboot", n.boot, "_trim.C", trim.C, ".rda", sep = ""))

bootse.cox = matrix(nrow = length(t0.list), ncol = 4)
colnames(bootse.cox) = c("se.DR", "se.IPQW", "se.IPSW", "se.IPSW2")
for(i in 1:length(t0.list)){
    bresult = bootresult.cox[[i]]
    bootse.cox[i,] = apply(bresult$t, 2, sd)
}


est = cbind(cf_RF = NA, 
            DR = result.cox$est_DR, 
            IPQW = result.cox$est_IPQW, 
            IPSW = result.cox$est_IPSW, 
            IPSW2 = result.cox$est_IPSW2,
            naive = est_naive)

se = cbind(se.cf_RF = NA, 
           bootse.cox,
           se.naive = NA)

se_bootse = cbind(se.cf_RF = NA, 
                  bootse.cox,
                  se.naive = NA)

alpha = 0.05
qz = qnorm(alpha/2, lower = FALSE)



### draw plots --------------------------------------------------------------------
# setup
names <- c("cf-RF-RF", "dr-Cox-Cox", "IPW.Q-Cox", "Reg.T1-Cox", "Reg.T2-Cox", "naive")  # naive = naive KM
colnames(est) <- names
rownames(est) <- t0.list
colnames(se) <- names
rownames(se) <- t0.list

xlim =  c(70,max(dat$X)+1)
ylim = c(0,1)
color.pl = 8  #"grey"
lty.pl = 6
lwd.pl = 2
color = c(2,1,"purple",4,3,"brown")
lty = c(2,1,4,5,3,3)
lwd = c(2, 2, 2, 2, 3, 3)



# survival curve
pdf(paste("./HAAS/c2.plots/c2HAAS_surv_n",nrow(dat), ".pdf", sep = ""), 
    width = 5.5, height = 5.5)
plot(t0.list, est.pl, type = "l", col = color.pl,
     lty = lty.pl, xlim = xlim, ylim = ylim, lwd = lwd.pl,
     xlab = "Age (years)", ylab = "Survival probability")
for(j in 2:(ncol(est)-1)){
    lines(t0.list, est[,j], col = color[j], lty = lty[j], lwd = lwd[j])
}
j = ncol(est)
lines(t0.list, est.km, col = color[j], lty = lty[j], lwd = lwd[j])
legend("topright", legend = c(names[2:(ncol(est)-1)], "PL", names[ncol(est)]), 
       lty = c(lty[2:(ncol(est)-1)], lty.pl, lty[ncol(est)]), 
       col = c(color[2:(ncol(est)-1)], color.pl, color[ncol(est)]), 
       lwd = c(lwd[2:(ncol(est)-1)], lwd.pl, lwd[ncol(est)]),
       bty = "n")
dev.off()


# survival curve with boot CIs
pdf(paste("./HAAS/c2.plots/c2HAAS_survCIboot", n.boot, "_n",nrow(dat), ".pdf", sep = ""), 
    width = 5.5, height = 5.5)
plot(t0.list, est.pl, type = "l", col = color.pl,
     lty = lty.pl, xlim = xlim, ylim = ylim, lwd = lwd.pl,
     xlab = "Age (years)", ylab = "Survival probability")
for(j in 2:(ncol(est)-1)){
    # shade the confidence intervals
    polygon(x = c(t0.list, rev(t0.list)),
            y = c(est[,j] - qz*se_bootse[,j], rev(est[,j] + qz*se_bootse[,j])),  # adjust the CI's to be between 0 and 1
            col =  adjustcolor(color[j], alpha.f = 0.3), border = NA)  
}
for(j in 2:(ncol(est)-1)){
    lines(t0.list, est[,j], col = color[j], lty = lty[j], lwd = lwd[j])
}
j = ncol(est)
lines(t0.list, est.km, col = color[j], lty = lty[j], lwd = lwd[j])
legend("topright", legend = c(names[2:(ncol(est)-1)], "PL", names[ncol(est)]), 
       lty = c(lty[2:(ncol(est)-1)], lty.pl, lty[ncol(est)]), 
       col = c(color[2:(ncol(est)-1)], color.pl, color[ncol(est)]), 
       lwd = c(lwd[2:(ncol(est)-1)], lwd.pl, lwd[ncol(est)]),
       bty = "n")
dev.off()



# zoom-in plot
pdf(paste("./HAAS/c2.plots/c2HAAS_zoomin_survCIboot", n.boot, "_n",nrow(dat), ".pdf", sep = ""), 
    width = 5.5, height = 5.5)
xlim = c(84,96)
ylim = c(0,0.7)
plot(t0.list, est.pl, type = "l", col = color.pl,
     lty = lty.pl, xlim = xlim, ylim = ylim, lwd = lwd.pl,
     xlab = "Age (years)", ylab = "Survival probability")
# main = "(b)", font.main = font.main, cex.lab = cex.lab, cex.axis = cex.axis)
for(j in 2:(ncol(est)-1)){
    # shade the confidence intervals
    polygon(x = c(t0.list, rev(t0.list)),
            y = c(est[,j] - qz*se_bootse[,j], rev(est[,j] + qz*se_bootse[,j])),  # adjust the CI's to be between 0 and 1
            col =  adjustcolor(color[j], alpha.f = 0.2), border = NA)  
}
for(j in 2:(ncol(est)-1)){
    lines(t0.list, est[,j], col = color[j], lty = lty[j], lwd = lwd[j])
}
j = ncol(est)
lines(t0.list, est.km, col = color[j], lty = lty[j], lwd = lwd[j])
legend("topright", legend = c(names[2:(ncol(est)-1)], "PL", names[ncol(est)]), 
       lty = c(lty[2:(ncol(est)-1)], lty.pl, lty[ncol(est)]), 
       col = c(color[2:(ncol(est)-1)], color.pl, color[ncol(est)]), 
       lwd = c(lwd[2:(ncol(est)-1)], lwd.pl, lwd[ncol(est)]),
       bty = "n")
dev.off()














### estimates for the survival at certain time points -------------------------------------
# estimates for the overall survival
t0.list = c(80, 85, 90, 95)
# t0.list = 80:93
cn = 7
n.boot = 100

# dr-Cox-Cox, IPW.Q-Cox, Reg.T1-Cox, Reg.T2-Cox 
result.cox = c.estSurv_DR.cox(t0.list, dat, covariates.T, covariates.Q, trim.C, type.C)

start_time <- Sys.time()
system.time({
    bootresult.cox = mclapply(t0.list, c.estSurv_bootse_DR.cox, mc.cores = getOption("mc.cores", cn),
                              dat = dat, covariates.T = covariates.T, covariates.Q = covariates.Q,
                              trim.C = trim.C, type.C = type.C, n.boot = n.boot)
})
now_time <- Sys.time()
print(paste("Time elapsed", now_time - start_time))
save(result.cox, bootresult.cox, t0.list, covariates.T, covariates.Q, trim.C, type.C,
     n.boot,
     file = paste("HAAS/c2.results/HAAS_table_n",nrow(dat),".coxboot", n.boot,
                  "_trim.C", trim.C, ".rda", sep = ""))


bootse.cox = matrix(nrow = 4, ncol = length(t0.list))
rownames(bootse.cox) = c("bootse.DR", "bootse.IPQW", "bootse.IPSW", "bootse.IPSW2")
colnames(bootse.cox) = t0.list
for(i in 1:length(t0.list)){
    bresult = bootresult.cox[[i]]
    bootse.cox[,i] = apply(bresult$t, 2, sd)
}


### PL and naive KM estimator
surv.pl = stepfun(plfit.T$time, c(1, plfit.T$surv))
surv.km = stepfun(kmfit.T$time, c(1, kmfit.T$surv))

est_pl = surv.pl(t0.list)
est_km = surv.km(t0.list)

bootresult_pl = boot(dat, estSurv_pl.boot, R = n.boot, t0 = t0.list)
bootresult_km = boot(dat, estSurv_KM.boot, R = n.boot, t0 = t0.list)

bootse_pl = apply(bootresult_pl$t, 2, sd)
bootse_km = apply(bootresult_km$t, 2, sd)


# combining the result
est = rbind(DR = result.cox$est_DR, 
            IPQW = result.cox$est_IPQW, 
            IPSW = result.cox$est_IPSW, 
            IPSW2 = result.cox$est_IPSW2,
            PL = est_pl,
            KM = est_km)
bootse = rbind(bootse.cox, bootse_pl = bootse_pl, bootse_km = bootse_km)
bootCI.l = est - qz*bootse
bootCI.u = est + qz*bootse

method.names = c("dr-Cox-Cox", "IPW.Q-Cox", "Reg.T1-Cox", "Reg.T2-Cox", "PL", "KM")
tab = NULL
for(i in 1:length(t0.list)){
    tabi = sprintf("%.3f (%.3f, %.3f)", est[,i], bootCI.l[,i], bootCI.u[,i])
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
dat$delta.1 = rep(1, nrow(dat))

dat$Alcohol = dat$AlcoholHL
dat$Cigarette = dat$PackYearsX3

# compute the survival curve Sd for the residual censoring time D:=C-Q from KM estimator
# If no censoring, Sd is the function that is constant 1
if(sum(1-dat$delta) == 0){
    Sd = function(x){return(1)}
}else{
    kmfit.D = survfit(Surv(X-Q, 1-delta)~1, data = dat, type = "kaplan-meier")
    Sd = stepfun(kmfit.D$time,  pmax(c(1, kmfit.D$surv), trim.C) )  # Bound Sc away from 0 so we do not get extreme weights
}
# IPCW weights
wd = rep(0, nrow(dat))
id1 = which(dat$delta == 1)
wd[id1] = 1/Sd(dat$X-dat$Q)[id1]
dat$wd = wd
dat1 = dat[id1,]

dim(dat)
sum(dat$wd == 0)

set.seed(123)
fit.T = cox.aalen(Surv(Q, X, delta) ~prop(Education)+prop(E4Positive)+prop(Alcohol)+prop(Cigarette),
                  data = dat)
summary(fit.T)  # test of proportionality: Education: 0.002; E4positive: 0.418; Alcohol: 0.064; Cigarette: 0.006



pdf("./HAAS/c2.plots/HAAS_T_Mresidual.pdf",width = 8, height = 6)
par(mfrow=c(2,2))
plot(fit.T, score = 1, xlim = c(min(dat$X), max(dat$X)), xlab = "Time")
dev.off()

# manually add the p-values
par(mfrow=c(1,1))
set.seed(123)
plot(fit.T, score = 1, xlim = c(min(dat$X), max(dat$X)), xlab = "Time")

# size 3.5*4.5
# Cigarette
text(0.94 * par('usr')[2], -0.87 * par('usr')[4], labels = paste("p-value:", 0.006))
# HAAS_T_Mresid_Cigarette

# Alcohol
text(0.94 * par('usr')[2], -0.87 * par('usr')[4], labels = paste("p-value:", 0.064))
# HAAS_T_Mresid_Alcohol

# E4positive
text(0.94 * par('usr')[2], -0.87 * par('usr')[4], labels = paste("p-value:", 0.418))
# HAAS_T_Mresid_E4positive

# Education
text(0.94 * par('usr')[2], -0.87 * par('usr')[4], labels = paste("p-value:", sprintf("%.3f", 0.002)))
# HAAS_T_Mresid_Education



dat$Q2 = tau - dat$Q
dat$X2 = tau - dat$X
dat1 = dat[id1,]
set.seed(123)
fit.Q2 = cox.aalen(Surv(X2, Q2, delta.1) ~prop(Education)+prop(E4Positive)+prop(Alcohol)+prop(Cigarette),
                   data = dat1, caseweight = dat1$wd)
summary(fit.Q2)  # test of proportionality: Education: 0.002; E4positive: 0.062; Alcohol: 0.392; Cigarette: 0.492


pdf("./HAAS/c2.plots/HAAS_Q_Mresidual.pdf",width = 8, height = 6)
par(mfrow=c(2,2))
plot(fit.Q2, score = 1, xlim = c(min(dat$Q2), max(dat$Q2)), xlab = "Reverse time")
dev.off()

par(mfrow=c(1,1))
set.seed(123)
plot(fit.Q2, score = 1, xlim = c(min(dat$Q2), max(dat$Q2)), xlab = "Reverse time")

# size 3.5*4.5
# Cigarette
text(0.45 * par('usr')[2], -0.87 * par('usr')[4], labels = paste("p-value:", sprintf("%.3f", 0.492)))
# HAAS_Q2_Mresid_Cigarette

# Alcohol
text(0.45 * par('usr')[2], -0.87 * par('usr')[4], labels = paste("p-value:", 0.392))
# HAAS_Q2_Mresid_Alcohol

# E4positive
text(0.45 * par('usr')[2], -0.87 * par('usr')[4], labels = paste("p-value:", 0.062))
# HAAS_Q2_Mresid_E4positive

# Education
text(0.45 * par('usr')[2], -0.87 * par('usr')[4], labels = paste("p-value:", sprintf("%.3f", 0.002)))
# HAAS_Q2_Mresid_Education



