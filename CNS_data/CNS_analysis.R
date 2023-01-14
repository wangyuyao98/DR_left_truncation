# analysis of the CNS lymphoma data (Vakulenko-Lagun, Qian, Chiou, Wang, Betensky, 
# "Nonparametric estimation of the survival distribution under covariate-induced dependent truncation")

rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(survival)
library(ggplot2)
library(ggalt) # for stepribbon in ggplot
library(depend.truncation) # copula methods
library(tranSurv) 
source("func_C_before_T_severalZ.R") 


d <- read.table(file="CNS_data.txt", sep="\t", header=TRUE)
#d$site.i <- ifelse(d$site==1, 1, 0) 
d$chemo.i <- ifelse(d$chemo==1, 1, 0)

s <- d[d$relapse==1,] # select only those who progressed to relapse

obs10 <- min(s$OS[s$chemo.i==1])
s <- s[s$OS!=obs10,] # omitting an observation that zeroes out the risk set in the beginning of follow-up 

# Z={chemo, RT} -----------------------------------------------
X.s <- s$OS
T.s <- s$PFS
delta <- s$death
o <- order(X.s)
Z <- cbind(s$chemo.i[o], s$radiation[o])
X.s.o <- X.s[o]
T.s.o <- T.s[o]
d.o <- delta[o]


# -----testing quasi-independence ----------------------------
cond.Kt <- cKendall(T.s, X.s, delta)
cond.Kt$PE
cond.Kt$p.value
# > cond.Kt$PE
# [1] 0.01591343
# > cond.Kt$p.value
# [1] 0.8443915

# weighted logrank test for QI based on copula alternatives:
# for residual censoring:
Logrank.stat(T.s, X.s, delta)
Logrank.stat.tie(T.s, X.s, delta)

run.test <- function(T.s, X.s, delta) 
  # jackknife variance estimation of the test statistics, and the pvalues of the test
{
  st <- Logrank.stat.tie(T.s, X.s, delta)
  n <- length(T.s)
  jk <- matrix(NA, n, 3)
  for (i in 1:n)
    jk[i,] <- Logrank.stat.tie(T.s[-i], X.s[-i], delta[-i])
  avg <- apply(jk, 2, mean, na.rm=TRUE)
  M <- matrix(avg, n, 3, byrow=TRUE)
  Var <- n/(n-1)*colSums((jk-M)^2)
  S <- st/sqrt(Var)
  pv <- 2*pnorm(-abs(S))
  list(logrank.test= st, Zvalue=S, s2=Var, pvalue=pv)
}
run.test(T.s, X.s, delta) 
# $logrank.test
# [1]  3.4374766 -0.1338007  0.6913417
# $s2
# [1] 46.668201  4.673401 11.649753
# $pvalue
# [1] 0.614833 0.950648 0.839486


# estimation -----------------------------------
KM.true <- survfit(Surv(time=d$OS, event=d$death) ~1) # unconditional survival function of all Complete Respondents (CR)
KM.indep <- survfit(Surv(time=T.s.o, time2=X.s.o, event=d.o) ~1) # Km estimator that assumes independent truncation
CW.est <- CW.deptr(X=X.s.o, T=T.s.o, delta=d.o, Z=Z, bs=TRUE, nbs.rep=500) # 
tvw <- surv.dep.trun(surv=X.s.o, trun=T.s.o, cens=d.o, zz=Z, bs=TRUE, nbs.rep=500)
# TVW alpha= 1.994624 0.7490019 

minX <- min(CW.est$time) # 8.815789 
min.true.S <-  get.S(minX, summary(KM.true)$time, summary(KM.true)$surv) #P(X>minX)=0.9883379

# an approach that tries to account for residual dependence
Z <- cbind(s$chemo.i[o], s$radiation[o], X.s[o]) # including truncation time as one of the covariates
CW.est.res <- CW.deptr(X=X.s.o, T=T.s.o, delta=d.o, Z=Z, bs=TRUE, nbs.rep=500)
tvw.res <- surv.dep.trun(surv=X.s.o, trun=T.s.o, cens=d.o, zz=Z, bs=TRUE, nbs.rep=500)
# TVW alpha= 2.010543 0.6705764 -0.002611981 


# estimators based on copula methods:
res.frank <- EMURA.Frank(T.s.o, X.s.o, d.o,
                         plotX = FALSE, plotY = FALSE)
res.frank[1:3] 
pr.frank <-  get.S(minX, X.s.o, res.frank$Sy)

res.clayton <- EMURA.Clayton(T.s.o, X.s.o, d.o,
                             plotX = FALSE, plotY = FALSE)
res.clayton[1:3]
pr.clayton <-  get.S(minX, X.s.o, res.clayton$Sy)

# > res.frank[1:3]
# $alpha
# [1] 0.8861669
# $tau
# [1] 0.01342582
# $c
# [1] 0.4602074
# > res.clayton <- EMURA.Clayton(T.s.o, X.s.o, d.o,
#                        +                              plotX = FALSE, plotY = FALSE)
# > res.clayton[1:3]
# $alpha
# [1] 1.084743
# $tau
# [1] -0.04064931
# $c
# [1] 0.3957038


# plotting the results -----------------------------------------------------

res <- rbind( data.frame(time=summary(KM.true)$time[summary(KM.true)$time>=minX], 
                         surv=summary(KM.true)$surv[summary(KM.true)$time>=minX]/min.true.S, 
                         L.95=summary(KM.true)$surv[summary(KM.true)$time>=minX]/min.true.S,
                         U.95=summary(KM.true)$surv[summary(KM.true)$time>=minX]/min.true.S,
                         method=1),
              data.frame(time=c(summary(KM.indep)$time[1], summary(KM.indep)$time), 
                         surv=c(1, summary(KM.indep)$surv),
                         L.95=c(1,summary(KM.indep)$lower), 
                         U.95=c(1, summary(KM.indep)$upper),
                         method=2),
              data.frame(time=c(CW.est$time[1], CW.est$time), 
                         surv=c(1, CW.est$surv), 
                         L.95=c(1,CW.est$CI.L), U.95=c(1,CW.est$CI.U),
                         method=3),
              data.frame(time=c(tvw$times[1], tvw$times), 
                         surv=c(1,tvw$surv.est), 
                         L.95=c(1,tvw$CI.L.bs), U.95=c(1,tvw$CI.U.bs),
                         method=4)
)
lev <- c("responders Pr(X>x|X>8.8)", "KM.indep.trunc", "CW", "TVW") #, "Frank's copula",  "Clayton's copula", "CW-residual")
res$Method <- factor(res$method, ordered=TRUE)
levels(res$Method) <- lev

p <- ggplot(res, aes(x=time, y=surv, 
                     color=Method, 
                     fill=Method,
                     linetype=Method)) +
  geom_step(size=1.2) +
  geom_ribbon(aes(ymin=L.95, ymax=U.95), alpha=0.2, stat="stepribbon") +
  scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1")
p <- p + 
  xlab("x=time in months")+
  ylab("Pr(X > x | X > 8.8 months)") +
  theme(axis.text.x = element_text(face="bold", angle=0),
        axis.text.y = element_text(face="bold"))+
  theme(legend.position = c(0.79, 0.75),
        panel.background=element_rect(fill="transparent"),
        plot.background=element_rect(fill="transparent"),
        legend.background=element_rect(fill="transparent"),
        legend.box.background = element_rect(fill="transparent", color=NA),
        legend.key = element_rect(fill="transparent"))  
p
ggsave("CNS_lymphoma_chemo_RT_w_CI.png", bg="transparent", width=5, height=5)


res <- rbind( data.frame(time=summary(KM.true)$time[summary(KM.true)$time>=minX], 
                         surv=summary(KM.true)$surv[summary(KM.true)$time>=minX]/min.true.S, 
                         L.95=summary(KM.true)$surv[summary(KM.true)$time>=minX]/min.true.S,
                         U.95=summary(KM.true)$surv[summary(KM.true)$time>=minX]/min.true.S,
                         method=1),
              data.frame(time=c(summary(KM.indep)$time[1], summary(KM.indep)$time), 
                         surv=c(1, summary(KM.indep)$surv),
                         L.95=c(1,summary(KM.indep)$lower), 
                         U.95=c(1, summary(KM.indep)$upper),
                         method=2),
              data.frame(time=c(CW.est$time[1], CW.est$time), 
                         surv=c(1, CW.est$surv), 
                         L.95=c(1,CW.est$CI.L), U.95=c(1,CW.est$CI.U),
                         method=3),
              data.frame(time=c(tvw$times[1], tvw$times), 
                         surv=c(1,tvw$surv.est), 
                         L.95=c(1,tvw$CI.L.bs), U.95=c(1,tvw$CI.U.bs),
                         method=4),
              data.frame(time=X.s.o[X.s.o>=minX], 
                         surv=res.frank$Sy[X.s.o>=minX]/pr.frank, 
                         L.95=rep(NA, length(d.o)),
                         U.95=rep(NA, length(d.o)),
                         method=5),
              data.frame(time=X.s.o[X.s.o>=minX], 
                         surv=res.clayton$Sy[X.s.o>=minX]/pr.clayton, 
                         L.95=rep(NA, length(d.o)),
                         U.95=rep(NA, length(d.o)),
                         method=6),
              data.frame(time=c(CW.est.res$time[1], CW.est.res$time), 
                         surv=c(1, CW.est.res$surv), 
                         L.95=c(1,CW.est.res$CI.L), U.95=c(1,CW.est.res$CI.U),
                         method=7)
)
lev <- c("responders Pr(X>x|X>8.8)", "KM.indep.trunc", "CW", "TVW", "Frank's copula", 
         "Clayton's copula", "CW-residual")
res$Method <- factor(res$method, ordered=TRUE)
levels(res$Method) <- lev


# without CI: 
p <- ggplot(res, aes(x=time, y=surv, color=Method, 
                     linetype=Method, fill=Method)) +
  scale_color_brewer(palette = "Set1")+
  geom_step(size=1.2) 
p <- p + 
  xlab("x=time in months")+
  ylab("Pr(X > x | X > 8.8 months)") +
  theme(axis.text.x = element_text(face="bold", angle=0),
        axis.text.y = element_text(face="bold"))+
  theme(legend.position = c(0.79, 0.75),
        panel.background=element_rect(fill="transparent"),
        plot.background=element_rect(fill="transparent"),
        legend.background=element_rect(fill="transparent"),
        legend.box.background = element_rect(fill="transparent", color=NA),
        legend.key = element_rect(fill="transparent")) 
p
ggsave("CNS_lymphoma_chemo_RT_wo_CI.png", bg="transparent", width=5, height=5)



