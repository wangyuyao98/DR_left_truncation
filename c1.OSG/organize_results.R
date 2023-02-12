# Organizing the results from OSG
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

itern = 500
n.boot = 100
nb.each = 20
ns = ceiling(n.boot/nb.each)

load("c1_Xcox.ltrcrrf/result.1_b1.rda")
bootse.mx = matrix(nrow = itern, ncol = ncol(bootresult.ltrcrrf$t))  # row - a data set;  col - an estimator
est.mx = matrix(nrow = itern, ncol = ncol(bootresult.ltrcrrf$t))  # row - a data set;  col - an estimator
colnames(bootse.mx) = names(bootresult.ltrcrrf$t0)
colnames(est.mx) = names(bootresult.ltrcrrf$t0)
for(i in 1:itern){
    load(paste("c1_Xcox.ltrcrrf/result.", i, "_b1.rda", sep = ""))
    est.mx[i,] = bootresult.ltrcrrf$t0
    bootse_i.mx = NULL
    for(j in 1:ns){
        load(paste("c1_Xcox.ltrcrrf/result.", i, "_b", j, ".rda", sep = ""))
        bootse_i.mx = rbind(bootse_i.mx, bootresult.ltrcrrf$t)
    }
    bootse_i.mx = bootse_i.mx[1:n.boot, ]
    bootse_i = apply(bootse_i.mx, 2, sd)
    bootse.mx[i,] = bootse_i
}

save(bootse.mx, file = paste("c1_Xcox_bootse_ltrcrrf_R", n.boot, ".rda", sep= ""))







