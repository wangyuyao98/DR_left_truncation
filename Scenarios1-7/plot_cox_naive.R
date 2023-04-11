# plot plots for bias, bootstrap SE, bootstrap CP for naive estimator and 
# estimators with Cox model to fit F and G

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

methods = c("dr","IPW.Q","Reg.T2","naive")

scenarios = c("Tweibull1_Qweibull1", "Tweibull1_Qweibull2",
              "Tweibull1_Qweibull2_AFT2", "Tweibull2_Qweibull1",
              "TAFT2_weibull2_Qweibull1", "Tweibull2_Qweibull2",
              "TAFT2_weibull2_Qweibull2_AFT2")

title = c("1: Both correct", "2: Q model moderately wrong", "3: Q model severely wrong",
          "4: T model moderately wrong", "5: T model severely wrong",
          "6: Both moderately wrong", "7: Both severely wrong")


# plot a 3*7 plot for (bias, boot se, boot CP) * senarios 1-7
xx = 1:length(methods)
id = c(1,2,4,5)    # plot DR, IPQW, IPW.OR2, and naive estimator
pdf("bias_sd_bootcp.pdf",width = 5.4*length(scenarios), height = 5*3)
par(mfrow = c(3,length(scenarios)), cex.lab = 3.5, cex.axis = 2.7, cex.main = 3,
    mai = c(0.3,0.85,0.6,0.2))
# bias
ii = 1
load(paste("results/", scenarios[ii], "_cox_naive_full_R100.rda", sep = ""))
bias = tab[,"bias"]
# plot DR, IPQW, IPW.OR2 and naive estimator
plot(xx, abs(bias)[id], type = "h", lwd = 7, col = 1:5, xaxt = "n",
     xlab = "", ylab = 'Bias', xlim = c(0.8,4.2), ylim = c(0,0.32), 
     main = title[ii])
abline(h = 0, lty = "dashed", col = 1, lwd = 3)
points(abs(bias)[id], pch = 19, col = 1:5, cex = 3)
axis(1, at = xx, labels = methods)
text(xx, abs(bias)[id]+0.025, labels = sprintf("%.3f", abs(bias)[id]), col = 1:5,
     cex = 2.5)

for(ii in 2:length(scenarios)){
    load(paste("results/", scenarios[ii], "_cox_naive_full_R100.rda", sep = ""))
    bias = tab[,"bias"]
    plot(xx, abs(bias)[id], type = "h", lwd = 7, col = 1:5, xaxt = "n",
         xlab = "", ylab = "", xlim = c(0.8,4.2), ylim = c(0,0.32), 
         main = title[ii])
    abline(h = 0, lty = "dashed", col = 1, lwd = 3)
    points(abs(bias)[id], pch = 19, col = 1:5, cex = 3)
    axis(1, at = xx, labels = methods)
    text(xx, abs(bias)[id]+0.025, labels = sprintf("%.3f", abs(bias)[id]), col = 1:5,
         cex = 2.5)
}

# sd
ii = 1
load(paste("results/", scenarios[ii], "_cox_naive_full_R100.rda", sep = ""))
sd = tab[,"sd"]
# plot DR, IPQW, IPW.OR2 and naive estimator
plot(xx, sd[id], type = "h", lwd = 7, col = 1:4, xaxt = "n",
     xlab = "", ylab = 'SD', xlim = c(0.8,4.2), ylim = c(0,0.1))
points(sd[id], pch = 19, col = 1:4, cex = 3)
axis(1, at = xx, labels = methods)
text(xx, sd[id]+0.01, labels = sprintf("%.3f", sd[id]), col = 1:4,
     cex = 2.5)

for(ii in 2:length(scenarios)){
    load(paste("results/", scenarios[ii], "_cox_naive_full_R100.rda", sep = ""))
    sd = tab[,"sd"]
    # plot DR, IPQW, IPW.OR2 and naive estimator
    plot(xx, sd[id], type = "h", lwd = 7, col = 1:4, xaxt = "n",
         xlab = "", ylab = "", xlim = c(0.8,4.2), ylim = c(0,0.1))
    points(sd[id], pch = 19, col = 1:4, cex = 3)
    axis(1, at = xx, labels = methods)
    text(xx, sd[id]+0.01, labels = sprintf("%.3f", sd[id]), col = 1:4,
         cex = 2.5)
}

# boot cp
ii = 1
load(paste("results/", scenarios[ii], "_cox_naive_full_R100.rda", sep = ""))
cp_boot = tab[,"cover.boot"]
# plot DR, IPQW, IPW.OR2 and naive estimator
plot(xx, cp_boot[id], type = "h", lwd = 7, col = 1:4, xaxt = "n",
     xlab = "", ylab = 'Coverage', xlim = c(0.8,4.2), ylim = c(0,1.1))
points(cp_boot[id], pch = 19, col = 1:4, cex = 3)
abline(h = 0.95, lty = "dashed", col = 1, lwd = 3)
axis(1, at = xx, labels = methods)
text(xx, cp_boot[id]+0.1, labels = sprintf("%.3f", cp_boot[id]), col = 1:4,
     cex = 2.5)

for(ii in 2:length(scenarios)){
    load(paste("results/", scenarios[ii], "_cox_naive_full_R100.rda", sep = ""))
    cp_boot = tab[,"cover.boot"]
    # plot DR, IPQW, IPW.OR2 and naive estimator
    plot(xx, cp_boot[id], type = "h", lwd = 7, col = 1:4, xaxt = "n",
         xlab = "", ylab = "", xlim = c(0.8,4.2), ylim = c(0,1.1))
    abline(h = 0.95, lty = "dashed", col = 1, lwd = 3)
    points(cp_boot[id], pch = 19, col = 1:4, cex = 3,)
    axis(1, at = xx, labels = methods)
    text(xx, cp_boot[id]+0.1, labels = sprintf("%.3f", cp_boot[id]), col = 1:4,
         cex = 2.5)
}
dev.off()







