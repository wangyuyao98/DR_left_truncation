## The main code for visualizing the simulation results
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


### load and reorganize the results
# results when both models are correct
load("results/Tweibull1_QCox1_coxZ1Z2_naive_full_R100.rda")
tab_Tcox1_Qcox1 = tab

# T Cox1, Q Cox2
load("results/Tweibull1_QCox1_coxTZ1Z2QZ1sqZ1Z2_naive_full_R100.rda")
tab_Tcox1_Qcox2 = tab

# T Cox2, Q Cox1
load("results/Tweibull1_QCox1_coxTZ1sqZ1Z2QZ1Z2_naive_full_R100.rda")
tab_Tcox2_Qcox1 = tab

# results when both models are wrong
load("results/Tweibull1_QCox1_coxZ1Z2sq_naive_full_R100.rda")
tab_Tcox2_Qcox2 = tab

# results for the PL estimate 
load("results/Tweibull1_QCox1_pl_R100.rda")
tab_PL = tab


# T RF, Q RF, ntree = 100, mtry = 2
load("OSG/ltrcrrf2_2/result.rda")
tab_TRF_QTF = tab.TRF_QRF

est_dr = rbind(tab_Tcox1_Qcox1["DR", ], 
               tab_Tcox1_Qcox2["DR", ], 
               tab_Tcox2_Qcox1["DR", ], 
               tab_Tcox2_Qcox2["DR", ],
               tab_TRF_QTF["DR", ])
rownames(est_dr) = c("dr-Cox1-Cox1","dr-Cox1-Cox2","dr-Cox2-Cox1","dr-Cox2-Cox2",
                     "cf-RF-RF")
est_dr

est_IPW.Q = rbind(tab_Tcox1_Qcox1["IPQW", ], 
                  tab_Tcox2_Qcox2["IPQW", ],
                  tab_TRF_QTF["IPQW", ])
rownames(est_IPW.Q) = c("IPW.Q-Cox1","IPW.Q-Cox2", "IPW.Q-RF")
est_IPW.Q

est_Reg.T1 = rbind(tab_Tcox1_Qcox1["IPSW", ], 
                   tab_Tcox2_Qcox2["IPSW", ], 
                   tab_TRF_QTF["IPSW", ])
rownames(est_Reg.T1) = c("Reg.T1-Cox1","Reg.T1-Cox2", "Ref.T1-RF")
est_Reg.T1

est_Reg.T2 = rbind(tab_Tcox1_Qcox1["IPSW2", ], 
                   tab_Tcox2_Qcox2["IPSW2", ],
                   tab_TRF_QTF["IPSW2", ])
rownames(est_Reg.T2) = c("Reg.T2-Cox1","Reg.T2-Cox2", "Ref.T2-RF")
est_Reg.T2

est_PL_naive_full = rbind(tab_PL, tab_Tcox1_Qcox1["naive", ], tab_Tcox1_Qcox1["full", ])
rownames(est_PL_naive_full) = c("PL","naive","full")
est_PL_naive_full

colnames(est_dr) = colnames(est_PL_naive_full) 
colnames(est_IPW.Q) = colnames(est_PL_naive_full) 
colnames(est_Reg.T1) = colnames(est_PL_naive_full) 
colnames(est_Reg.T2) = colnames(est_PL_naive_full) 


# save(est_dr, est_IPW.Q, est_Reg.T1, est_Reg.T2, est_PL_naive_full,
#      file = "simu_results_Tweibull1_Qcox1.rda")




### visualize the results
library(ggplot2)
tab = rbind(est_dr, est_IPW.Q, est_Reg.T1, est_Reg.T2, est_PL_naive_full)
bias = abs(tab[,"bias"])
SD = tab[,"sd"]
bootCP = tab[,"bootCP"]
xx = c(1:nrow(est_dr), 
       nrow(est_dr)+(1:nrow(est_IPW.Q))+1, 
       nrow(est_dr)+nrow(est_IPW.Q)+(1:nrow(est_Reg.T1))+2,
       nrow(est_dr)+nrow(est_IPW.Q)+nrow(est_Reg.T1)+(1:nrow(est_Reg.T2))+3,
       nrow(est_dr)+nrow(est_IPW.Q)+nrow(est_Reg.T1)+nrow(est_Reg.T2)+(1:nrow(est_PL_naive_full))+4)
methods = rownames(tab)
# m.colors = c(1,1,1,2,1,2,1,2,1,2,2,2,1)
ccolor = as.factor(c(rep(1,nrow(est_dr)), rep(2,nrow(est_IPW.Q)), 
                     rep(3,nrow(est_Reg.T1)), rep(4,nrow(est_Reg.T2)), rep(6,3)))
# ccolor = c(c(rep("red",nrow(est_dr)), rep("blue",nrow(est_IPW.Q)), 
#              rep("purple",nrow(est_Reg.T1)), rep("green",nrow(est_Reg.T2)), rep("black",3)))
data = data.frame(xx = xx, bias = bias, SD = SD, bootCP = bootCP)

data3 = rbind(data.frame(xx = xx, yy = bias, quantity = "Bias", yintercept = 0),
              data.frame(xx = xx, yy = SD, quantity = "SD", yintercept = 0),
              data.frame(xx = xx, yy = bootCP, quantity = "Coverage", yintercept = 0.95))



# One figure that contains the plots of the bias, SD, and bootCP
data3$ymax = rep(c(0.16,0.04,1.2), each = length(xx))
yy = data3$yy
yintercept = data3$yintercept
ccolor3 = rep(ccolor, 3)

data_hline <- data.frame(quantity = unique(data3$quantity),  # Create data for lines
                         hline = c(0, 0, 0.95))

pp = ggplot(data3, aes(x=xx, y=yy)) +
    geom_hline(data = data_hline, aes(yintercept = hline),  # Add different lines to facet plot
                 colour=8, lty = "dashed") +
    # geom_hline(yintercept = yintercept,colour="black", lty = "dashed") +
    geom_segment(aes(x = xx, xend = xx, y = 0, yend = yy), lwd = 1.5, color = ccolor3) +
    geom_point(size = 3, color = ccolor3) +
    labs(x="", y = "") + 
    theme_light() +
    theme(legend.position="none", axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
    scale_x_continuous(labels = methods, breaks = xx, limits = c(1,max(xx))) +
    #scale_y_continuous(labels = seq(0,0.15, by = 0.02), breaks = seq(0,0.15, by = 0.02), limits = c(0,0.15)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1)) +
    geom_text(aes(size=1), label = sprintf("%.3f", yy), vjust = -1, size = 3)+ 
    facet_grid(quantity ~ ., scales = "free_y")+
    geom_blank(aes(y = ymax))

pdf("plots/simu_bias_CP_SD.pdf",width = 9, height = 4.5)
pp
dev.off()





# bias
pdf("plots/simu_bias.pdf",width = 9, height = 3)
ggplot(data, aes(x=xx, y=bias)) +
    geom_hline(yintercept = 0,colour="black", lty = "dashed") +
    geom_segment(aes(x = xx, xend = xx, y = 0, yend = bias), lwd = 1.5, color = ccolor) +
    geom_point(size = 3, color = ccolor) +
    labs(x="", y = "Bias") + 
    theme_light() +
    theme(legend.position="none", axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
    scale_x_continuous(labels = methods, breaks = xx, limits = c(0,max(xx)+1)) +
    scale_y_continuous(labels = seq(0,0.15, by = 0.02), breaks = seq(0,0.15, by = 0.02), limits = c(0,0.15)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1)) +
    geom_text(aes(size=1), label = sprintf("%.3f", bias), vjust = -1, size = 3)
dev.off()

# SD
pdf("plots/simu_SD.pdf",width = 9, height = 3)
ggplot(data, aes(x=xx, y=SD)) +
    geom_hline(yintercept = 0, colour="black", lty = "dashed") +
    geom_segment(aes(x = xx, xend = xx, y = 0, yend = SD), lwd = 1.5, color = ccolor) +
    geom_point(size = 3, color = ccolor) +
    labs(x="", y = "SD") + 
    theme_light() +
    theme(legend.position="none", axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
    scale_x_continuous(labels = methods, breaks = xx, limits = c(0,max(xx)+1)) +
    scale_y_continuous(labels = seq(0,0.05,by=0.01), breaks = seq(0,0.05,by=0.01), limits = c(0,0.05)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1)) +
    geom_text(aes(size=1), label = sprintf("%.3f", SD), vjust = -1, size = 3)
dev.off()

# bootCP
pdf("plots/simu_bootCP.pdf",width = 9, height = 3)
ggplot(data, aes(x=xx, y=bootCP)) +
    geom_hline(yintercept = 0.95, colour="black", lty = "dashed") +
    geom_segment(aes(x = xx, xend = xx, y = 0, yend = bootCP), lwd = 1.5, color = ccolor) +
    geom_point(size = 3, color = ccolor) +
    labs(x="", y = "Coverage") + 
    theme_light() +
    theme(legend.position="none", axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
    scale_x_continuous(labels = methods, breaks = xx, limits = c(0,max(xx)+1)) +
    scale_y_continuous(labels = seq(0,1, by = 0.1), breaks = seq(0,1, by = 0.1), limits = c(0,1.1)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1)) +
    geom_text(aes(size=1), label = sprintf("%.3f", bootCP), vjust = -1, size = 3)
dev.off()



