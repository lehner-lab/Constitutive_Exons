
library(ggplot2)
source("185_Functions.R")
Adamson <- read.table(file = "../Analyses/Data/Adamson_DeltaPSI.tsv",
                      header = T)
Adamson.Annotations <- read.delim(file = "147_Vex-seq_positions.tsv")

Adamson$WT.PSI.K <- Adamson$mean_PSI_K - Adamson$delta_psi_K
Adamson$WT.PSI.H <- Adamson$mean_PSI_H - Adamson$delta_psi_H
Adamson$WT.PSI.Mean <- apply(X = Adamson[,c("WT.PSI.K","WT.PSI.H")],
                             MARGIN = 1,
                             FUN = mean, na.rm = T)

pdf(file = "207_Adamson_DPSI_H_Distribution.pdf", width = 8, height = 7, useDingbats = F)
plot(x = density(Adamson$delta_psi_H, bw = 5, from = -100)$x,
     y = density(Adamson$delta_psi_H, bw = 5, from = -100)$y,
     type = "l",
     xlim = c(-100,100),
     ylim = c(0,0.05),
     xlab = "DPSI",
     ylab = "density",
     las = 1,
     lwd = 3,
     col = "gray70",
     main = "DPSI Distribution (Vex-seq Library in HepG2)")
dev.off()


pdf(file = "207_Adamson_DPSI_K_Distribution.pdf", width = 8, height = 7, useDingbats = F)
plot(x = density(Adamson$delta_psi_K, bw = 5, from = -100)$x,
     y = density(Adamson$delta_psi_K, bw = 5, from = -100)$y,
     type = "l",
     xlim = c(-100,100),
     ylim = c(0,0.05),
     xlab = "DPSI",
     ylab = "density",
     las = 1,
     lwd = 3,
     col = "gray70",
     main = "DPSI Distribution (Vex-seq Library in K562)")
dev.off()


