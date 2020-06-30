
# SRE Dataset
Kosuri <- read.table(file = "../Analyses/Data/sre_data_clean.txt", header = T)
Kosuri <- Kosuri[complete.cases(Kosuri),]
Kosuri$mean_nat_index <- apply(X = Kosuri[,c("nat_index_smn1", "nat_index_dhfr")],
                               MARGIN = 1,
                               FUN = mean, na.rm = T)

mean(Kosuri$mean_nat_index >=0.9)



hist(Kosuri$dpsi_dhfr)
hist(Kosuri$dpsi_smn1)

pdf(file = "208_SRE_DPSI_DHFR_Distribution.pdf", width = 8, height = 7, useDingbats = F)
plot(x = density(Kosuri$dpsi_dhfr[-which(is.na(Kosuri$dpsi_dhfr))]*100, bw = 5, from = -100, to = 100)$x,
     y = density(Kosuri$dpsi_dhfr[-which(is.na(Kosuri$dpsi_dhfr))]*100, bw = 5, from = -100, to = 100)$y,
     type = "l",
     xlim = c(-100,100),
     ylim = c(0,0.06),
     xlab = "DPSI",
     ylab = "density",
     las = 1,
     lwd = 3,
     col = "gray70",
     main = "DPSI Distribution (SRE Library in DHFR Construct)")
dev.off()


pdf(file = "208_SRE_DPSI_SMN1_Distribution.pdf", width = 8, height = 7, useDingbats = F)
plot(x = density(Kosuri$dpsi_smn1[-which(is.na(Kosuri$dpsi_smn1))]*100, bw = 5, from = -100, to = 100)$x,
     y = density(Kosuri$dpsi_smn1[-which(is.na(Kosuri$dpsi_smn1))]*100, bw = 5, from = -100, to = 100)$y,
     type = "l",
     xlim = c(-100,100),
     ylim = c(0,0.05),
     xlab = "DPSI",
     ylab = "density",
     las = 1,
     lwd = 3,
     col = "gray70",
     main = "DPSI Distribution (SRE Library in SMN1 Construct)")
dev.off()






Kosuri <- read.table(file = "../Analyses/Data/snv_data_clean.txt", header = T)

Kosuri$mean_dpsi <- apply(X = Kosuri[,c("v1_dpsi", "v2_dpsi")],
                          MARGIN = 1,
                          FUN = mean, na.rm = T)

Kosuri$mean_nat_index <- apply(X = Kosuri[,c("nat_v1_index","nat_v2_index")],
                               MARGIN = 1,
                               FUN = mean, na.rm = T)


mean(Kosuri$mean_nat_index >=0.9)


pdf(file = "208_SNV_DPSI_Distribution.pdf", width = 8, height = 7, useDingbats = F)
plot(x = density(Kosuri$mean_dpsi[-which(is.na(Kosuri$mean_dpsi))]*100, bw = 5, from = -100)$x,
     y = density(Kosuri$mean_dpsi[-which(is.na(Kosuri$mean_dpsi))]*100, bw = 5, from = -100)$y,
     type = "l",
     xlim = c(-100,100),
     ylim = c(0,0.06),
     xlab = "DPSI",
     ylab = "density",
     las = 1,
     lwd = 3,
     col = "gray70",
     main = "DPSI Distribution (SNV Library)")
dev.off()










