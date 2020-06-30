Adamson <- read.table(file = "../Analyses/Data/Adamson_DeltaPSI.tsv",
                      header = T)
Adamson.Annotations <- read.delim(file = "147_Vex-seq_positions.tsv")

Adamson$WT.PSI.K <- Adamson$mean_PSI_K - Adamson$delta_psi_K
Adamson$WT.PSI.H <- Adamson$mean_PSI_H - Adamson$delta_psi_H


Adamson$Annotation <- sapply(X = as.character(Adamson$variant),
                             FUN = function(x){
                               
                               annotation <- "none"
                               idx <- which(as.character(Adamson.Annotations$variant) == x)
                               
                               if (length(idx)>0){
                                 
                                 if (Adamson.Annotations$X3SS[idx] %in% c(-1,-2)) {
                                   annotation <- "acceptor"
                                 }
                                 
                                 if (Adamson.Annotations$X5SS[idx] %in% c(1,2)) {
                                   annotation <- "donor"
                                 }
                                 
                                 if (Adamson.Annotations$X3SS[idx] %in% c(1)) {
                                   annotation <- "exon.1"
                                 }
                                 
                                 if (Adamson.Annotations$X5SS[idx] %in% c(-1)) {
                                   annotation <- "exon.2"
                                 }
                                 
                                 annotation
                               }
                               
                               annotation
                             })

table(Adamson$Annotation)








SpliceSite.Mutations <- Adamson
SpliceSite.Mutations <- SpliceSite.Mutations[which(SpliceSite.Mutations$Annotation %in% c("acceptor", "donor")),]


HepG2.Dataset <- data.frame(Delta.PSI = SpliceSite.Mutations$delta_psi_H,
                            Starting.PSI = SpliceSite.Mutations$WT.PSI.H,
                            ID = SpliceSite.Mutations$variant)


Plot.DF <- data.frame(Starting.PSI = HepG2.Dataset$Starting.PSI,
                      Delta.PSI = HepG2.Dataset$Delta.PSI,
                      Final.PSI = HepG2.Dataset$Starting.PSI + HepG2.Dataset$Delta.PSI) 










Plot.DF$Starting.PSI.Group <- sapply(X = Plot.DF$Starting.PSI,
                                     FUN = function(x){
                                       group <- NA
                                       
                                       if (x <= 37.48543 ) {
                                         group <- 1
                                       } else if (x <= 66.18387 ) {
                                         group <- 2
                                       } else if (x <= 78.35756 ) {
                                         group <- 3
                                       } else if (x <= 84.62438 ) {
                                         group <- 4
                                       } else if (x <= 100) {
                                         group <- 5
                                       }
                                       
                                       group
                                       
                                     })

table(Plot.DF$Starting.PSI.Group, useNA = "always")
Starting.PSI.Means <- sapply(1:max(Plot.DF$Starting.PSI.Group),
                             function(x){
                               mean(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group == x)])
                             })

Plot.DF$Starting.PSI.Group <- factor(Plot.DF$Starting.PSI.Group, levels = 1:length(Starting.PSI.Means))




map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}


My.Colour.Palette <- colorRampPalette(c(colorRampPalette(c("dodgerblue1","white"))(n = 7)[c(1,3,5)],"white", colorRampPalette(c("firebrick1","white"))(n = 7)[c(5,3,1)]))(n = 999)

set.seed(0)
Violin.Plots <- ggplot(data = Plot.DF,
                       mapping = aes(x = Starting.PSI.Group,
                                     y = Final.PSI,
                                     fill = Starting.PSI.Group,
                                     colour = Starting.PSI.Group)) +
  geom_violin(scale = "width",trim = F, adjust = 4) +
  geom_jitter(width = 0.2, height = 0) +
  scale_fill_manual(values = map2color(x = Starting.PSI.Means[1:5], pal = My.Colour.Palette, limits = c(0,100))) +
  scale_colour_manual(values = rep("black", length(Starting.PSI.Means))) +
  theme_bw() +
  ylim(c(0,100)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.7,
        legend.position = "none") +
  xlab("starting PSI") +
  ylab("final PSI") +
  ggtitle("DPSI distributions at different starting PSIs")+
  scale_x_discrete(labels=c(as.character(round(Starting.PSI.Means, 4)))) +
  annotate(geom = "text",
           x = max(as.numeric(Plot.DF$Starting.PSI.Group)),
           y = 0,
           label = "n = 2-7")

ggsave(filename = "196_Violin_Plots_HepG2.pdf", plot = Violin.Plots, width = 8, height = 8, useDingbats = F)












































































































K562.Dataset <- data.frame(Delta.PSI = SpliceSite.Mutations$delta_psi_K,
                           Starting.PSI = SpliceSite.Mutations$WT.PSI.K,
                           ID = SpliceSite.Mutations$variant)


Plot.DF <- data.frame(Starting.PSI = K562.Dataset$Starting.PSI,
                      Delta.PSI = K562.Dataset$Delta.PSI,
                      Final.PSI = K562.Dataset$Starting.PSI + K562.Dataset$Delta.PSI)





Plot.DF$Starting.PSI.Group <- sapply(X = Plot.DF$Starting.PSI,
                                     FUN = function(x){
                                       group <- NA
                                       
                                       if (x <= 24.41629  ) {
                                         group <- 1
                                       } else if (x <= 57.81404 ) {
                                         group <- 2
                                       } else if (x <= 67.02732 ) {
                                         group <- 3
                                       } else if (x <= 81.20527 ) {
                                         group <- 4
                                       } else if (x <= 100) {
                                         group <- 5
                                       }
                                       
                                       group
                                       
                                     })

table(Plot.DF$Starting.PSI.Group, useNA = "always")
Starting.PSI.Means <- sapply(1:max(Plot.DF$Starting.PSI.Group),
                             function(x){
                               mean(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group == x)])
                             })

Plot.DF$Starting.PSI.Group <- factor(Plot.DF$Starting.PSI.Group, levels = 1:length(Starting.PSI.Means))
sapply(1:5, function(x){ max(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group==x)])})




map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}


My.Colour.Palette <- colorRampPalette(c(colorRampPalette(c("dodgerblue1","white"))(n = 7)[c(1,3,5)],"white", colorRampPalette(c("firebrick1","white"))(n = 7)[c(5,3,1)]))(n = 999)

set.seed(0)
Violin.Plots <- ggplot(data = Plot.DF,
                       mapping = aes(x = Starting.PSI.Group,
                                     y = Final.PSI,
                                     fill = Starting.PSI.Group,
                                     colour = Starting.PSI.Group)) +
  geom_violin(scale = "width",trim = F, adjust = 4) +
  geom_jitter(width = 0.2, height = 0) +
  scale_fill_manual(values = map2color(x = Starting.PSI.Means[1:5], pal = My.Colour.Palette, limits = c(0,100))) +
  scale_colour_manual(values = rep("black", length(Starting.PSI.Means))) +
  theme_bw() +
  ylim(c(0,100)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.7,
        legend.position = "none") +
  xlab("starting PSI") +
  ylab("final PSI") +
  ggtitle("DPSI distributions at different starting PSIs")+
  scale_x_discrete(labels=c(as.character(round(Starting.PSI.Means, 4)))) +
  annotate(geom = "text",
           x = max(as.numeric(Plot.DF$Starting.PSI.Group)),
           y = 0,
           label = "n = 1-7")

ggsave(filename = "196_Violin_Plots_K562.pdf", plot = Violin.Plots, width = 8, height = 8, useDingbats = F)



