
library(ggplot2)
source("185_Functions.R")

Adamson <- read.table(file = "../Analyses/Data/Adamson_DeltaPSI.tsv",
                      header = T)
Adamson.Annotations <- read.delim(file = "../Analyses/Data/Adamson_VariantAnnotation.tsv")

Adamson$WT.PSI.K <- Adamson$mean_PSI_K - Adamson$delta_psi_K
Adamson$WT.PSI.H <- Adamson$mean_PSI_H - Adamson$delta_psi_H


Adamson$Annotation <- sapply(X = as.character(Adamson$variant),
                             FUN = function(x){
                               
                               annotation <- "none"
                               idx <- which(as.character(Adamson.Annotations$variant) == x)
                               
                               if (length(idx)>0){
                                 annotation <- as.character(Adamson.Annotations$annotation)[idx]
                               }
                               
                               annotation
                             })


rm(Adamson.Annotations)










Exonic.Class <- c("frameshift",
                  "inframe deletion",
                  "missense",
                  "synonymous",
                  "stop gained",
                  "stop lost",
                  "non coding transcript exon")


Exonic.Mutations <- Adamson
Exonic.Mutations <- Exonic.Mutations[which(Exonic.Mutations$Annotation %in% Exonic.Class),]


HepG2.Dataset <- data.frame(Delta.PSI = Exonic.Mutations$delta_psi_H,
                            Starting.PSI = Exonic.Mutations$WT.PSI.H,
                            ID = Exonic.Mutations$variant)



Plot.DF <- data.frame(Starting.PSI = HepG2.Dataset$Starting.PSI,
                      Delta.PSI = HepG2.Dataset$Delta.PSI,
                      Final.PSI = HepG2.Dataset$Starting.PSI + HepG2.Dataset$Delta.PSI) 










Plot.DF$Starting.PSI.Group <- EquallyPopulatedBins(vector = Plot.DF$Starting.PSI, k = 5)
table(Plot.DF$Starting.PSI.Group, useNA = "always")
Starting.PSI.Means <- sapply(1:max(Plot.DF$Starting.PSI.Group),
                             function(x){
                               mean(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group == x)])
                             })

Plot.DF$Starting.PSI.Group <- factor(Plot.DF$Starting.PSI.Group, levels = 1:length(Starting.PSI.Means))
sapply(1:5, function(x){ max(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group==x)])})

mean(abs(Plot.DF$Delta.PSI[which(Plot.DF$Starting.PSI.Group == 5)]) > 10)
median(abs(Plot.DF$Delta.PSI[which(Plot.DF$Starting.PSI.Group == 5)]))

mean(abs(Plot.DF$Delta.PSI[which(Plot.DF$Starting.PSI.Group == 2)]) > 10)
median(abs(Plot.DF$Delta.PSI[which(Plot.DF$Starting.PSI.Group == 2)]))



map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

My.Colour.Palette <- colorRampPalette(c(colorRampPalette(c("dodgerblue1","white"))(n = 7)[c(1,3,5)],"white", colorRampPalette(c("firebrick1","white"))(n = 7)[c(5,3,1)]))(n = 999)



Violin.Plots <- ggplot(data = Plot.DF,
                       mapping = aes(x = Starting.PSI.Group,
                                     y = Final.PSI,
                                     fill = Starting.PSI.Group,
                                     colour = Starting.PSI.Group)) +
  geom_violin(scale = "width",trim = F, adjust = 4) +
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
  scale_x_discrete(labels=c(as.character(round(Starting.PSI.Means, 2)))) +
  annotate(geom = "text",
           x = max(as.numeric(Plot.DF$Starting.PSI.Group)),
           y = 0,
           label = "n = 203-204")

ggsave(filename = "186_Violin_Plots_HepG2.pdf", plot = Violin.Plots, width = 8, height = 8, useDingbats = F)












































































































K562.Dataset <- data.frame(Delta.PSI = Exonic.Mutations$delta_psi_K,
                           Starting.PSI = Exonic.Mutations$WT.PSI.K,
                           ID = Exonic.Mutations$variant)


Plot.DF <- data.frame(Starting.PSI = K562.Dataset$Starting.PSI,
                      Delta.PSI = K562.Dataset$Delta.PSI,
                      Final.PSI = K562.Dataset$Starting.PSI + K562.Dataset$Delta.PSI)




Plot.DF$Starting.PSI.Group <- EquallyPopulatedBins(vector = Plot.DF$Starting.PSI, k = 5)
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


Violin.Plots <- ggplot(data = Plot.DF,
                       mapping = aes(x = Starting.PSI.Group,
                                     y = Final.PSI,
                                     fill = Starting.PSI.Group,
                                     colour = Starting.PSI.Group)) +
  geom_violin(scale = "width",trim = F, adjust = 4) +
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
  scale_x_discrete(labels=c(as.character(round(Starting.PSI.Means, 2)))) +
  annotate(geom = "text",
           x = max(as.numeric(Plot.DF$Starting.PSI.Group)),
           y = 0,
           label = "n = 203-204")

ggsave(filename = "186_Violin_Plots_K562.pdf", plot = Violin.Plots, width = 8, height = 8, useDingbats = F)






Legend.Colours <- map2color(x = seq(0,100,0.1),
                            pal = My.Colour.Palette,
                            limits = c(0,100))
plot(NULL,
     ylim = c(0,100),
     xlim = c(0,1))
segments(x0 = 0,
         x1 = 1,
         y0 = seq(0,100,0.1),
         y1 = seq(0,100,0.1),
         col = Legend.Colours)
