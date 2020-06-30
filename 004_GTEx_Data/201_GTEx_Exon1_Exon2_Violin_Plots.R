# all of GTEx



library(ggplot2)
source("185_Functions.R")

Tissues <- c("Adipose Tissue",
             "Adrenal Gland",
             "Bladder",
             "Blood",
             "Blood Vessel",
             "Bone Marrow",
             "Brain",
             "Breast",
             "Cervix Uteri",
             "Colon",
             "Esophagus",
             "Fallopian Tube",
             "Heart",
             "Kidney",
             "Liver",
             "Lung",
             "Muscle",
             "Nerve",
             "Ovary",
             "Pancreas",
             "Pituitary",
             "Prostate",
             "Salivary Gland",
             "Skin",
             "Small Intestine",
             "Spleen",
             "Stomach",
             "Testis",
             "Thyroid",
             "Uterus",
             "Vagina")


All.Starting.PSIs <- c()
All.Mutation.Effects <- c()


for (each.tissue in Tissues) {
  print(each.tissue)
  
  tissue.no.spaces <- gsub(pattern = " ",
                           replacement = "",
                           x = each.tissue,
                           fixed = TRUE)
  
  file.location <- paste("./187_GTEx_Exon_Mutation_Effects/",
                         tissue.no.spaces,
                         ".RData",
                         sep = "")
  
  load(file.location)
  
  if (length(Starting.PSIs)==0){
    next
  }
  
  All.Starting.PSIs <- c(All.Starting.PSIs, Starting.PSIs)
  All.Mutation.Effects <- c(All.Mutation.Effects, Mutation.Effects)
  
}

All.Mutation.Effects <- All.Mutation.Effects*2


# negative mutations
Plot.DF <- data.frame(Starting.PSI = All.Starting.PSIs,
                      Delta.PSI = All.Mutation.Effects,
                      Final.PSI = All.Starting.PSIs + All.Mutation.Effects)
Plot.DF <- Plot.DF[-which(Plot.DF$Starting.PSI < 0),]
Plot.DF <- Plot.DF[-which(Plot.DF$Starting.PSI > 1),]




Plot.DF$Starting.PSI.Group <- EquallyPopulatedBins(vector = Plot.DF$Starting.PSI, k = 15)
table(Plot.DF$Starting.PSI.Group, useNA = "always")
Starting.PSI.Means <- sapply(1:max(Plot.DF$Starting.PSI.Group),
                             function(x){
                               mean(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group == x)])
                             })

Plot.DF$Starting.PSI.Group <- factor(Plot.DF$Starting.PSI.Group, levels = 1:length(Starting.PSI.Means))
sapply(1:15, function(x){ max(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group==x)])})

mean(abs(Plot.DF$Delta.PSI)[which(Plot.DF$Starting.PSI.Group==15)]<0.1)*100


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
  scale_fill_manual(values = map2color(x = Starting.PSI.Means[1:25], pal = My.Colour.Palette, limits = c(0,1))) +
  scale_colour_manual(values = rep("black", length(Starting.PSI.Means))) +
  theme_bw() +
  ylim(c(0,1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.3,
        legend.position = "none") +
  xlab("starting PSI") +
  ylab("final PSI") +
  ggtitle("DPSI distributions at different starting PSIs")+
  scale_x_discrete(labels=c(as.character(round(Starting.PSI.Means, 4)*100))) +
  annotate(geom = "text",
           x = max(as.numeric(Plot.DF$Starting.PSI.Group)),
           y = 0,
           label = "n = 7975-7976")

ggsave(filename = "201_Violin_Plots_GTEx_15_Groups_Exonwide.pdf", plot = Violin.Plots, width = 12, height = 12, useDingbats = F)












































































































library(ggplot2)
Tissues <- c("Adipose Tissue",
             "Adrenal Gland",
             "Bladder",
             "Blood",
             "Blood Vessel",
             "Bone Marrow",
             "Brain",
             "Breast",
             "Cervix Uteri",
             "Colon",
             "Esophagus",
             "Fallopian Tube",
             "Heart",
             "Kidney",
             "Liver",
             "Lung",
             "Muscle",
             "Nerve",
             "Ovary",
             "Pancreas",
             "Pituitary",
             "Prostate",
             "Salivary Gland",
             "Skin",
             "Small Intestine",
             "Spleen",
             "Stomach",
             "Testis",
             "Thyroid",
             "Uterus",
             "Vagina")


All.Starting.PSIs <- c()
All.Mutation.Effects <- c()


for (each.tissue in Tissues) {
  print(each.tissue)
  
  tissue.no.spaces <- gsub(pattern = " ",
                           replacement = "",
                           x = each.tissue,
                           fixed = TRUE)
  
  file.location <- paste("./156_Exon1_Effects_vs_Starting_PSI/",
                         tissue.no.spaces,
                         ".RData",
                         sep = "")
  
  load(file.location)
  
  if (length(Starting.PSIs)==0){
    next
  }
  
  All.Starting.PSIs <- c(All.Starting.PSIs, Starting.PSIs)
  All.Mutation.Effects <- c(All.Mutation.Effects, Mutation.Effects)
  
}

All.Mutation.Effects <- All.Mutation.Effects*2



# negative mutations
Plot.DF <- data.frame(Starting.PSI = All.Starting.PSIs,
                      Delta.PSI = All.Mutation.Effects,
                      Final.PSI = All.Starting.PSIs + All.Mutation.Effects)
Plot.DF <- Plot.DF[-which(Plot.DF$Starting.PSI < 0),]
Plot.DF <- Plot.DF[-which(Plot.DF$Starting.PSI > 1),]

# Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI < 0)] <- 0
# Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI > 1)] <- 1




Plot.DF$Starting.PSI.Group <- sapply(X = Plot.DF$Starting.PSI,
                                     FUN = function(x){
                                       group <- NA
                                       
                                       if (x <= 0.001594883  ) {
                                         group <- 1
                                       } else if (x <= 0.004624439) {
                                         group <- 2
                                       } else if (x <= 0.010841853 ) {
                                         group <- 3
                                       } else if (x <= 0.023742260 ) {
                                         group <- 4
                                       } else if (x <= 0.049948917) {
                                         group <- 5
                                       } else if (x <= 0.112350678) {
                                         group <- 6
                                       } else if (x <= 0.295190228) {
                                         group <- 7
                                       } else if (x <= 0.657180986) {
                                         group <- 8
                                       } else if (x <= 0.867587457) {
                                         group <- 9
                                       } else if (x <= 0.940297361) {
                                         group <- 10
                                       } else if (x <= 0.969923057) {
                                         group <- 11
                                       } else if (x <= 0.985507246) {
                                         group <- 12
                                       } else if (x <= 0.993644307) {
                                         group <- 13
                                       } else if (x <= 0.997765643) {
                                         group <- 14
                                       } else if (x <= 1) {
                                         group <- 15
                                       } 
                                       
                                       group
                                       
                                     })







table(Plot.DF$Starting.PSI.Group, useNA = "always")
Starting.PSI.Means <- sapply(1:max(Plot.DF$Starting.PSI.Group),
                             function(x){
                               mean(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group == x)])
                             })

Plot.DF$Starting.PSI.Group <- factor(Plot.DF$Starting.PSI.Group, levels = 1:length(Starting.PSI.Means))

mean(abs(Plot.DF$Delta.PSI)[which(Plot.DF$Starting.PSI.Group==15)]<0.1)*100

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
  scale_fill_manual(values = map2color(x = Starting.PSI.Means[1:25], pal = My.Colour.Palette, limits = c(0,1))) +
  scale_colour_manual(values = rep("black", length(Starting.PSI.Means))) +
  theme_bw() +
  ylim(c(0,1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.3,
        legend.position = "none") +
  xlab("starting PSI") +
  ylab("final PSI") +
  ggtitle("DPSI distributions at different starting PSIs")+
  scale_x_discrete(labels=c(as.character(round(Starting.PSI.Means, 4)*100))) +
  annotate(geom = "text",
           x = max(as.numeric(Plot.DF$Starting.PSI.Group)),
           y = 0,
           label = "n = 12-110")

ggsave(filename = "201_Violin_Plots_GTEx_25_Groups_Exon1.pdf", plot = Violin.Plots, width = 12, height = 12, useDingbats = F)
































































































library(ggplot2)
Tissues <- c("Adipose Tissue",
             "Adrenal Gland",
             "Bladder",
             "Blood",
             "Blood Vessel",
             "Bone Marrow",
             "Brain",
             "Breast",
             "Cervix Uteri",
             "Colon",
             "Esophagus",
             "Fallopian Tube",
             "Heart",
             "Kidney",
             "Liver",
             "Lung",
             "Muscle",
             "Nerve",
             "Ovary",
             "Pancreas",
             "Pituitary",
             "Prostate",
             "Salivary Gland",
             "Skin",
             "Small Intestine",
             "Spleen",
             "Stomach",
             "Testis",
             "Thyroid",
             "Uterus",
             "Vagina")


All.Starting.PSIs <- c()
All.Mutation.Effects <- c()


for (each.tissue in Tissues) {
  print(each.tissue)
  
  tissue.no.spaces <- gsub(pattern = " ",
                           replacement = "",
                           x = each.tissue,
                           fixed = TRUE)
  
  file.location <- paste("./157_Exon2_Effects_vs_Starting_PSI/",
                         tissue.no.spaces,
                         ".RData",
                         sep = "")
  
  load(file.location)
  
  if (length(Starting.PSIs)==0){
    next
  }
  
  All.Starting.PSIs <- c(All.Starting.PSIs, Starting.PSIs)
  All.Mutation.Effects <- c(All.Mutation.Effects, Mutation.Effects)
  
}

All.Mutation.Effects <- All.Mutation.Effects*2

# negative mutations
Plot.DF <- data.frame(Starting.PSI = All.Starting.PSIs,
                      Delta.PSI = All.Mutation.Effects,
                      Final.PSI = All.Starting.PSIs + All.Mutation.Effects)
Plot.DF <- Plot.DF[-which(Plot.DF$Starting.PSI < 0),]
Plot.DF <- Plot.DF[-which(Plot.DF$Starting.PSI > 1),]

# Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI < 0)] <- 0
# Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI > 1)] <- 1




Plot.DF$Starting.PSI.Group <- sapply(X = Plot.DF$Starting.PSI,
                                     FUN = function(x){
                                       group <- NA
                                       
                                       if (x <= 0.001594883  ) {
                                         group <- 1
                                       } else if (x <= 0.004624439) {
                                         group <- 2
                                       } else if (x <= 0.010841853 ) {
                                         group <- 3
                                       } else if (x <= 0.023742260 ) {
                                         group <- 4
                                       } else if (x <= 0.049948917) {
                                         group <- 5
                                       } else if (x <= 0.112350678) {
                                         group <- 6
                                       } else if (x <= 0.295190228) {
                                         group <- 7
                                       } else if (x <= 0.657180986) {
                                         group <- 8
                                       } else if (x <= 0.867587457) {
                                         group <- 9
                                       } else if (x <= 0.940297361) {
                                         group <- 10
                                       } else if (x <= 0.969923057) {
                                         group <- 11
                                       } else if (x <= 0.985507246) {
                                         group <- 12
                                       } else if (x <= 0.993644307) {
                                         group <- 13
                                       } else if (x <= 0.997765643) {
                                         group <- 14
                                       } else if (x <= 1) {
                                         group <- 15
                                       } 
                                       
                                       group
                                       
                                     })







table(Plot.DF$Starting.PSI.Group, useNA = "always")
Starting.PSI.Means <- sapply(1:max(Plot.DF$Starting.PSI.Group),
                             function(x){
                               mean(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group == x)])
                             })

Plot.DF$Starting.PSI.Group <- factor(Plot.DF$Starting.PSI.Group, levels = 1:length(Starting.PSI.Means))

mean(abs(Plot.DF$Delta.PSI)[which(Plot.DF$Starting.PSI.Group==15)]<0.1)*100


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
  scale_fill_manual(values = map2color(x = Starting.PSI.Means[1:25], pal = My.Colour.Palette, limits = c(0,1))) +
  scale_colour_manual(values = rep("black", length(Starting.PSI.Means))) +
  theme_bw() +
  ylim(c(0,1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.3,
        legend.position = "none") +
  xlab("starting PSI") +
  ylab("final PSI") +
  ggtitle("DPSI distributions at different starting PSIs")+
  scale_x_discrete(labels=c(as.character(round(Starting.PSI.Means, 4)*100))) +
  annotate(geom = "text",
           x = max(as.numeric(Plot.DF$Starting.PSI.Group)),
           y = 0,
           label = "n = 3-107")

ggsave(filename = "201_Violin_Plots_GTEx_25_Groups_Exon2.pdf", plot = Violin.Plots, width = 12, height = 12, useDingbats = F)














