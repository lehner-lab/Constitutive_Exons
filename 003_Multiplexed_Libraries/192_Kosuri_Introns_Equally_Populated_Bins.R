
source("185_Functions.R")
library(ggplot2)

# SRE Dataset
Kosuri <- read.table(file = "../Analyses/Data/sre_data_clean.txt", header = T)

table(Kosuri$category)

Kosuri.Subtable <- Kosuri[,c("id",
                             "sub_id",
                             "category",
                             "dpsi_smn1",
                             "dpsi_dhfr",
                             "nat_index_smn1",
                             "nat_index_dhfr")]

Kosuri.Subtable <- Kosuri.Subtable[complete.cases(Kosuri.Subtable),]

Kosuri.Subtable$id <- sapply(X = as.character(Kosuri.Subtable$id),
                             FUN = function(x){
                               strsplit(x, "_")[[1]][1]
                             })

Number.Of.Variants.Per.Exon <- table(Kosuri.Subtable$id) - 1

IDs.To.Keep <- names(Number.Of.Variants.Per.Exon)[which(Number.Of.Variants.Per.Exon > 0)]
Kosuri.Subtable <- Kosuri.Subtable[which(Kosuri.Subtable$id %in% IDs.To.Keep),]







# define exonic class mutations
Exonic.Class <- c("clst_Vlkr07_AICS",
                  "clst_Vlkr07_DICS",
                  "rmv_Vlkr07_AICS",
                  "rmv_Vlkr07_DICS",
                  "rnd_intron_1nt",
                  "rnd_intron_2nt",
                  "rnd_intron_3nt",
                  "rnd_intron_5nt",
                  "aggr_intron",
                  "p_aggr_intr")


# subset the data so it only includes exonic mutations that decrease splicing (SMN1 dataset)
Exonic.Mutations <- Kosuri.Subtable
Exonic.Mutations <- Exonic.Mutations[which(Exonic.Mutations$category %in% Exonic.Class),]





















DHFR.Dataset <- data.frame(Delta.PSI = Exonic.Mutations$dpsi_dhfr,
                           Starting.PSI = Exonic.Mutations$nat_index_dhfr,
                           ID = Exonic.Mutations$id)


List.Mutation.Effects <- vector(mode = "list",
                                length = length(unique(DHFR.Dataset$ID)))
names(List.Mutation.Effects) <- as.character(unique(DHFR.Dataset$ID))

for (i in 1:nrow(DHFR.Dataset)){
  this.id <- as.character(DHFR.Dataset$ID)[i]
  List.Mutation.Effects[[this.id]] <- c(List.Mutation.Effects[[this.id]],
                                        DHFR.Dataset$Delta.PSI[i])
}

Median.Effects <- sapply(List.Mutation.Effects, median)

Starting.PSIs <- vector(mode = "numeric",
                        length = length(Median.Effects))
names(Starting.PSIs) <- as.character(names(Median.Effects))

for (each.name in as.character(names(Starting.PSIs))){
  idx <- which(DHFR.Dataset$ID == each.name)[1]
  Starting.PSIs[each.name] <- DHFR.Dataset$Starting.PSI[idx]
}


Plot.DF <- data.frame(Starting.PSI = Starting.PSIs,
                      Delta.PSI = Median.Effects,
                      Final.PSI = Starting.PSIs + Median.Effects)















#Plot.DF$Starting.PSI.Group <- as.numeric(cut(x = Plot.DF$Starting.PSI, breaks = seq(0,1,0.2), include.lowest = T))


Plot.DF$Starting.PSI.Group <- EquallyPopulatedBins(vector = Plot.DF$Starting.PSI, k = 5)
table(Plot.DF$Starting.PSI.Group, useNA = "always")
Starting.PSI.Means <- sapply(1:max(Plot.DF$Starting.PSI.Group),
                             function(x){
                               mean(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group == x)])
                             })

Plot.DF$Starting.PSI.Group <- factor(Plot.DF$Starting.PSI.Group, levels = 1:length(Starting.PSI.Means))
sapply(1:5, function(x){ max(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group==x)])})



mean(abs(Plot.DF$Delta.PSI[which(Plot.DF$Starting.PSI.Group == 5)]) > 0.1)*100
median(abs(Plot.DF$Delta.PSI[which(Plot.DF$Starting.PSI.Group == 5)]))*100

mean(abs(Plot.DF$Delta.PSI[which(Plot.DF$Starting.PSI.Group == 2)]) > 0.1)*100
median(abs(Plot.DF$Delta.PSI[which(Plot.DF$Starting.PSI.Group == 2)]))*100



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
  scale_fill_manual(values = map2color(x = Starting.PSI.Means[1:5], pal = My.Colour.Palette, limits = c(0,1))) +
  scale_colour_manual(values = rep("black", length(Starting.PSI.Means))) +
  theme_bw() +
  ylim(c(0,1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.7,
        legend.position = "none") +
  xlab("starting PSI") +
  ylab("final PSI") +
  ggtitle("DPSI distributions at different starting PSIs")+
  scale_x_discrete(labels=c(as.character(round(Starting.PSI.Means, 4)*100))) +
  annotate(geom = "text",
           x = max(as.numeric(Plot.DF$Starting.PSI.Group)),
           y = 0,
           label = "n = 27")

ggsave(filename = "192_Violin_Plots_DHFR.pdf", plot = Violin.Plots, width = 8, height = 8, useDingbats = F)

























































































SMN1.Dataset <- data.frame(Delta.PSI = Exonic.Mutations$dpsi_smn1,
                           Starting.PSI = Exonic.Mutations$nat_index_smn1,
                           ID = Exonic.Mutations$id)


List.Mutation.Effects <- vector(mode = "list",
                                length = length(unique(SMN1.Dataset$ID)))
names(List.Mutation.Effects) <- as.character(unique(SMN1.Dataset$ID))

for (i in 1:nrow(SMN1.Dataset)){
  this.id <- as.character(SMN1.Dataset$ID)[i]
  List.Mutation.Effects[[this.id]] <- c(List.Mutation.Effects[[this.id]],
                                        SMN1.Dataset$Delta.PSI[i])
}

Median.Effects <- sapply(List.Mutation.Effects, median)

Starting.PSIs <- vector(mode = "numeric",
                        length = length(Median.Effects))
names(Starting.PSIs) <- as.character(names(Median.Effects))

for (each.name in as.character(names(Starting.PSIs))){
  idx <- which(SMN1.Dataset$ID == each.name)[1]
  Starting.PSIs[each.name] <- SMN1.Dataset$Starting.PSI[idx]
}


Plot.DF <- data.frame(Starting.PSI = Starting.PSIs,
                      Delta.PSI = Median.Effects,
                      Final.PSI = Starting.PSIs + Median.Effects)















#Plot.DF$Starting.PSI.Group <- as.numeric(cut(x = Plot.DF$Starting.PSI, breaks = seq(0,1,0.2), include.lowest = T))


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
  scale_fill_manual(values = map2color(x = Starting.PSI.Means[1:5], pal = My.Colour.Palette, limits = c(0,1))) +
  scale_colour_manual(values = rep("black", length(Starting.PSI.Means))) +
  theme_bw() +
  ylim(c(0,1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.7,
        legend.position = "none") +
  xlab("starting PSI") +
  ylab("final PSI") +
  ggtitle("DPSI distributions at different starting PSIs")+
  scale_x_discrete(labels=c(as.character(round(Starting.PSI.Means, 4)*100))) +
  annotate(geom = "text",
           x = max(as.numeric(Plot.DF$Starting.PSI.Group)),
           y = 0,
           label = "n = 27")

ggsave(filename = "192_Violin_Plots_SMN1.pdf", plot = Violin.Plots, width = 8, height = 8, useDingbats = F)


















































































































































































library(ggplot2)
Kosuri <- read.table(file = "../Analyses/Data/snv_data_clean.txt", header = T)

Kosuri$mean_dpsi <- apply(X = Kosuri[,c("v1_dpsi", "v2_dpsi")],
                          MARGIN = 1,
                          FUN = mean, na.rm = T)

Kosuri$mean_nat_index <- apply(X = Kosuri[,c("nat_v1_index","nat_v2_index")],
                               MARGIN = 1,
                               FUN = mean, na.rm = T)

Kosuri.Subtable <- Kosuri[,c("id",
                             "sub_id",
                             "mean_dpsi",
                             "mean_nat_index",
                             "label")]

Kosuri.Subtable <- Kosuri.Subtable[complete.cases(Kosuri.Subtable),]


Kosuri.Subtable$id <- sapply(X = as.character(Kosuri.Subtable$id),
                             FUN = function(x){
                               strsplit(x, "_")[[1]][1]
                             })

Number.Of.Variants.Per.Exon <- table(Kosuri.Subtable$id) - 1


IDs.To.Keep <- names(Number.Of.Variants.Per.Exon)[which(Number.Of.Variants.Per.Exon > 0)]
Kosuri.Subtable <- Kosuri.Subtable[which(Kosuri.Subtable$id %in% IDs.To.Keep),]


Exonic.Mutations <- Kosuri.Subtable
Exonic.Mutations <- Exonic.Mutations[which(Exonic.Mutations$label %in% c("downstr_intron", "upstr_intron")),]




SNV.Dataset <- data.frame(Delta.PSI = Exonic.Mutations$mean_dpsi,
                          Starting.PSI = Exonic.Mutations$mean_nat_index,
                          ID = Exonic.Mutations$id)



List.Mutation.Effects <- vector(mode = "list",
                                length = length(unique(SNV.Dataset$ID))) 
names(List.Mutation.Effects) <- as.character(unique(SNV.Dataset$ID))

for (i in 1:nrow(SNV.Dataset)){
  this.id <- as.character(SNV.Dataset$ID)[i]
  List.Mutation.Effects[[this.id]] <- c(List.Mutation.Effects[[this.id]],
                                        SNV.Dataset$Delta.PSI[i])
}

Median.Effects <- sapply(List.Mutation.Effects, median)

Starting.PSIs <- vector(mode = "numeric",
                        length = length(Median.Effects))
names(Starting.PSIs) <- as.character(names(Median.Effects))

for (each.name in as.character(names(Starting.PSIs))){
  idx <- which(SNV.Dataset$ID == each.name)[1]
  Starting.PSIs[each.name] <- SNV.Dataset$Starting.PSI[idx]
}


Plot.DF <- data.frame(Starting.PSI = Starting.PSIs,
                      Delta.PSI = Median.Effects,
                      Final.PSI = Starting.PSIs + Median.Effects)












Plot.DF$Starting.PSI.Group <- EquallyPopulatedBins(vector = Plot.DF$Starting.PSI, k = 5)
table(Plot.DF$Starting.PSI.Group, useNA = "always")
Starting.PSI.Means <- sapply(1:max(Plot.DF$Starting.PSI.Group),
                             function(x){
                               mean(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group == x)])
                             })

Plot.DF$Starting.PSI.Group <- factor(Plot.DF$Starting.PSI.Group, levels = 1:length(Starting.PSI.Means))
sapply(1:5, function(x){ max(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group==x)])})


mean(abs(Plot.DF$Delta.PSI[which(Plot.DF$Starting.PSI.Group == 5)]) > 0.1)*100
median(abs(Plot.DF$Delta.PSI[which(Plot.DF$Starting.PSI.Group == 5)]))*100

mean(abs(Plot.DF$Delta.PSI[which(Plot.DF$Starting.PSI.Group == 1)]) > 0.1)*100
median(abs(Plot.DF$Delta.PSI[which(Plot.DF$Starting.PSI.Group == 1)]))*100



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
  scale_fill_manual(values = map2color(x = Starting.PSI.Means[1:5], pal = My.Colour.Palette, limits = c(0,1))) +
  scale_colour_manual(values = rep("black", length(Starting.PSI.Means))) +
  theme_bw() +
  ylim(c(0,1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.7,
        legend.position = "none") +
  xlab("starting PSI") +
  ylab("final PSI") +
  ggtitle("DPSI distributions at different starting PSIs")+
  scale_x_discrete(labels=c(as.character(round(Starting.PSI.Means, 4)*100))) +
  annotate(geom = "text",
           x = max(as.numeric(Plot.DF$Starting.PSI.Group)),
           y = 0,
           label = "n = 415-416")

ggsave(filename = "192_Violin_Plots_SNV.pdf", plot = Violin.Plots, width = 8, height = 8, useDingbats = F)












