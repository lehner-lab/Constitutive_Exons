


library(ggplot2)
Kosuri <- read.table(file = "../Analyses/Data/snv_data_clean.txt", header = T)

Kosuri$mean_dpsi <- apply(X = Kosuri[,c("v1_dpsi", "v2_dpsi")],
                          MARGIN = 1,
                          FUN = mean, na.rm = T)

Kosuri$mean_nat_index <- apply(X = Kosuri[,c("nat_v1_index","nat_v2_index")],
                               MARGIN = 1,
                               FUN = mean, na.rm = T)


Kosuri$mutation.type <- apply(X = data.frame(Kosuri$strand,
                                             Kosuri$intron1_len,
                                             Kosuri$exon_len,
                                             Kosuri$intron2_len,
                                             Kosuri$rel_position),
                              MARGIN = 1,
                              FUN = function(x){
                                
                                strand <- as.character(x[1])
                                intron1.length <- as.numeric(x[2])
                                exon.length <- as.numeric(x[3])
                                intron2.length <- as.numeric(x[4])
                                relative.position <- as.numeric(x[5])
                                
                                if (strand == "-") {
                                  intron2.length <- as.numeric(x[2])
                                  intron1.length <- as.numeric(x[4])
                                }
                                
                                mutation.class <- "other"
                                
                                if (! is.na(relative.position)) {
                                  if (relative.position %in% c(intron1.length, (intron1.length-1))){
                                    mutation.class <- "acceptor.site"
                                  }
                                  
                                  if (relative.position %in% c((intron1.length+exon.length+1),(intron1.length+exon.length+2))) {
                                    mutation.class <- "donor.site"
                                  }
                                  
                                  if (relative.position == (intron1.length+1)){
                                    mutation.class <- "exon.1"
                                  }
                                  
                                  if (relative.position == (intron1.length+exon.length)){
                                    mutation.class <- "exon.2"
                                  }
                                }
                                
                                mutation.class
                                
                              })

table(Kosuri$mutation.type)

Kosuri.Subtable <- Kosuri[,c("id",
                             "sub_id",
                             "mean_dpsi",
                             "mean_nat_index",
                             "mutation.type")]












SpliceSite.Mutations <- Kosuri.Subtable
SpliceSite.Mutations <- SpliceSite.Mutations[which(SpliceSite.Mutations$mutation.type == "exon.1"),]



SNV.Dataset <- data.frame(Delta.PSI = SpliceSite.Mutations$mean_dpsi,
                          Starting.PSI = SpliceSite.Mutations$mean_nat_index,
                          ID = SpliceSite.Mutations$id)








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














# same bins as intron dataset
Plot.DF$Starting.PSI.Group <- sapply(X = Plot.DF$Starting.PSI,
                                     FUN = function(x){
                                       group <- NA
                                       
                                       if (x <= 0.8926642) {
                                         group <- 1
                                       } else if (x <= 0.9440571) {
                                         group <- 2
                                       } else if (x <= 0.9664795) {
                                         group <- 3
                                       } else if (x <= 0.9835943) {
                                         group <- 4
                                       } else if (x <= 1) {
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
           label = "n = 24-51")

ggsave(filename = "198_Violin_Plots_SNV_Exon1.pdf", plot = Violin.Plots, width = 8, height = 8, useDingbats = F)
































































































SpliceSite.Mutations <- Kosuri.Subtable
SpliceSite.Mutations <- SpliceSite.Mutations[which(SpliceSite.Mutations$mutation.type == "exon.2"),]



SNV.Dataset <- data.frame(Delta.PSI = SpliceSite.Mutations$mean_dpsi,
                          Starting.PSI = SpliceSite.Mutations$mean_nat_index,
                          ID = SpliceSite.Mutations$id)








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














# same bins as intron dataset
Plot.DF$Starting.PSI.Group <- sapply(X = Plot.DF$Starting.PSI,
                                     FUN = function(x){
                                       group <- NA
                                       
                                       if (x <= 0.8926642) {
                                         group <- 1
                                       } else if (x <= 0.9440571) {
                                         group <- 2
                                       } else if (x <= 0.9664795) {
                                         group <- 3
                                       } else if (x <= 0.9835943) {
                                         group <- 4
                                       } else if (x <= 1) {
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
           label = "n = 27-50")

ggsave(filename = "198_Violin_Plots_SNV_Exon2.pdf", plot = Violin.Plots, width = 8, height = 8, useDingbats = F)

