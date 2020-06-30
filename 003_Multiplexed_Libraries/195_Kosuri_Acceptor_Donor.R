
# SRE Dataset
Kosuri <- read.table(file = "../Analyses/Data/sre_data_clean.txt", header = T)

Kosuri <- Kosuri[which(Kosuri$category %in% c("rnd_intron_1nt",
                                              "rnd_exon_1nt")),]



Kosuri$mutation.type <- apply(X = data.frame(Kosuri$strand,
                                             Kosuri$intron1_len,
                                             Kosuri$exon_len,
                                             Kosuri$intron2_len,
                                             as.character(Kosuri$mutations)),
                              MARGIN = 1,
                              FUN = function(x){
                                strand <- as.numeric(x[1])
                                intron1.length <- as.numeric(x[2])
                                exon.length <- as.numeric(x[3])
                                intron2.length <- as.numeric(x[4])
                                mutations <- as.character(x[5])
                                mutation.position <- as.numeric(strsplit(strsplit(mutations, ",")[[1]][2], ")")[[1]][1])
                                
                                if (strand == -1) {
                                  intron2.length <- as.numeric(x[2])
                                  intron1.length <- as.numeric(x[4])
                                }
                                
                                mutation.type = "other"
                                
                                if (mutation.position %in% c(intron1.length,(intron1.length-1))){
                                  mutation.type <- "acceptor"
                                }
                                
                                if (mutation.position %in% c((intron1.length+exon.length+1),(intron1.length+exon.length+2))){
                                  mutation.type <- "donor"
                                }
                                
                                # if (mutation.position == (intron1.length+1)) {
                                #   mutation.type <- "exon.1"
                                # }
                                # 
                                # if (mutation.position == (intron1.length+exon.length)){
                                #   mutation.type <- "exon.2"
                                # }
                                
                                mutation.type
                                
                                
                              })


table(Kosuri$mutation.type)

hist(Kosuri$dpsi_dhfr[Kosuri$mutation.type == "exon.2"])

Kosuri.So.Far <- Kosuri







Kosuri <- read.table(file = "../Analyses/Data/sre_data_clean.txt", header = T)

Kosuri <- Kosuri[which(Kosuri$category %in% c("rnd_intron_2nt",
                                              "rnd_exon_2nt")),]

# 
# as.numeric(strsplit(strsplit("[((9,10),C),((38,39),T)]", ",")[[1]][2], ")")[[1]][1])
# as.numeric(strsplit(strsplit("[((9,10),C),((38,39),T)]", ",")[[1]][5], ")")[[1]][1])

Kosuri$mutation.type <- apply(X = data.frame(Kosuri$strand,
                                             Kosuri$intron1_len,
                                             Kosuri$exon_len,
                                             Kosuri$intron2_len,
                                             as.character(Kosuri$mutations)),
                              MARGIN = 1,
                              FUN = function(x){
                                strand <- as.numeric(x[1])
                                intron1.length <- as.numeric(x[2])
                                exon.length <- as.numeric(x[3])
                                intron2.length <- as.numeric(x[4])
                                mutations <- as.character(x[5])
                                mutation.position.1 <- as.numeric(strsplit(strsplit(mutations, ",")[[1]][2], ")")[[1]][1])
                                mutation.position.2 <- as.numeric(strsplit(strsplit(mutations, ",")[[1]][5], ")")[[1]][1])
                                mutation.position <- c(mutation.position.1, mutation.position.2)
                                if (strand == -1) {
                                  intron2.length <- as.numeric(x[2])
                                  intron1.length <- as.numeric(x[4])
                                }
                                
                                mutation.type = "other"
                                
                                if (any(mutation.position %in% c(intron1.length,(intron1.length-1)))){
                                  mutation.type <- "acceptor"
                                }
                                
                                if (any(mutation.position %in% c((intron1.length+exon.length+1),(intron1.length+exon.length+2)))){
                                  mutation.type <- "donor"
                                }
                                
                                # if (any(mutation.position == (intron1.length+1))) {
                                #   mutation.type <- "exon.1"
                                # }
                                # 
                                # if (any(mutation.position == (intron1.length+exon.length))){
                                #   mutation.type <- "exon.2"
                                # }
                                
                                mutation.type
                                
                                
                              })


Kosuri.So.Far <- rbind(Kosuri.So.Far,
                       Kosuri)






Kosuri <- read.table(file = "../Analyses/Data/sre_data_clean.txt", header = T)

Kosuri <- Kosuri[which(Kosuri$category %in% c("rnd_intron_3nt",
                                              "rnd_exon_3nt")),]

# 
# as.numeric(strsplit(strsplit("[((138,139),C),((140,141),G),((146,147),A)]", ",")[[1]][2], ")")[[1]][1])
# as.numeric(strsplit(strsplit("[((138,139),C),((140,141),G),((146,147),A)]", ",")[[1]][5], ")")[[1]][1])
# as.numeric(strsplit(strsplit("[((138,139),C),((140,141),G),((146,147),A)]", ",")[[1]][8], ")")[[1]][1])

Kosuri$mutation.type <- apply(X = data.frame(Kosuri$strand,
                                             Kosuri$intron1_len,
                                             Kosuri$exon_len,
                                             Kosuri$intron2_len,
                                             as.character(Kosuri$mutations)),
                              MARGIN = 1,
                              FUN = function(x){
                                strand <- as.numeric(x[1])
                                intron1.length <- as.numeric(x[2])
                                exon.length <- as.numeric(x[3])
                                intron2.length <- as.numeric(x[4])
                                mutations <- as.character(x[5])
                                mutation.position.1 <- as.numeric(strsplit(strsplit(mutations, ",")[[1]][2], ")")[[1]][1])
                                mutation.position.2 <- as.numeric(strsplit(strsplit(mutations, ",")[[1]][5], ")")[[1]][1])
                                mutation.position.3 <- as.numeric(strsplit(strsplit(mutations, ",")[[1]][8], ")")[[1]][1])
                                mutation.position <- c(mutation.position.1, mutation.position.2, mutation.position.3)
                                if (strand == -1) {
                                  intron2.length <- as.numeric(x[2])
                                  intron1.length <- as.numeric(x[4])
                                }
                                
                                mutation.type = "other"
                                
                                if (any(mutation.position %in% c(intron1.length,(intron1.length-1)))){
                                  mutation.type <- "acceptor"
                                }
                                
                                if (any(mutation.position %in% c((intron1.length+exon.length+1),(intron1.length+exon.length+2)))){
                                  mutation.type <- "donor"
                                }
                                
                                # if (any(mutation.position == (intron1.length+1))) {
                                #   mutation.type <- "exon.1"
                                # }
                                # 
                                # if (any(mutation.position == (intron1.length+exon.length))){
                                #   mutation.type <- "exon.2"
                                # }
                                
                                mutation.type
                                
                                
                              })


Kosuri.So.Far <- rbind(Kosuri.So.Far,
                       Kosuri)


























Kosuri <- read.table(file = "../Analyses/Data/sre_data_clean.txt", header = T)

Kosuri <- Kosuri[which(Kosuri$category %in% c("rnd_intron_5nt",
                                              "rnd_exon_5nt")),]

# 
# as.numeric(strsplit(strsplit("[((1,2),T),((14,15),G),((28,29),G),((37,38),T),((38,39),G)]", ",")[[1]][2], ")")[[1]][1])
# as.numeric(strsplit(strsplit("[((1,2),T),((14,15),G),((28,29),G),((37,38),T),((38,39),G)]", ",")[[1]][5], ")")[[1]][1])
# as.numeric(strsplit(strsplit("[((1,2),T),((14,15),G),((28,29),G),((37,38),T),((38,39),G)]", ",")[[1]][8], ")")[[1]][1])
# as.numeric(strsplit(strsplit("[((1,2),T),((14,15),G),((28,29),G),((37,38),T),((38,39),G)]", ",")[[1]][11], ")")[[1]][1])
# as.numeric(strsplit(strsplit("[((1,2),T),((14,15),G),((28,29),G),((37,38),T),((38,39),G)]", ",")[[1]][14], ")")[[1]][1])

Kosuri$mutation.type <- apply(X = data.frame(Kosuri$strand,
                                             Kosuri$intron1_len,
                                             Kosuri$exon_len,
                                             Kosuri$intron2_len,
                                             as.character(Kosuri$mutations)),
                              MARGIN = 1,
                              FUN = function(x){
                                strand <- as.numeric(x[1])
                                intron1.length <- as.numeric(x[2])
                                exon.length <- as.numeric(x[3])
                                intron2.length <- as.numeric(x[4])
                                mutations <- as.character(x[5])
                                mutation.position.1 <- as.numeric(strsplit(strsplit(mutations, ",")[[1]][2], ")")[[1]][1])
                                mutation.position.2 <- as.numeric(strsplit(strsplit(mutations, ",")[[1]][5], ")")[[1]][1])
                                mutation.position.3 <- as.numeric(strsplit(strsplit(mutations, ",")[[1]][8], ")")[[1]][1])
                                mutation.position.4 <- as.numeric(strsplit(strsplit(mutations, ",")[[1]][11], ")")[[1]][1])
                                mutation.position.5 <- as.numeric(strsplit(strsplit(mutations, ",")[[1]][14], ")")[[1]][1])
                                mutation.position <- c(mutation.position.1,
                                                       mutation.position.2,
                                                       mutation.position.3,
                                                       mutation.position.4,
                                                       mutation.position.5)
                                if (strand == -1) {
                                  intron2.length <- as.numeric(x[2])
                                  intron1.length <- as.numeric(x[4])
                                }
                                
                                mutation.type = "other"
                                
                                if (any(mutation.position %in% c(intron1.length,(intron1.length-1)))){
                                  mutation.type <- "acceptor"
                                }
                                
                                if (any(mutation.position %in% c((intron1.length+exon.length+1),(intron1.length+exon.length+2)))){
                                  mutation.type <- "donor"
                                }
                                
                                # if (any(mutation.position == (intron1.length+1))) {
                                #   mutation.type <- "exon.1"
                                # }
                                # 
                                # if (any(mutation.position == (intron1.length+exon.length))){
                                #   mutation.type <- "exon.2"
                                # }
                                
                                mutation.type
                                
                                
                              })


Kosuri.So.Far <- rbind(Kosuri.So.Far,
                       Kosuri)
































Kosuri.Subtable <- Kosuri.So.Far[,c("id",
                                    "sub_id",
                                    "mutation.type",
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


SpliceSite.Class <- c("acceptor",
                      "donor")

# subset the data so it only includes splice site mutations that decrease splicing (SMN1 dataset)
SpliceSite.Mutations <- Kosuri.Subtable
SpliceSite.Mutations <- SpliceSite.Mutations[which(SpliceSite.Mutations$mutation.type %in% SpliceSite.Class),]






SMN1.Dataset <- data.frame(Delta.PSI = SpliceSite.Mutations$dpsi_smn1,
                           Starting.PSI = SpliceSite.Mutations$nat_index_smn1,
                           ID = SpliceSite.Mutations$id)



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

















# same bins as intron dataset
Plot.DF$Starting.PSI.Group <- sapply(X = Plot.DF$Starting.PSI,
                                     FUN = function(x){
                                       group <- NA
                                       
                                       if (x <= 0.2363344) {
                                         group <- 1
                                       } else if (x <= 0.8242058) {
                                         group <- 2
                                       } else if (x <= 0.9062518) {
                                         group <- 3
                                       } else if (x <= 0.9467333) {
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
           label = "n = 9-16")

ggsave(filename = "195_Violin_Plots_SMN1.pdf", plot = Violin.Plots, width = 8, height = 8, useDingbats = F)























































































DHFR.Dataset <- data.frame(Delta.PSI = SpliceSite.Mutations$dpsi_dhfr,
                           Starting.PSI = SpliceSite.Mutations$nat_index_dhfr,
                           ID = SpliceSite.Mutations$id)


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













# same bins as intron dataset
Plot.DF$Starting.PSI.Group <- sapply(X = Plot.DF$Starting.PSI,
                                     FUN = function(x){
                                       group <- NA
                                       
                                       if (x <= 0.04178355 ) {
                                         group <- 1
                                       } else if (x <= 0.80554968 ) {
                                         group <- 2
                                       } else if (x <= 0.92305210) {
                                         group <- 3
                                       } else if (x <= 0.96636764) {
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
           label = "n = 10-13")

ggsave(filename = "195_Violin_Plots_DHFR.pdf", plot = Violin.Plots, width = 8, height = 8, useDingbats = F)



















































































































































































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








Kosuri.Subtable <- Kosuri.Subtable[complete.cases(Kosuri.Subtable),]


Kosuri.Subtable$id <- sapply(X = as.character(Kosuri.Subtable$id),
                             FUN = function(x){
                               strsplit(x, "_")[[1]][1]
                             })

Number.Of.Variants.Per.Exon <- table(Kosuri.Subtable$id) - 1


IDs.To.Keep <- names(Number.Of.Variants.Per.Exon)[which(Number.Of.Variants.Per.Exon > 0)]
Kosuri.Subtable <- Kosuri.Subtable[which(Kosuri.Subtable$id %in% IDs.To.Keep),]


SpliceSite.Mutations <- Kosuri.Subtable
SpliceSite.Mutations <- SpliceSite.Mutations[which(SpliceSite.Mutations$mutation.type == "acceptor.site" | SpliceSite.Mutations$mutation.type == "donor.site"),]



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
                                       
                                       if (x <= 0.8921605  ) {
                                         group <- 1
                                       } else if (x <= 0.9440571  ) {
                                         group <- 2
                                       } else if (x <= 0.9663807) {
                                         group <- 3
                                       } else if (x <= 0.9834099) {
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
           label = "n = 35-59")

ggsave(filename = "195_Violin_Plots_SNV.pdf", plot = Violin.Plots, width = 8, height = 8, useDingbats = F)





















