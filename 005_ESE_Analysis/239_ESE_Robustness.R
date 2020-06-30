
library(data.table)
load("239bis_Hexamer_Table.RData")

# reverse complement function
rev.comp<-function(x,rev=TRUE)
{
  x<-toupper(x)
  y<-rep("N",nchar(x))
  xx<-unlist(strsplit(x,NULL))
  for (bbb in 1:nchar(x))
  {
    if(xx[bbb]=="A") y[bbb]<-"T"    
    if(xx[bbb]=="C") y[bbb]<-"G"    
    if(xx[bbb]=="G") y[bbb]<-"C"    
    if(xx[bbb]=="T") y[bbb]<-"A"
  }
  if(rev==FALSE) 
  {
    for(ccc in (1:nchar(x)))
    {
      if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
    }
  }
  if(rev==T)
  {
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x)))
    {
      zz[ccc]<-y[nchar(x)+1-ccc]
      if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
    }
  }
  return(yy)  
}


# load table
Bed.Table <- read.table(file = "238_GTEx_Exon_Regions_With_Sequences.bed",
                        header = F,
                        sep = "\t",
                        col.names = c("chromosome",
                                      "start",
                                      "end",
                                      "id",
                                      "strand",
                                      "gene",
                                      "length",
                                      "sequence"),
                        stringsAsFactors = F)


# reverse complement any sequence in the reverse strand
Bed.Table$sequence <- apply(X = Bed.Table[,c("sequence", "strand")],
                            MARGIN = 1,
                            FUN = function(x){
                              mySequence <- x[1]
                              myStrand <- x[2]
                              
                              if (myStrand == "-") {
                                mySequence <- rev.comp(mySequence)
                              }
                              
                              mySequence
                            })





# function to count number of overlapping matches
find_overlaps <- function(p,s) {
  gg <- gregexpr(paste0("(?=",p,")"),s,perl=TRUE)[[1]]
  if (length(gg)==1 && gg==-1) 0 else length(gg)
}






Bed.Table$ese.count <- sapply(X = Bed.Table$sequence,
                              FUN = function(mySequence){
                                sum(sapply(X = ESE$Hexamer,
                                           FUN = function(ese){
                                             find_overlaps(p = ese,
                                                           s = mySequence)
                                           }))
                              })



Bed.Table$ese.matches <- sapply(X = Bed.Table$sequence,
                                FUN = function(mySequence){
                                  matches <- sapply(X = ESE$Hexamer,
                                                    FUN = function(ese){
                                                      find_overlaps(p = ese,
                                                                    s = mySequence)
                                                    })
                                  matches <- names(matches[which(matches>0)])
                                  paste(matches,
                                        sep = ",",
                                        collapse = ",")
                                })




# find what % of hexamers one SNP away from each ESE are still ESEs
PercentESEs1HammingDistanceAway <- function(this.hexamer, hexamer.network = Hexamer.Network, ese = ESE){
  hexamers.1.hamming.distance.away <- hexamer.network[Main == this.hexamer, Neighbour]
  number.ese.1.hamming.distance.away <- length(which(hexamers.1.hamming.distance.away %in% ese$Hexamer))
  number.ese.1.hamming.distance.away/length(hexamers.1.hamming.distance.away)
}



Bed.Table$neighbourhood.eses <- sapply(X = Bed.Table$ese.matches,
                                       FUN = function(x){
                                         
                                         
                                         if (x != "") {
                                           ese.here <- strsplit(x = x, split = ",")[[1]]
                                           results <- sapply(X = ese.here,
                                                            FUN = function(y){
                                                              PercentESEs1HammingDistanceAway(this.hexamer = y)
                                                            })
                                         } else {
                                           results <- NA
                                         }
                                         
                                         results
                                       })








































































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



First.Iteration <- TRUE

for (each.tissue in Tissues) {
  
  print(each.tissue)
  
  tissue.no.spaces <- gsub(pattern = " ",
                           replacement = "",
                           x = each.tissue,
                           fixed = TRUE)
  
  file.location <- paste("182_GTEx_PSI_Distributions_Plus_IDs/ALL_PSIandID_",
                         tissue.no.spaces,
                         ".RData",
                         sep = "")
  
  load(file.location)
  
  if (length(Mean.PSI.Values)==0){
    next
  }
  
  
  if (First.Iteration) {
    
    Constitutive.Alternative.Exon.Matrix <- matrix(data = 0,
                                                   nrow = length(Splicing.Events),
                                                   ncol = 3,
                                                   dimnames = list(Splicing.Events,
                                                                   c("Number.Of.Tissues",
                                                                     "Times.Under.60",
                                                                     "Times.Over.90"))
    )
    
    First.Iteration <- FALSE
  }
  
  
  if (any(is.nan(Mean.PSI.Values))) {
    idx <- which(is.nan(Mean.PSI.Values))
    Mean.PSI.Values <- Mean.PSI.Values[-idx]
    Splicing.Events <- Splicing.Events[-idx]
  } 
  
  Constitutive.Alternative.Exon.Matrix[Splicing.Events,"Number.Of.Tissues"] <- Constitutive.Alternative.Exon.Matrix[Splicing.Events,"Number.Of.Tissues"] + 1
  Constitutive.Alternative.Exon.Matrix[Splicing.Events[which(Mean.PSI.Values<0.6)],"Times.Under.60"] <- Constitutive.Alternative.Exon.Matrix[Splicing.Events[which(Mean.PSI.Values<0.6)],"Times.Under.60"] + 1
  Constitutive.Alternative.Exon.Matrix[Splicing.Events[which(Mean.PSI.Values>0.9)],"Times.Over.90"] <- Constitutive.Alternative.Exon.Matrix[Splicing.Events[which(Mean.PSI.Values>0.9)],"Times.Over.90"] + 1
  
}


Constitutive.Alternative.Exon.DF <- as.data.frame(Constitutive.Alternative.Exon.Matrix)
Constitutive.Alternative.Exon.DF$Exon.Class <- apply(X = Constitutive.Alternative.Exon.Matrix,
                                                     MARGIN = 1,
                                                     FUN = function(x){
                                                       result <- "other"
                                                       if (x[2] > 0) {
                                                         result <- "alternative"
                                                       } else if (x[3] == x[1]) {
                                                         result <- "constitutive"
                                                       }
                                                       result
                                                     })








Bed.Table$Exon.Class <- sapply(X = Bed.Table$id,
                               FUN = function(x){
                                 Constitutive.Alternative.Exon.DF[x,"Exon.Class"]
                               })




save(Bed.Table, file = "239_ESE_Robustness_Bed_Table.RData")
load("239_ESE_Robustness_Bed_Table.RData")




Bed.Table.Alternative <- Bed.Table[Bed.Table$Exon.Class == "alternative",]
Bed.Table.Constitutive <- Bed.Table[Bed.Table$Exon.Class == "constitutive",]



Percent.Neighbourhood.ESEs.Alternative <- unname(unlist(Bed.Table.Alternative$neighbourhood.eses))
Percent.Neighbourhood.ESEs.Constitutive <- unname(unlist(Bed.Table.Constitutive$neighbourhood.eses))


mean(Percent.Neighbourhood.ESEs.Alternative, na.rm = T)
mean(Percent.Neighbourhood.ESEs.Constitutive, na.rm = T)

plot(density(Percent.Neighbourhood.ESEs.Alternative, na.rm = T, bw = 0.05))
plot(density(Percent.Neighbourhood.ESEs.Constitutive, na.rm = T, bw = 0.05))






















# calculate robustness of each ESE hexamer
Robustness.Each.ESE.Hexamer <- sapply(X = ESE$Hexamer,
       FUN = function(x){
         PercentESEs1HammingDistanceAway(this.hexamer = x, hexamer.network = Hexamer.Network, ese = ESE)
       })
table(Robustness.Each.ESE.Hexamer)

Unique.Percentages <- rev(sort(unique(Robustness.Each.ESE.Hexamer)))

Hexamers.Sorted.Into.Robustness.Groups <- vector(mode = "list", length = 18)
names(Hexamers.Sorted.Into.Robustness.Groups) <- 1:18



for (i in 1:length(Robustness.Each.ESE.Hexamer)) {
  Robustness.Group <- which(Unique.Percentages == Robustness.Each.ESE.Hexamer[i])
  Hexamers.Sorted.Into.Robustness.Groups[[Robustness.Group]] <- c(Hexamers.Sorted.Into.Robustness.Groups[[Robustness.Group]],
                                                                  names(Robustness.Each.ESE.Hexamer[i]))

}







Bed.Table.Alternative$ese.matches[1]



EseInSequence <- function(ese, ese.matches) {
  ese.matches <- strsplit(x = ese.matches,
                          split = ",")[[1]]
  ese %in% ese.matches
}




mean(sapply(X = Bed.Table.Constitutive$ese.matches,
            FUN = function(y){
              EseInSequence(ese = "GACGTC", ese.matches = y)
            }))


mean(sapply(X = Bed.Table.Alternative$ese.matches,
            FUN = function(y){
              EseInSequence(ese = "GACGTC", ese.matches = y)
            }))






EseGroupInSequence <- function(ese.group, ese.matches) {
  ese.matches <- strsplit(x = ese.matches,
                          split = ",")[[1]]
  any(Hexamers.Sorted.Into.Robustness.Groups[[ese.group]] %in% ese.matches)
}







mean(sapply(X = Bed.Table.Constitutive$ese.matches,
            FUN = function(y){
              EseGroupInSequence(ese.group = 1, ese.matches = y)
            }))

mean(sapply(X = Bed.Table.Alternative$ese.matches,
            FUN = function(y){
              EseGroupInSequence(ese.group = 1, ese.matches = y)
            }))


mean(sapply(X = Bed.Table.Constitutive$ese.matches,
            FUN = function(y){
              EseGroupInSequence(ese.group = 16, ese.matches = y)
            }))

mean(sapply(X = Bed.Table.Alternative$ese.matches,
            FUN = function(y){
              EseGroupInSequence(ese.group = 16, ese.matches = y)
            }))




Plotting.Matrix <- matrix(data = NA, nrow = 18, ncol = 3)

for (i in 1:18) {
  print(i)
  pcnt.constitutive.exons.with.hexamers.group.i <- mean(sapply(X = Bed.Table.Constitutive$ese.matches,
                                                               FUN = function(y){
                                                                 EseGroupInSequence(ese.group = i, ese.matches = y)
                                                               }))
  
  pcnt.alternative.exons.with.hexamers.group.i <- mean(sapply(X = Bed.Table.Alternative$ese.matches,
                                                              FUN = function(y){
                                                                EseGroupInSequence(ese.group = i, ese.matches = y)
                                                              }))
  
  Plotting.Matrix[i,] <- c(18 - i + 1,
                           pcnt.constitutive.exons.with.hexamers.group.i,
                           pcnt.alternative.exons.with.hexamers.group.i)
  
}



Plotting.Matrix <- as.data.frame(Plotting.Matrix)
names(Plotting.Matrix) <- c("Neighbourhood.ESEs",
                            "Constitutive",
                            "Alternative")

library(tidyr)

Plotting.Matrix <- gather(data = Plotting.Matrix,
                          key = Exon.Type,
                          value = Percent.Exons.With.Hexamers.This.Group,
                          Constitutive, Alternative)
Plotting.Matrix$Neighbourhood.ESEs <- factor(x = Plotting.Matrix$Neighbourhood.ESEs,
                                             levels = as.character(18:1))

library(ggplot2)

Plot <- ggplot(data = Plotting.Matrix,
               aes(x = Exon.Type, y = Percent.Exons.With.Hexamers.This.Group, fill = Exon.Type)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~Neighbourhood.ESEs, scale = "free", ncol = 6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values=c(rgb(255/255,
                                 255/255,
                                 255/255,
                                 0.25),
                             rgb(51/255,
                                 134/255,
                                 255/255,
                                 0.5))) +
  ylab("Percent of exons with at least one hexamer in group") +
  xlab("")
ggsave(filename = "239_Robustness_Barplots.pdf", plot = Plot, width = 15, height = 9)









# normalised by exon length

SumEseGroupInSequence <- function(ese.group, ese.matches) {
  ese.matches <- strsplit(x = ese.matches,
                          split = ",")[[1]]
  sum(Hexamers.Sorted.Into.Robustness.Groups[[ese.group]] %in% ese.matches)
}


Plotting.Matrix <- matrix(data = NA, nrow = 18, ncol = 3)

for (i in 1:18) {
  print(i)
  pcnt.constitutive.exons.with.hexamers.group.i <- sapply(X = Bed.Table.Constitutive$ese.matches,
                                                          FUN = function(y){
                                                            SumEseGroupInSequence(ese.group = i, ese.matches = y)
                                                          })
  pcnt.constitutive.exons.with.hexamers.group.i <- pcnt.constitutive.exons.with.hexamers.group.i/Bed.Table.Constitutive$length
  
  pcnt.constitutive.exons.with.hexamers.group.i <- mean(pcnt.constitutive.exons.with.hexamers.group.i)
  
  pcnt.alternative.exons.with.hexamers.group.i <- sapply(X = Bed.Table.Alternative$ese.matches,
                                                         FUN = function(y){
                                                           SumEseGroupInSequence(ese.group = i, ese.matches = y)
                                                         })
  
  pcnt.alternative.exons.with.hexamers.group.i <- pcnt.alternative.exons.with.hexamers.group.i/Bed.Table.Alternative$length
  
  pcnt.alternative.exons.with.hexamers.group.i <- mean(pcnt.alternative.exons.with.hexamers.group.i)
  
  Plotting.Matrix[i,] <- c(18 - i + 1,
                           pcnt.constitutive.exons.with.hexamers.group.i,
                           pcnt.alternative.exons.with.hexamers.group.i)
  
}



Plotting.Matrix <- as.data.frame(Plotting.Matrix)
names(Plotting.Matrix) <- c("Neighbourhood.ESEs",
                            "Constitutive",
                            "Alternative")

library(tidyr)

Plotting.Matrix <- gather(data = Plotting.Matrix,
                          key = Exon.Type,
                          value = Percent.Exons.With.Hexamers.This.Group,
                          Constitutive, Alternative)
Plotting.Matrix$Neighbourhood.ESEs <- factor(x = Plotting.Matrix$Neighbourhood.ESEs,
                                             levels = as.character(18:1))

library(ggplot2)

Plot <- ggplot(data = Plotting.Matrix,
               aes(x = Exon.Type, y = Percent.Exons.With.Hexamers.This.Group, fill = Exon.Type)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~Neighbourhood.ESEs, scale = "free", ncol = 6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values=c(rgb(255/255,
                                 255/255,
                                 255/255,
                                 0.25),
                             rgb(51/255,
                                 134/255,
                                 255/255,
                                 0.5))) +
  ylab("# sites per nucleotide") +
  xlab("")
ggsave(filename = "239_Robustness_Barplots_Normalised_By_Exon_Length.pdf", plot = Plot, width = 15, height = 9)
















# 

SumEseGroupInSequence <- function(ese.group, ese.matches) {
  ese.matches <- strsplit(x = ese.matches,
                          split = ",")[[1]]
  sum(Hexamers.Sorted.Into.Robustness.Groups[[ese.group]] %in% ese.matches)
}

# 
# Plotting.Matrix <- matrix(data = NA, nrow = 18, ncol = 5)
# 
# for (i in 1:18) {
#   print(i)
#   pcnt.constitutive.exons.with.hexamers.group.i <- sapply(X = Bed.Table.Constitutive$ese.matches,
#                                                           FUN = function(y){
#                                                             SumEseGroupInSequence(ese.group = i, ese.matches = y)
#                                                           })
#   pcnt.constitutive.exons.with.hexamers.group.i <- pcnt.constitutive.exons.with.hexamers.group.i/Bed.Table.Constitutive$length
#   
#   const.sd <- sd(pcnt.constitutive.exons.with.hexamers.group.i)
#   
#   pcnt.constitutive.exons.with.hexamers.group.i <- mean(pcnt.constitutive.exons.with.hexamers.group.i)
#   
#   pcnt.alternative.exons.with.hexamers.group.i <- sapply(X = Bed.Table.Alternative$ese.matches,
#                                                          FUN = function(y){
#                                                            SumEseGroupInSequence(ese.group = i, ese.matches = y)
#                                                          })
#   
#   pcnt.alternative.exons.with.hexamers.group.i <- pcnt.alternative.exons.with.hexamers.group.i/Bed.Table.Alternative$length
#   
#   alt.sd <- sd(pcnt.alternative.exons.with.hexamers.group.i)
#   
#   pcnt.alternative.exons.with.hexamers.group.i <- mean(pcnt.alternative.exons.with.hexamers.group.i)
#   
#   Plotting.Matrix[i,] <- c(18 - i + 1,
#                            pcnt.constitutive.exons.with.hexamers.group.i,
#                            pcnt.alternative.exons.with.hexamers.group.i,
#                            const.sd,
#                            alt.sd)
#   
# }
# 
# 
# Plotting.Matrix <- as.data.frame(Plotting.Matrix)
# names(Plotting.Matrix) <- c("Neighbourhood.ESEs",
#                             "Constitutive.Mean",
#                             "Alternative.Mean",
#                             "Constitutive.SD",
#                             "Alternative.SD")
# 
# Plotting.Matrix$Ratio <- Plotting.Matrix[,2]/Plotting.Matrix[,3]
# 
# Plotting.Matrix$Ratio.SD <- abs(Plotting.Matrix$Ratio) * sqrt((Plotting.Matrix$Constitutive.SD/Plotting.Matrix$Constitutive.Mean)^2 + (Plotting.Matrix$Alternative.SD/Plotting.Matrix$Alternative.Mean)^2)
# 
# 
# Plotting.Matrix$Ratio.Upper <- Plotting.Matrix$Ratio + 1.96 * Plotting.Matrix$Ratio.SD
# Plotting.Matrix$Ratio.Lower <- Plotting.Matrix$Ratio - 1.96 * Plotting.Matrix$Ratio.SD
# 
# 
# 
# library(epitools)
# 
# ?oddsratio
# oddsratio(x = Plotting.Matrix$Constitutive.Mean[1],
#           y = Plotting.Matrix$Alternative.Mean[1])
# 
# 
# 
# 
# 
# 
# Summary.Table <- Summary.Table/(apply(Summary.Table,1,sum))*100
# 
# Plot.DF <- as.data.frame(Summary.Table)
# 


















Plotting.Matrix <- matrix(data = NA, nrow = 18, ncol = 2)

for (i in 1:18) {
  print(i)
  pcnt.constitutive.exons.with.hexamers.group.i <- mean(sapply(X = Bed.Table.Constitutive$ese.matches,
                                                               FUN = function(y){
                                                                 EseGroupInSequence(ese.group = i, ese.matches = y)
                                                               }))
  
  pcnt.alternative.exons.with.hexamers.group.i <- mean(sapply(X = Bed.Table.Alternative$ese.matches,
                                                              FUN = function(y){
                                                                EseGroupInSequence(ese.group = i, ese.matches = y)
                                                              }))
  
  Plotting.Matrix[i,] <- c(18 - i + 1,
                           pcnt.constitutive.exons.with.hexamers.group.i/pcnt.alternative.exons.with.hexamers.group.i)
  
}

plot(Plotting.Matrix,
     type = "l",
     ylim = c(0,2),
     las = 1,
     xlab = "number of neighbouring ESEs",
     ylab = "% of constitutive exons / % of alternative exons",
     xaxt = "n",
     col = rgb(51/255,
               134/255,
               255/255,
               0.995),
     lwd = 3)
axis(side = 1, at = 1:18)




























Bed.Table.Alternative <- Bed.Table[Bed.Table$Exon.Class == "alternative",]
Bed.Table.Constitutive <- Bed.Table[Bed.Table$Exon.Class == "constitutive",]



Plotting.Matrix <- matrix(data = NA, nrow = 18, ncol = 2)

for (i in 1:18) {
  
  print(i)
  pcnt.constitutive.exons.with.hexamers.group.i <- sapply(X = Bed.Table.Constitutive$ese.matches,
                                                          FUN = function(y){
                                                            SumEseGroupInSequence(ese.group = i, ese.matches = y)
                                                          })
  pcnt.constitutive.exons.with.hexamers.group.i <- pcnt.constitutive.exons.with.hexamers.group.i/Bed.Table.Constitutive$length
  
  pcnt.constitutive.exons.with.hexamers.group.i <- mean(pcnt.constitutive.exons.with.hexamers.group.i)
  
  pcnt.alternative.exons.with.hexamers.group.i <- sapply(X = Bed.Table.Alternative$ese.matches,
                                                         FUN = function(y){
                                                           SumEseGroupInSequence(ese.group = i, ese.matches = y)
                                                         })
  
  pcnt.alternative.exons.with.hexamers.group.i <- pcnt.alternative.exons.with.hexamers.group.i/Bed.Table.Alternative$length
  
  pcnt.alternative.exons.with.hexamers.group.i <- mean(pcnt.alternative.exons.with.hexamers.group.i)
  
  
  
  Plotting.Matrix[i,] <- c(18 - i + 1,
                           pcnt.constitutive.exons.with.hexamers.group.i/pcnt.alternative.exons.with.hexamers.group.i)
  
}

plot(Plotting.Matrix,
     type = "l",
     ylim = c(0,2),
     las = 1,
     xlab = "number of neighbouring ESEs",
     ylab = "# in constitutive exons / # in alternative exons",
     xaxt = "n",
     col = rgb(51/255,
               134/255,
               255/255,
               0.995),
     lwd = 3)
axis(side = 1, at = 1:18)
abline(h=1)



Plot.DF <- as.data.frame(Plotting.Matrix)
names(Plot.DF) <- c("near.ESEs", "Ratio")

My.Plot <- ggplot(data = Plot.DF,
                  mapping = aes(x = near.ESEs, y = Ratio)) +
  geom_hline(yintercept = 1, colour = "gray40", linetype = "dashed") +
  geom_line(colour = "dodgerblue1", size = 1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        # axis.line.y = element_line(colour = "black"),
        aspect.ratio = 0.65,
        legend.position="none") +
  scale_x_continuous(breaks = 1:18,labels=c(1:18)) +
  ggtitle("") +
  ylab("ratio")+
  xlab("near ESEs")
library(Cairo)
ggsave(plot = My.Plot, filename = "239a_line_plot.pdf", width = 8, height = 8, device = cairo_pdf)






































load("238_Alternative_High_In_One_Tissue.RData")


Bed.Table$Exon.Class <- sapply(X = Bed.Table$id,
                               FUN = function(x){
                                 Constitutive.Alternative.Exon.DF[x,"Exon.Class"]
                               })

Bed.Table.Alternative <- Bed.Table[Bed.Table$Exon.Class == "alternative",]
Bed.Table.Constitutive <- Bed.Table[Bed.Table$Exon.Class == "constitutive",]




Plotting.Matrix <- matrix(data = NA, nrow = 18, ncol = 2)

for (i in 1:18) {
  
  print(i)
  pcnt.constitutive.exons.with.hexamers.group.i <- sapply(X = Bed.Table.Constitutive$ese.matches,
                                                          FUN = function(y){
                                                            SumEseGroupInSequence(ese.group = i, ese.matches = y)
                                                          })
  pcnt.constitutive.exons.with.hexamers.group.i <- pcnt.constitutive.exons.with.hexamers.group.i/Bed.Table.Constitutive$length
  
  pcnt.constitutive.exons.with.hexamers.group.i <- mean(pcnt.constitutive.exons.with.hexamers.group.i)
  
  pcnt.alternative.exons.with.hexamers.group.i <- sapply(X = Bed.Table.Alternative$ese.matches,
                                                         FUN = function(y){
                                                           SumEseGroupInSequence(ese.group = i, ese.matches = y)
                                                         })
  
  pcnt.alternative.exons.with.hexamers.group.i <- pcnt.alternative.exons.with.hexamers.group.i/Bed.Table.Alternative$length
  
  pcnt.alternative.exons.with.hexamers.group.i <- mean(pcnt.alternative.exons.with.hexamers.group.i)
  
  
  
  Plotting.Matrix[i,] <- c(18 - i + 1,
                           pcnt.constitutive.exons.with.hexamers.group.i/pcnt.alternative.exons.with.hexamers.group.i)
  
}

plot(Plotting.Matrix,
     type = "l",
     ylim = c(0,2),
     las = 1,
     xlab = "number of neighbouring ESEs",
     ylab = "# in constitutive exons / # in alternative exons",
     xaxt = "n",
     col = rgb(51/255,
               134/255,
               255/255,
               0.995),
     lwd = 3)
axis(side = 1, at = 1:18)
abline(h=1)

Plot.DF <- as.data.frame(Plotting.Matrix)
names(Plot.DF) <- c("near.ESEs", "Ratio")

My.Plot <- ggplot(data = Plot.DF,
       mapping = aes(x = near.ESEs, y = Ratio)) +
  geom_hline(yintercept = 1, colour = "gray40", linetype = "dashed") +
  geom_line(colour = "dodgerblue1", size = 1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        # axis.line.y = element_line(colour = "black"),
        aspect.ratio = 0.65,
        legend.position="none") +
  scale_x_continuous(breaks = 1:18,labels=c(1:18)) +
  ggtitle("alternative exons with a PSI > 90% in at least one tissue") +
  ylab("ratio")+
  xlab("near ESEs")

ggsave(plot = My.Plot, filename = "239b_Alt90_at_least_one_tissue.pdf", width = 8, height = 8, device = cairo_pdf)


#
























































































library(ggplot2)
Tissues <- c("Adipose Tissue",
             "Adrenal Gland",
             "Bladder",
             "Blood",
             "Blood Vessel",
             #"Bone Marrow",
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



Plot.List <- vector(mode = "list", length = length(Tissues))
names(Plot.List) <- Tissues

for (each.tissue in Tissues) {
  
  print(each.tissue)
  
  tissue.no.spaces <- gsub(pattern = " ",
                           replacement = "",
                           x = each.tissue,
                           fixed = TRUE)
  
  file.location <- paste("182_GTEx_PSI_Distributions_Plus_IDs/ALL_PSIandID_",
                         tissue.no.spaces,
                         ".RData",
                         sep = "")
  
  load(file.location)
  
  if (length(Mean.PSI.Values)==0){
    next
  }
  
  
  if (any(is.nan(Mean.PSI.Values))) {
    idx <- which(is.nan(Mean.PSI.Values))
    Mean.PSI.Values <- Mean.PSI.Values[-idx]
    Splicing.Events <- Splicing.Events[-idx]
  } 
  
  
  # keep everything above 90% inclusion only
  
  idx_above_90 <- which(Mean.PSI.Values > 0.9)
  if (length(idx_above_90)>0){
    Mean.PSI.Values <- Mean.PSI.Values[idx_above_90]
    Splicing.Events <- Splicing.Events[idx_above_90]
  }
  
  
  
  bed_table_idx <- which(Bed.Table$id %in% Splicing.Events)
  Bed.Table.This.Tissue <- Bed.Table[bed_table_idx,]
  
  
  Bed.Table.Alternative <- Bed.Table.This.Tissue[Bed.Table.This.Tissue$Exon.Class == "alternative",]
  Bed.Table.Constitutive <- Bed.Table.This.Tissue[Bed.Table.This.Tissue$Exon.Class == "constitutive",]
  
  
  
  
  Plotting.Matrix <- matrix(data = NA, nrow = 18, ncol = 2)
  
  for (i in 1:18) {
    
    print(i)
    pcnt.constitutive.exons.with.hexamers.group.i <- sapply(X = Bed.Table.Constitutive$ese.matches,
                                                            FUN = function(y){
                                                              SumEseGroupInSequence(ese.group = i, ese.matches = y)
                                                            })
    pcnt.constitutive.exons.with.hexamers.group.i <- pcnt.constitutive.exons.with.hexamers.group.i/Bed.Table.Constitutive$length
    
    pcnt.constitutive.exons.with.hexamers.group.i <- mean(pcnt.constitutive.exons.with.hexamers.group.i)
    
    pcnt.alternative.exons.with.hexamers.group.i <- sapply(X = Bed.Table.Alternative$ese.matches,
                                                           FUN = function(y){
                                                             SumEseGroupInSequence(ese.group = i, ese.matches = y)
                                                           })
    
    pcnt.alternative.exons.with.hexamers.group.i <- pcnt.alternative.exons.with.hexamers.group.i/Bed.Table.Alternative$length
    
    pcnt.alternative.exons.with.hexamers.group.i <- mean(pcnt.alternative.exons.with.hexamers.group.i)
    
    
    
    Plotting.Matrix[i,] <- c(18 - i + 1,
                             pcnt.constitutive.exons.with.hexamers.group.i/pcnt.alternative.exons.with.hexamers.group.i)
    
  }
  
  Plot.DF <- as.data.frame(Plotting.Matrix)
  names(Plot.DF) <- c("near.ESEs", "Ratio")
  
  My.Plot <- ggplot(data = Plot.DF,
                    mapping = aes(x = near.ESEs, y = Ratio)) +
    geom_hline(yintercept = 1, colour = "gray40", linetype = "dashed") +
    geom_line(colour = "dodgerblue1", size = 1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # panel.border = element_blank(),
          # axis.line.y = element_line(colour = "black"),
          aspect.ratio = 0.65,
          legend.position="none") +
    scale_x_continuous(breaks = 1:18,labels=c(1:18)) +
    ggtitle(each.tissue) +
    ylab("ratio")+
    xlab("near ESEs")
  
  
  
  Plot.List[[each.tissue]] <- My.Plot
  
  
  
}








library(gridExtra)
library(Cairo)
All.Plots <- do.call(what = "grid.arrange",
                         args = c(Plot.List,
                                  nrow = 6,
                                  ncol = 5,
                                  left = "number of sites per nucleotide in constitutive exons / number of sites per nucleotide in alternative exons",
                                  bottom = "ESE robustness"))

ggsave(filename = "239b_Plots_Robust_ESE_Per_Tissue_Highly_Included_Only.pdf",
       width = 2*8.27,
       height = 15,
       #useDingbats = F,
       plot = All.Plots,
       device=cairo_pdf)




