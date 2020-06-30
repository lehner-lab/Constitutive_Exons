


# function to count number of overlapping matches
find_overlaps <- function(p,s) {
  gg <- gregexpr(paste0("(?=",p,")"),s,perl=TRUE)[[1]]
  if (length(gg)==1 && gg==-1) 0 else length(gg)
}
# function to count number of overlapping matches
find_overlaps_position <- function(p,s) {
  gg <- gregexpr(paste0("(?=",p,")"),s,perl=TRUE)[[1]]
  as.numeric(gg)
}



load("239_ESE_Robustness_Bed_Table.RData")
load("239bis_Hexamer_Table.RData")



Bed.Table$ese.matches.positions <- sapply(X = Bed.Table$sequence,
                                          FUN = function(mySequence){
                                            matches <- sapply(X = ESE$Hexamer,
                                                              FUN = function(ese){
                                                                find_overlaps_position(p = ese,
                                                                                       s = mySequence)
                                                              })
                                            
                                            matches <- unlist(matches)
                                            
                                            idx.to.remove <- which(matches < 1)
                                            
                                            if (length(idx.to.remove)>0){
                                              matches <- matches[-idx.to.remove]
                                            }
                                            
                                            if (length(matches) > 0) {
                                              matches <- sort(matches)
                                            }
                                            
                                            matches
                                          })


CalculateDistances <- function(myPositionsVector){
  vector.length <- length(myPositionsVector)
  
  distances <- c()
  
  for (i in 1:(vector.length-1)) {
    d <- myPositionsVector[i+1] - myPositionsVector[i]
    distances <- c(distances,d)
  }
  
  distances
  
}





Bed.Table$ese.matches.distances <- sapply(X = Bed.Table$ese.matches.positions,
                                          FUN = function(x){
                                            if (length(x)>1){
                                              CalculateDistances(x)
                                            } else {
                                              NA
                                            }
                                          })







Bed.Table.Alternative <- Bed.Table[which(Bed.Table$Exon.Class=="alternative"),]
Bed.Table.Constitutive <- Bed.Table[which(Bed.Table$Exon.Class=="constitutive"),]




Distances.Alternative <- unlist(Bed.Table.Alternative$ese.matches.distances)
Distances.Alternative <- Distances.Alternative[-which(is.na(Distances.Alternative))]
Distances.Constitutive <- unlist(Bed.Table.Constitutive$ese.matches.distances)
Distances.Constitutive <- Distances.Constitutive[-which(is.na(Distances.Constitutive))]




Distances.Alternative[which(Distances.Alternative>30)] <- 30
Distances.Constitutive[which(Distances.Constitutive>30)] <- 30




Alternative.Exon.Numbers <- as.numeric(table(Distances.Alternative))
Constitutive.Exon.Numbers <- as.numeric(table(Distances.Constitutive))

Alternative.Exon.Numbers.Pcnt <- Alternative.Exon.Numbers/(Constitutive.Exon.Numbers + Alternative.Exon.Numbers)
Constitutive.Exon.Numbers.Pcnt <- Constitutive.Exon.Numbers/(Constitutive.Exon.Numbers + Alternative.Exon.Numbers)


Plot.DF <- data.frame(Sites.Created = rep(c(1:30),2),
                      Exon.Class = c(rep("Alternative",30),
                                     rep("Constitutive",30)),
                      Frequency = c(Alternative.Exon.Numbers.Pcnt,
                                    Constitutive.Exon.Numbers.Pcnt))




library(ggplot2)

My.Plot <- ggplot(data = Plot.DF,
                  mapping = aes(x = factor(Sites.Created, levels = as.character(1:30)),
                                y = Frequency*100,
                                fill = Exon.Class)) +
  geom_bar(position="stack", stat="identity", aes(color = "black")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(colour = "black"),
        #axis.text.x=element_text(angle=45,hjust=1),
        aspect.ratio = 0.65) +
  scale_fill_manual(values=c(rgb(255/255,
                                 255/255,
                                 255/255,
                                 0.25),
                             rgb(51/255,
                                 134/255,
                                 255/255,
                                 0.5))) +
  scale_color_manual(values = "black") +
  scale_y_continuous(breaks = seq(0, 100, 25),
                     limits = c(-0.1, 100.1),
                     expand = c(0, 0)) +
  ylab("percent") +
  xlab("nucleotide distance between consecutive ESEs")

ggsave(plot = My.Plot, filename = "242_ESE_Distances.pdf", width = 8, height = 8)

























































Distances.Alternative <- unlist(Bed.Table.Alternative$ese.matches.distances)
Distances.Alternative <- Distances.Alternative[-which(is.na(Distances.Alternative))]
Distances.Constitutive <- unlist(Bed.Table.Constitutive$ese.matches.distances)
Distances.Constitutive <- Distances.Constitutive[-which(is.na(Distances.Constitutive))]



Alternative.Exon.DistanceGroup <- sapply(X = Distances.Alternative,
                                         FUN = function(x){
                                           if (x <11){
                                             "overlapping"
                                           } else {
                                             "not_overlapping"
                                           }
                                         })

Constitutive.Exon.DistanceGroup <- sapply(X = Distances.Constitutive,
                                         FUN = function(x){
                                           if (x <11){
                                             "overlapping"
                                           } else {
                                             "not_overlapping"
                                           }
                                         })


Alternative.Exon.Numbers <- as.numeric(table(Alternative.Exon.DistanceGroup))
Constitutive.Exon.Numbers <- as.numeric(table(Constitutive.Exon.DistanceGroup))
Alternative.Exon.Numbers.Pcnt <- Alternative.Exon.Numbers/(Constitutive.Exon.Numbers + Alternative.Exon.Numbers)
Constitutive.Exon.Numbers.Pcnt <- Constitutive.Exon.Numbers/(Constitutive.Exon.Numbers + Alternative.Exon.Numbers)


Plot.DF <- data.frame(Distance.Group = rep(c("not_overlapping",
                                             "overlapping"),2),
                      Exon.Class = c(rep("Alternative",2),
                                     rep("Constitutive",2)),
                      Frequency = c(Alternative.Exon.Numbers.Pcnt,
                                    Constitutive.Exon.Numbers.Pcnt))




library(ggplot2)

My.Plot <- ggplot(data = Plot.DF,
                  mapping = aes(x = Distance.Group,
                                y = Frequency*100,
                                fill = Exon.Class)) +
  geom_bar(position="stack", stat="identity", aes(color = "black")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(colour = "black"),
        #axis.text.x=element_text(angle=45,hjust=1),
        aspect.ratio = 1.5) +
  scale_fill_manual(values=c(rgb(255/255,
                                 255/255,
                                 255/255,
                                 0.25),
                             rgb(51/255,
                                 134/255,
                                 255/255,
                                 0.5))) +
  scale_color_manual(values = "black") +
  scale_y_continuous(breaks = seq(0, 100, 25),
                     limits = c(-0.1, 100.1),
                     expand = c(0, 0)) +
  ylab("percent") +
  xlab("")


ggsave(plot = My.Plot, filename = "242_ESE_OverlappingGroup_10nt.pdf", width = 8, height = 8)

