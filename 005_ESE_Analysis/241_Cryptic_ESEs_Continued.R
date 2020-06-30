




load("239_ESE_Robustness_Bed_Table.RData")

Alternative.Idx <- which(Bed.Table$Exon.Class == "alternative")

List.Alternative.Exons <- vector(mode = "list", length = length(Alternative.Idx))
names(List.Alternative.Exons) <- as.character(Alternative.Idx)

for (i in 1:246) {
  print(i)
  
  File.Name <- paste0("~/Desktop/mount2/project/prj004631/pbaeza/2020_Constitutive_Exons_Paper_Dimsum/005_Cryptic_ESEs/r_output/",
                      as.character(i),
                      ".RData")
  
  load(File.Name)
  
  
  if (any(names(delta.ese.counts) %in% as.character(Alternative.Idx))) {
    idx <- which(names(delta.ese.counts) %in% as.character(Alternative.Idx))
    
    for (each.idx in idx) {
      each.name <- names(delta.ese.counts)[each.idx]
      List.Alternative.Exons[[each.name]] <- delta.ese.counts[each.name]
    }
    
  }
  
}



hist(unlist(List.Alternative.Exons))
table(unlist(List.Alternative.Exons))
Alternative.Exon.Numbers <- as.numeric(table(unlist(List.Alternative.Exons)))
Alternative.Exon.Numbers.Pcnt <- Alternative.Exon.Numbers/(Constitutive.Exon.Numbers + Alternative.Exon.Numbers)






Constitutive.Idx <- which(Bed.Table$Exon.Class == "constitutive")

List.Constitutive.Exons <- vector(mode = "list", length = length(Constitutive.Idx))
names(List.Constitutive.Exons) <- as.character(Constitutive.Idx)

for (i in 1:246) {
  print(i)
  
  File.Name <- paste0("~/Desktop/mount2/project/prj004631/pbaeza/2020_Constitutive_Exons_Paper_Dimsum/005_Cryptic_ESEs/r_output/",
                      as.character(i),
                      ".RData")
  
  load(File.Name)
  
  
  if (any(names(delta.ese.counts) %in% as.character(Constitutive.Idx))) {
    idx <- which(names(delta.ese.counts) %in% as.character(Constitutive.Idx))
    
    for (each.idx in idx) {
      each.name <- names(delta.ese.counts)[each.idx]
      List.Constitutive.Exons[[each.name]] <- delta.ese.counts[each.name]
    }
    
  }
  
}


hist(unlist(List.Constitutive.Exons))
table(unlist(List.Constitutive.Exons))

Constitutive.Exon.Numbers <- as.numeric(table(unlist(List.Constitutive.Exons)))
Constitutive.Exon.Numbers.Pcnt <- Constitutive.Exon.Numbers/(Constitutive.Exon.Numbers + Alternative.Exon.Numbers)



Plot.DF <- data.frame(Sites.Created = rep(c(-6:4),2),
                      Exon.Class = c(rep("Alternative",11),
                                     rep("Constitutive",11)),
                      Frequency = c(Alternative.Exon.Numbers.Pcnt,
                                    Constitutive.Exon.Numbers.Pcnt))





library(ggplot2)

My.Plot <- ggplot(data = Plot.DF,
       mapping = aes(x = factor(Sites.Created, levels = as.character(-6:4)),
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
  xlab("number of ESEs lost/gained after destroying an existing ESE")

ggsave(plot = My.Plot, filename = "241_Cryptic_Sites.pdf", width = 8, height = 8)
#





Fold.Changes <- as.numeric(table(unlist(List.Constitutive.Exons))) / as.numeric(table(unlist(List.Alternative.Exons)))

plot(y = Fold.Changes,
     x = 1:11,
     type = "l",
     ylim = c(0,1.1))
abline(h=1)
