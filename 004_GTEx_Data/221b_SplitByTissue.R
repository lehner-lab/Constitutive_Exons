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
                                                   ncol = 2,
                                                   dimnames = list(Splicing.Events,
                                                                   c("Number.Of.Tissues",
                                                                     "Times.Under.90"))
    )
    
    First.Iteration <- FALSE
  }
  
  
  if (any(is.nan(Mean.PSI.Values))) {
    idx <- which(is.nan(Mean.PSI.Values))
    Mean.PSI.Values <- Mean.PSI.Values[-idx]
    Splicing.Events <- Splicing.Events[-idx]
  } 
  
  Constitutive.Alternative.Exon.Matrix[Splicing.Events,"Number.Of.Tissues"] <- Constitutive.Alternative.Exon.Matrix[Splicing.Events,"Number.Of.Tissues"] + 1
  Constitutive.Alternative.Exon.Matrix[Splicing.Events[which(Mean.PSI.Values<0.9)],"Times.Under.90"] <- Constitutive.Alternative.Exon.Matrix[Splicing.Events[which(Mean.PSI.Values<0.9)],"Times.Under.90"] + 1
  
}


Constitutive.Alternative.Exon.DF <- as.data.frame(Constitutive.Alternative.Exon.Matrix)
Constitutive.Alternative.Exon.DF$Exon.Class <- sapply(X = Constitutive.Alternative.Exon.Matrix[,"Times.Under.90"],
                                                      FUN = function(x){
                                                        result <- "alternative"
                                                        if (x == 0) {
                                                          result <- "constitutive"
                                                        }
                                                        result
                                                      })







Tissues <- Tissues[-which(Tissues=="Bone Marrow")]


Plot.List <- lapply(X = Tissues,
                    FUN = function(each.tissue){
                      print(each.tissue)
                      
                      tissue.no.spaces <- gsub(pattern = " ",
                                               replacement = "",
                                               x = each.tissue,
                                               fixed = TRUE)
                      
                      
                      file.location <- paste("./220_GTEx_Exon_Mutation_Effects_DoubleSlope/Also_ID_",
                                             tissue.no.spaces,
                                             ".RData",
                                             sep = "")
                      
                      load(file.location)
                      
                      
                      
                      
                      Plot.DF <- data.frame(Exon.ID = Exon.IDs,
                                            Mutation.Effect = Mutation.Effects * 100,
                                            Starting.PSI = Starting.PSIs * 100,
                                            Exon.Class = Constitutive.Alternative.Exon.DF[Exon.IDs,"Exon.Class"])
                      
                      
                      Plot.DF <- Plot.DF[which(Plot.DF$Mutation.Effect<0),]
                      Plot.DF <- Plot.DF[which(Plot.DF$Starting.PSI>=90),]
                      Plot.DF <- Plot.DF[which(Plot.DF$Starting.PSI<=100),]
                      
                      
                      P.Value <- signif(x = summary(lm(Mutation.Effect ~ Starting.PSI + Exon.Class,
                                                       data = Plot.DF))$coefficients[3,4],
                                        digits = 2)
                      
                      My.Plot <- ggplot(data = Plot.DF,
                                        mapping = aes(x = Starting.PSI,
                                                      y = Mutation.Effect,
                                                      colour = Exon.Class)) +
                        geom_point() +
                        scale_colour_manual(values = c("dodgerblue1", "gray30")) +
                        geom_smooth(level=0.95, method = "loess") +
                        coord_cartesian(ylim = c(-100,0))+
                        theme_bw() +
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              legend.position = "none",
                              axis.title = element_blank()) +
                        ggtitle(each.tissue) +
                        annotate(geom = "text",
                                 x = 90,
                                 y = -100,
                                 label = paste("p = ",
                                               as.character(P.Value),
                                               sep = ""),
                                 hjust = 0,
                                 vjust = 0)
                      
                      My.Plot
                      
                    })


library(gridExtra)
All.Plots <- do.call(what = "grid.arrange",
                            args = c(Plot.List,
                                     nrow = 6,
                                     ncol = 5,
                                     left = "delta PSI",
                                     bottom = "starting PSI",
                                     top = "Alternative vs Constitutive exons"))

ggsave(filename = "221b_Split_By_Tissues.pdf",
       width = 15,
       height = 15,
       useDingbats = F,
       plot = All.Plots)
