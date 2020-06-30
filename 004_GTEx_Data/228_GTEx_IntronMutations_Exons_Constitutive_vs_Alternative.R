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














for (each.tissue in Tissues) {
  
  print(each.tissue)
  
  tissue.no.spaces <- gsub(pattern = " ",
                           replacement = "",
                           x = each.tissue,
                           fixed = TRUE)
  
  
  file.location <- paste("./227_GTEx_Intron_Mutation_Effects_DoubleSlope/Also_ID_",
                         tissue.no.spaces,
                         ".RData",
                         sep = "")
  
  load(file.location)
  
  if (length(Starting.PSIs)==0){
    next
  }
  
  
  
  Plot.DF <- data.frame(Exon.ID = Exon.IDs,
                        Mutation.Effect = Mutation.Effects * 100 * 2, # ifx this when I get proper mutation effects
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
    scale_colour_manual(values = c("firebrick1", "gray30")) +
    geom_smooth(level=0.95, method = "loess") +
    coord_cartesian(ylim = c(-100,0))+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 1,
          legend.position = "bottom") +
    xlab("starting PSI") +
    ylab("delta PSI") +
    ggtitle(each.tissue) +
    annotate(geom = "text",
             x = 90,
             y = -100,
             label = paste("p = ",
                           as.character(P.Value),
                           sep = ""),
             hjust = 0,
             vjust = 0)
  
  File.Name <- paste("228_GTEx_IntronMutations_Exons_Constitutive_vs_Alternative/228_Scatter_Plot_ConstitutiveVsAlternative_",
                     tissue.no.spaces,
                     ".pdf",
                     sep = "",
                     collapse = "")
  ggsave(filename = File.Name, plot = My.Plot, width = 12, height = 12, useDingbats = F)
  
  
  
  
}















# Now do all exons in all tissues combined
Combined.Plot.DF <- data.frame(Exon.ID = c(),
                               Mutation.Effect = c(),
                               Starting.PSI = c(),
                               Exon.Class = c())

for (each.tissue in Tissues) {
  
  print(each.tissue)
  
  tissue.no.spaces <- gsub(pattern = " ",
                           replacement = "",
                           x = each.tissue,
                           fixed = TRUE)
  
  
  file.location <- paste("./227_GTEx_Intron_Mutation_Effects_DoubleSlope/Also_ID_",
                         tissue.no.spaces,
                         ".RData",
                         sep = "")
  
  load(file.location)
  
  if (length(Starting.PSIs)==0){
    next
  }
  
  
  
  Plot.DF <- data.frame(Exon.ID = Exon.IDs,
                        Mutation.Effect = Mutation.Effects * 100 * 2, # fix this when I get the proper DoubleSlope data
                        Starting.PSI = Starting.PSIs * 100,
                        Exon.Class = Constitutive.Alternative.Exon.DF[Exon.IDs,"Exon.Class"])
  
  
  Plot.DF <- Plot.DF[which(Plot.DF$Mutation.Effect<0),]
  Plot.DF <- Plot.DF[which(Plot.DF$Starting.PSI>=90),]
  Plot.DF <- Plot.DF[which(Plot.DF$Starting.PSI<=100),]
  
  
  Combined.Plot.DF <- rbind(Combined.Plot.DF,
                            Plot.DF)
  
}








P.Value <- signif(x = summary(lm(Mutation.Effect ~ Starting.PSI + Exon.Class,
                                 data = Combined.Plot.DF))$coefficients[3,4],
                  digits = 2)

My.Plot <- ggplot(data = Combined.Plot.DF,
                  mapping = aes(x = Starting.PSI,
                                y = Mutation.Effect,
                                colour = Exon.Class)) +
  geom_point(size = 0.01) +
  #stat_density_2d(geom = "polygon",aes(alpha = log2(..level..), fill = Exon.Class),bins = 4) +
  #geom_raster()+
  scale_colour_manual(values = c("firebrick1", "gray30")) +
  geom_smooth(level=0.95, method = "loess") +
  coord_cartesian(ylim = c(-100,0))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        legend.position = "bottom") +
  xlab("starting PSI") +
  ylab("delta PSI") +
  ggtitle("all tissues") +
  annotate(geom = "text",
           x = 90,
           y = -100,
           label = paste("p = ",
                         as.character(P.Value),
                         sep = ""),
           hjust = 0,
           vjust = 0)
ggsave(filename = "228_Pooled_Data_IntronicMutations_Exons_ConstitutiveVsAlternative.pdf", plot = My.Plot, width = 12, height = 12, useDingbats = F)
