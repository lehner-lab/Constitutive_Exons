



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



Parameter.DF <- data.frame(Constitutive.Exon.Threshold = c(),
                           Delta.Effect = c())

for (n in seq(0.50,0.99,0.01)) {
  
  print(n)
  
  First.Iteration <- TRUE
  
  for (each.tissue in Tissues) {
    
    #print(each.tissue)
    
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
                                                                       "Times.Under.85"))
      )
      
      First.Iteration <- FALSE
    }
    
    
    if (any(is.nan(Mean.PSI.Values))) {
      idx <- which(is.nan(Mean.PSI.Values))
      Mean.PSI.Values <- Mean.PSI.Values[-idx]
      Splicing.Events <- Splicing.Events[-idx]
    } 
    
    Constitutive.Alternative.Exon.Matrix[Splicing.Events,"Number.Of.Tissues"] <- Constitutive.Alternative.Exon.Matrix[Splicing.Events,"Number.Of.Tissues"] + 1
    Constitutive.Alternative.Exon.Matrix[Splicing.Events[which(Mean.PSI.Values<n)],"Times.Under.85"] <- Constitutive.Alternative.Exon.Matrix[Splicing.Events[which(Mean.PSI.Values<n)],"Times.Under.85"] + 1
    
  }
  
  
  Constitutive.Alternative.Exon.DF <- as.data.frame(Constitutive.Alternative.Exon.Matrix)
  Constitutive.Alternative.Exon.DF$Exon.Class <- sapply(X = Constitutive.Alternative.Exon.Matrix[,c("Times.Under.85")],
                                                        FUN = function(x){
                                                          
                                                          result <- "alternative"
                                                          
                                                          if (x == 0) {
                                                            result <- "constitutive"
                                                          }
                                                          result
                                                        })
  
  
  
  
  
  
  
  
  
  
  
  
  # Now do all exons in all tissues combined
  Combined.Plot.DF <- data.frame(Exon.ID = c(),
                                 Mutation.Effect = c(),
                                 Starting.PSI = c(),
                                 Exon.Class = c())
  
  for (each.tissue in Tissues) {
    
    #print(each.tissue)
    
    tissue.no.spaces <- gsub(pattern = " ",
                             replacement = "",
                             x = each.tissue,
                             fixed = TRUE)
    
    
    file.location <- paste("./220_GTEx_Exon_Mutation_Effects_DoubleSlope/Also_ID_",
                           tissue.no.spaces,
                           ".RData",
                           sep = "")
    
    load(file.location)
    
    if (length(Starting.PSIs)==0){
      next
    }
    
    
    
    Plot.DF <- data.frame(Exon.ID = Exon.IDs,
                          Mutation.Effect = Mutation.Effects * 100,
                          Starting.PSI = Starting.PSIs * 100,
                          Exon.Class = Constitutive.Alternative.Exon.DF[Exon.IDs,"Exon.Class"])
    
    
    Plot.DF <- Plot.DF[which(Plot.DF$Mutation.Effect<0),]
    Plot.DF <- Plot.DF[which(Plot.DF$Starting.PSI>=(n*100)),]
    Plot.DF <- Plot.DF[which(Plot.DF$Starting.PSI<=100),]
    
    
    Combined.Plot.DF <- rbind(Combined.Plot.DF,
                              Plot.DF)
    
  }
  
  
  
  
  Value <- summary(lm(Mutation.Effect ~ Starting.PSI + Exon.Class,
                      data = Combined.Plot.DF))$coefficients[3,1]
  
  
  
  
  Temp.DF <- data.frame(Constitutive.Exon.Threshold = n*100,
                             Delta.Effect = Value)
  
  
  Parameter.DF <- rbind(Parameter.DF,
                        Temp.DF)
}



plot(Parameter.DF[which(Parameter.DF$Constitutive.Exon.Threshold>60),])

save(Parameter.DF, file = "225_Simultaneously_Change_Constitutive_and_Alternative_Exon_Threshold.RData")

load("225_Simultaneously_Change_Constitutive_and_Alternative_Exon_Threshold.RData")

Parameter.DF <- Parameter.DF[which(Parameter.DF$Constitutive.Exon.Threshold>60),]


My.Plot <- ggplot(data = Parameter.DF,
       mapping = aes(x = Constitutive.Exon.Threshold,
                     y = Delta.Effect)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



ggsave(filename = "225_Simultaneous_Parameter_Shift.pdf", plot = My.Plot, width = 8, height = 8, useDingbats = F)







#


















































































































































































































Parameter.DF <- data.frame(Constitutive.Exon.Threshold = c(),
                           Alternative.Exon.Threshold = c(),
                           Delta.Effect = c())

for (n in seq(0.60,0.99,0.01)) {
  
  print(n)
  
  
  for (m in seq(0.60,0.99,0.01)) {
    First.Iteration <- TRUE
    
    for (each.tissue in Tissues) {
      
      #print(each.tissue)
      
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
                                                                         "Times.Under.Constitutive.Threshold",
                                                                         "Times.Under.Alternative.Threshold"))
        )
        
        First.Iteration <- FALSE
      }
      
      
      if (any(is.nan(Mean.PSI.Values))) {
        idx <- which(is.nan(Mean.PSI.Values))
        Mean.PSI.Values <- Mean.PSI.Values[-idx]
        Splicing.Events <- Splicing.Events[-idx]
      } 
      
      Constitutive.Alternative.Exon.Matrix[Splicing.Events,"Number.Of.Tissues"] <- Constitutive.Alternative.Exon.Matrix[Splicing.Events,"Number.Of.Tissues"] + 1
      Constitutive.Alternative.Exon.Matrix[Splicing.Events[which(Mean.PSI.Values<n)],"Times.Under.Constitutive.Threshold"] <- Constitutive.Alternative.Exon.Matrix[Splicing.Events[which(Mean.PSI.Values<n)],"Times.Under.Constitutive.Threshold"] + 1
      
      
      Constitutive.Alternative.Exon.Matrix[Splicing.Events[which(Mean.PSI.Values<m)],"Times.Under.Alternative.Threshold"] <- Constitutive.Alternative.Exon.Matrix[Splicing.Events[which(Mean.PSI.Values<m)],"Times.Under.Alternative.Threshold"] + 1
      
    }
    
    
    Constitutive.Alternative.Exon.DF <- as.data.frame(Constitutive.Alternative.Exon.Matrix)
    Constitutive.Alternative.Exon.DF$Exon.Class <- apply(X = Constitutive.Alternative.Exon.Matrix[,c("Times.Under.Constitutive.Threshold", "Times.Under.Alternative.Threshold")],
                                                         MARGIN = 1,
                                                         FUN = function(x){
                                                           
                                                           result <- "other"
                                                           
                                                           if (x[1] == 0) {
                                                             result <- "constitutive"
                                                           }
                                                           
                                                           if (x[2] >= 1) {
                                                             result <- "alternative"
                                                           }
                                                           result
                                                         })
    
    
    Constitutive.Alternative.Exon.DF$Exon.Class <- as.character(Constitutive.Alternative.Exon.DF$Exon.Class)
    
    if (any(Constitutive.Alternative.Exon.DF$Exon.Class == "other")) {
      Constitutive.Alternative.Exon.DF <- Constitutive.Alternative.Exon.DF[-which(Constitutive.Alternative.Exon.DF$Exon.Class=="other"),]
    }
    
    
    
    
    
    
    
    
    # Now do all exons in all tissues combined
    Combined.Plot.DF <- data.frame(Exon.ID = c(),
                                   Mutation.Effect = c(),
                                   Starting.PSI = c(),
                                   Exon.Class = c())
    
    for (each.tissue in Tissues) {
      
      #print(each.tissue)
      
      tissue.no.spaces <- gsub(pattern = " ",
                               replacement = "",
                               x = each.tissue,
                               fixed = TRUE)
      
      
      file.location <- paste("./220_GTEx_Exon_Mutation_Effects_DoubleSlope/Also_ID_",
                             tissue.no.spaces,
                             ".RData",
                             sep = "")
      
      load(file.location)
      
      if (length(Starting.PSIs)==0){
        next
      }
      
      
      
      Plot.DF <- data.frame(Exon.ID = Exon.IDs,
                            Mutation.Effect = Mutation.Effects * 100,
                            Starting.PSI = Starting.PSIs * 100,
                            Exon.Class = Constitutive.Alternative.Exon.DF[Exon.IDs,"Exon.Class"])
      
      
      Plot.DF <- Plot.DF[which(Plot.DF$Mutation.Effect<0),]
      Plot.DF <- Plot.DF[which(Plot.DF$Starting.PSI>=(n*100)),]
      Plot.DF <- Plot.DF[which(Plot.DF$Starting.PSI<=100),]
      
      
      Combined.Plot.DF <- rbind(Combined.Plot.DF,
                                Plot.DF)
      
    }
    
    
    
    
    Value <- summary(lm(Mutation.Effect ~ Starting.PSI + Exon.Class,
                        data = Combined.Plot.DF))$coefficients[3,1]
    
    
    
    
    Temp.DF <- data.frame(Constitutive.Exon.Threshold = n*100,
                          Alternative.Exon.Threshold = m*100,
                          Delta.Effect = Value)
    
    
    Parameter.DF <- rbind(Parameter.DF,
                          Temp.DF) 
  }
  
}



library(metR) # for geom_contour_fill

save(Parameter.DF, file = "225_Independently_Change_Constitutive_and_Alternative_Exon_Threshold.RData")


# Triangle.Parameter.DF <- Parameter.DF
# Triangle.Parameter.DF$Delta.Effect[which(Triangle.Parameter.DF$Constitutive.Exon.Threshold <= Triangle.Parameter.DF$Alternative.Exon.Threshold)] <- 0

ggplot(data = Parameter.DF,
       mapping = aes(x = Constitutive.Exon.Threshold,
                     y = Alternative.Exon.Threshold,
                     z = Delta.Effect)) +
  geom_contour_fill(na.fill = T)+
  geom_contour(color = "black")



Loess.Model <- loess(Delta.Effect ~ Constitutive.Exon.Threshold + Alternative.Exon.Threshold,
                     data = Parameter.DF,
                     span = 0.75)

Parameter.DF$Loess.Smoother <- predict(object = Loess.Model,
                                       newdata = data.frame(Constitutive.Exon.Threshold = Parameter.DF$Constitutive.Exon.Threshold,
                                                            Alternative.Exon.Threshold = Parameter.DF$Alternative.Exon.Threshold))


My.Plot <- ggplot(data = Parameter.DF,
       mapping = aes(x = Constitutive.Exon.Threshold,
                     y = Alternative.Exon.Threshold,
                     z = Loess.Smoother)) +
  geom_contour_fill(na.fill = T)+
  geom_contour(color = "black") +
  theme_bw() +
  scale_x_continuous(limits = c(60,99), expand = c(0, 0)) +
  scale_y_continuous(limits = c(60,99), expand = c(0, 0)) +
  theme(aspect.ratio = 1) +
  scale_fill_gradient(low = "white", high = "dodgerblue1")
  scale_fill_manual("Z",values=colorRampPalette(c("white","dodgerblue1"))(n=12))


ggsave(filename = "225_2D_Parameter_Plot.pdf", plot = My.Plot, width = 8, height = 8)
