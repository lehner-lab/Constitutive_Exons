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
  
  file.location <- paste("./220_GTEx_Exon_Mutation_Effects_DoubleSlope/Also_ID_",
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



# negative mutations
Plot.DF <- data.frame(Starting.PSI = All.Starting.PSIs,
                      Delta.PSI = All.Mutation.Effects,
                      Final.PSI = All.Starting.PSIs + All.Mutation.Effects)
Plot.DF <- Plot.DF[-which(Plot.DF$Starting.PSI < 0),]
Plot.DF <- Plot.DF[-which(Plot.DF$Starting.PSI > 1),]




Plot.DF$Starting.PSI.Group <- EquallyPopulatedBins(vector = Plot.DF$Starting.PSI, k = 25)
table(Plot.DF$Starting.PSI.Group, useNA = "always")
Starting.PSI.Means <- sapply(1:max(Plot.DF$Starting.PSI.Group),
                             function(x){
                               mean(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group == x)])
                             })

Plot.DF$Starting.PSI.Group <- factor(Plot.DF$Starting.PSI.Group, levels = 1:length(Starting.PSI.Means))
sapply(1:25, function(x){ max(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group==x)])})

mean(abs(Plot.DF$Delta.PSI[which(Plot.DF$Starting.PSI.Group == 25)]) < 0.1)*100
median(abs(Plot.DF$Delta.PSI[which(Plot.DF$Starting.PSI.Group == 25)]))*100

mean(abs(Plot.DF$Delta.PSI[which(Plot.DF$Starting.PSI.Group == 13)]) < 0.1)*100
median(abs(Plot.DF$Delta.PSI[which(Plot.DF$Starting.PSI.Group == 13)]))*100



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
           label = "n = 4785-4786")

ggsave(filename = "220_Violin_Plots_GTEx_25_Groups_DoubleSlope.pdf", plot = Violin.Plots, width = 12, height = 12, useDingbats = F)






































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
                      
                      
                      # negative mutations
                      Plot.DF <- data.frame(Starting.PSI = Starting.PSIs,
                                            Delta.PSI = Mutation.Effects,
                                            Final.PSI = Starting.PSIs + Mutation.Effects)
                      Plot.DF <- Plot.DF[-which(Plot.DF$Starting.PSI < 0),]
                      Plot.DF <- Plot.DF[-which(Plot.DF$Starting.PSI > 1),]
                      
                      
                      Plot.DF$Starting.PSI.Group <- EquallyPopulatedBins(vector = Plot.DF$Starting.PSI, k = 25)
                      table(Plot.DF$Starting.PSI.Group, useNA = "always")
                      Starting.PSI.Means <- sapply(1:max(Plot.DF$Starting.PSI.Group),
                                                   function(x){
                                                     mean(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group == x)])
                                                   })
                      
                      Plot.DF$Starting.PSI.Group <- factor(Plot.DF$Starting.PSI.Group, levels = 1:length(Starting.PSI.Means))
                      sapply(1:25, function(x){ max(Plot.DF$Starting.PSI[which(Plot.DF$Starting.PSI.Group==x)])})
                      
                      
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
                              legend.position = "none",
                              axis.title = element_blank(),
                              axis.text.x = element_text(hjust = 1, angle = 45)) +
                        xlab("starting PSI") +
                        ylab("final PSI") +
                        ggtitle(each.tissue)+
                        scale_x_discrete(labels=c(as.character(round(Starting.PSI.Means, 4)*100)))# +
                        # annotate(geom = "text",
                        #          x = max(as.numeric(Plot.DF$Starting.PSI.Group)),
                        #          y = 0,
                        #          label = "n = 4785-4786")
                        
                      Violin.Plots
                      
                    })


library(gridExtra)
All.Violin.Plots <- do.call(what = "grid.arrange",
                            args = c(Plot.List,
                                     nrow = 10,
                                     ncol = 3,
                                     left = "final PSI",
                                     bottom = "starting PSI",
                                     top = "DPSI distributions at different starting PSI"))

ggsave(filename = "220b_Supplementary_Violin_Plots_Exons.pdf",
       width = 2*8.27,
       height = 2*11.69,
       useDingbats = F,
       plot = All.Violin.Plots)
