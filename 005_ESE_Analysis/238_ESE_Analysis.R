

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








# count how many (overlapping) ESEs we find, using
ESEseq <- readxl::read_excel(path = "238_Supplemental_Table_1.xls",
                             sheet = 2,
                             range = "A2:A1184",
                             col_names = TRUE)

# function to count number of overlapping matches
find_overlaps <- function(p,s) {
  gg <- gregexpr(paste0("(?=",p,")"),s,perl=TRUE)[[1]]
  if (length(gg)==1 && gg==-1) 0 else length(gg)
}


Bed.Table$ese.count <- sapply(X = Bed.Table$sequence,
                              FUN = function(mySequence){
                                sum(sapply(X = ESEseq$ESEseq,
                                           FUN = function(ese){
                                             find_overlaps(p = ese,
                                                           s = mySequence)
                                           }))
                              })


Bed.Table$strong.ese.count <- sapply(X = Bed.Table$sequence,
                                     FUN = function(mySequence){
                                       sum(sapply(X = ESEseq$ESEseq[1:200],
                                                  FUN = function(ese){
                                                    find_overlaps(p = ese,
                                                                  s = mySequence)
                                                  }))
                                     })


Bed.Table$medium.ese.count <- sapply(X = Bed.Table$sequence,
                                     FUN = function(mySequence){
                                       sum(sapply(X = ESEseq$ESEseq[1:500],
                                                  FUN = function(ese){
                                                    find_overlaps(p = ese,
                                                                  s = mySequence)
                                                  }))
                                     })


hist(Bed.Table$medium.ese.count/(Bed.Table$length/6), breaks = 40)


save(Bed.Table, file = "238_Bed_Table_ESE_Counts.RData")





library(Biostrings)

countPattern(pattern = "AA",
             subject = DNAString("GCATNAATGGG"), max.mismatch = 1,
             with.indels = F)


# this takes AGES to run
Bed.Table$fuzzy.ese.count <- sapply(X = Bed.Table$sequence,
                              FUN = function(mySequence){
                                sum(sapply(X = ESEseq$ESEseq,
                                           FUN = function(ese){
                                             countPattern(pattern = ese,
                                                          subject = DNAString(mySequence),
                                                          max.mismatch = 1,
                                                          with.indels = F)
                                           }))
                              })

save(Bed.Table, file = "238_Bed_Table_ESE_Counts.RData")

load("238_Bed_Table_ESE_Counts.RData")




































































load("238_Bed_Table_ESE_Counts.RData")


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
                                                       if ((x[2] > 0) & (x[3] >= 1)) { # in at least one tissue > 90%
                                                         result <- "alternative"
                                                       } else if (x[3] == x[1]) {
                                                         result <- "constitutive"
                                                       }
                                                       result
                                                     })
save(Constitutive.Alternative.Exon.DF, file = "238_Alternative_High_In_One_Tissue.RData")

table(Constitutive.Alternative.Exon.DF$Exon.Class)





Bed.Table$Exon.Class <- sapply(X = Bed.Table$id,
                               FUN = function(x){
                                 Constitutive.Alternative.Exon.DF[x,"Exon.Class"]
                               })


Bed.Table$ese.density <- Bed.Table$ese.count/Bed.Table$length


Bed.Table$ese.density.group <- findInterval(x = Bed.Table$ese.density*6,
                                                   vec = seq(0,3,0.3), all.inside = T)



ggplot(data = Bed.Table,
       mapping = aes(x = ese.density*6, colour = Exon.Class)) +
  geom_density()

t.test(x = (Bed.Table$ese.density*6)[which(Bed.Table$Exon.Class == "alternative")],
       y = (Bed.Table$ese.density*6)[which(Bed.Table$Exon.Class == "constitutive")])

wilcox.test(x = (Bed.Table$ese.density*6)[which(Bed.Table$Exon.Class == "alternative")],
            y = (Bed.Table$ese.density*6)[which(Bed.Table$Exon.Class == "constitutive")])

Summary.Table <- table(Bed.Table$ese.density.group, Bed.Table$Exon.Class)
Summary.Table <- Summary.Table[,-3]
Summary.Table <- Summary.Table/(apply(Summary.Table,1,sum))*100

Plot.DF <- as.data.frame(Summary.Table)


My.Plot <- ggplot(data = Plot.DF,
       mapping = aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(position="stack", stat="identity", aes(color = "black")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.text.x=element_text(angle=45,hjust=1),
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
  scale_x_discrete("ESE hexamers per nucleotide",
                   labels=c("x < 0.3",
                            "0.3 \u2264 x < 0.6",
                            "0.6 \u2264 x < 0.9",
                            "0.9 \u2264 x < 1.2",
                            "1.2 \u2264 x < 1.5",
                            "1.5 \u2264 x < 1.8",
                            "1.8 \u2264 x < 2.1",
                            "2.1 \u2264 x < 2.4",
                            "2.4 \u2264 x < 2.7",
                            "2.7 \u2264 x")) +
  ylab("percent")
library(Cairo)
ggsave(plot = My.Plot, filename = "238_ESE_Alt90_at_least_one_tissue.pdf", width = 8, height = 8, device = cairo_pdf)
#ggsave(plot = My.Plot, filename = "238_ESE.pdf", width = 8, height = 8)




Bed.Table$fuzzy.ese.density <- Bed.Table$fuzzy.ese.count/Bed.Table$length


Bed.Table$fuzzy.ese.density.group <- findInterval(x = Bed.Table$fuzzy.ese.density*6,
                                      vec = seq(0,60,6), all.inside = T)


Summary.Table <- table(Bed.Table$fuzzy.ese.density.group, Bed.Table$Exon.Class)
Summary.Table <- Summary.Table[,-3]
Summary.Table <- Summary.Table/(apply(Summary.Table,1,sum))*100

Plot.DF <- as.data.frame(Summary.Table)


My.Plot <- ggplot(data = Plot.DF,
       mapping = aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(position="stack", stat="identity", aes(color = "black")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.text.x=element_text(angle=45,hjust=1),
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
  scale_x_discrete("ESE hexamers per nucleotide (one Hamming distance away)",
                   labels=c("x < 6",
                            "6 \u2264 x < 12",
                            "12 \u2264 x < 18",
                            "18 \u2264 x < 24",
                            "24 \u2264 x < 30",
                            "30 \u2264 x < 36",
                            "36 \u2264 x < 42",
                            "42 \u2264 x < 48",
                            "48 \u2264 x < 54",
                            "54 \u2264 x")) +
  ylab("percent")
library(Cairo)
ggsave(plot = My.Plot, filename = "238_1HammAway_Alt90_at_least_one_tissue.pdf", width = 8, height = 8, device = cairo_pdf)
#ggsave(plot = My.Plot, filename = "238_1HammAway.pdf", width = 8, height = 8)




















































# global definition of alternative vs constitutive
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
  
  
  Summary.Table <- table(Bed.Table.This.Tissue$ese.density.group, Bed.Table.This.Tissue$Exon.Class)
  Summary.Table <- Summary.Table[,-3]
  Summary.Table <- Summary.Table/(apply(Summary.Table,1,sum))*100
  
  Plot.DF <- as.data.frame(Summary.Table)
  
  
  My.Plot <- ggplot(data = Plot.DF,
                    mapping = aes(x = Var1, y = Freq, fill = Var2)) +
    geom_bar(position="stack", stat="identity", aes(color = "black")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.y = element_line(colour = "black"),
          axis.text.x=element_text(angle=45,hjust=1),
          aspect.ratio = 0.65,
          legend.position="none") +
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
    scale_x_discrete(labels=c("x < 0.3",
                              "0.3 \u2264 x < 0.6",
                              "0.6 \u2264 x < 0.9",
                              "0.9 \u2264 x < 1.2",
                              "1.2 \u2264 x < 1.5",
                              "1.5 \u2264 x < 1.8",
                              "1.8 \u2264 x < 2.1",
                              "2.1 \u2264 x < 2.4",
                              "2.4 \u2264 x < 2.7",
                              "2.7 \u2264 x")) +
    ggtitle(each.tissue) +
    ylab("")+
    xlab("")
  
  
  
  Plot.List[[each.tissue]] <- My.Plot
  
  
  
}




library(gridExtra)
library(Cairo)
All.Bar.Plots <- do.call(what = "grid.arrange",
                         args = c(Plot.List,
                                  nrow = 6,
                                  ncol = 5,
                                  left = "percent",
                                  bottom = "ESE hexamers per nucleotide"))

ggsave(filename = "238b_Supplementary_Bar_Plots_ESE_Per_Tissue_Highly_Included_Only.pdf",
       width = 2*8.27,
       height = 17.5,
       #useDingbats = F,
       plot = All.Bar.Plots,
       device=cairo_pdf)
#


















































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
  
  
  Summary.Table <- table(Bed.Table.This.Tissue$fuzzy.ese.density.group, Bed.Table.This.Tissue$Exon.Class)
  Summary.Table <- Summary.Table[,-3]
  Summary.Table <- Summary.Table/(apply(Summary.Table,1,sum))*100
  
  Plot.DF <- as.data.frame(Summary.Table)
  
  
  My.Plot <- ggplot(data = Plot.DF,
                    mapping = aes(x = Var1, y = Freq, fill = Var2)) +
    geom_bar(position="stack", stat="identity", aes(color = "black")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.y = element_line(colour = "black"),
          axis.text.x=element_text(angle=45,hjust=1),
          aspect.ratio = 0.65,
          legend.position="none") +
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
    scale_x_discrete(labels=c("x < 6",
                              "6 \u2264 x < 12",
                              "12 \u2264 x < 18",
                              "18 \u2264 x < 24",
                              "24 \u2264 x < 30",
                              "30 \u2264 x < 36",
                              "36 \u2264 x < 42",
                              "42 \u2264 x < 48",
                              "48 \u2264 x < 54",
                              "54 \u2264 x")) +
    ggtitle(each.tissue) +
    ylab("")+
    xlab("")
  
  
  
  Plot.List[[each.tissue]] <- My.Plot
  
  
  
}




library(gridExtra)
library(Cairo)
All.Bar.Plots <- do.call(what = "grid.arrange",
                         args = c(Plot.List,
                                  nrow = 6,
                                  ncol = 5,
                                  left = "percent",
                                  bottom = "ESE hexamers per nucleotide (one Hamming distance away)"))

ggsave(filename = "238c_Supplementary_Bar_Plots_ESE_1HammDist_Per_Tissue_Highly_Included_Only.pdf",
       width = 2*8.27,
       height = 17.5,
       #useDingbats = F,
       plot = All.Bar.Plots,
       device=cairo_pdf)




