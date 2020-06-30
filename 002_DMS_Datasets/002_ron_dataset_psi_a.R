# library to read excel files
library(readxl)

# load data
Braun.Data <- read_excel(path = "Data/Braun_Supplementary_Table_3.xlsx",
                         sheet = 5)

# function that converts delta PSI values into delta splicing efficiency A
# (from equation 1 in figure 1)
FromDeltaPsiToA <- function(delta.psi, starting.psi){
  y <- delta.psi
  x <- starting.psi
  A <- (y + x - ((y*x + x^2)/(100)))/(x - ((y*x + x^2)/(100)))
  A
}

# the delta PSI depends on the initial (wild-type) PSI
WT.PSI <- Braun.Data$`AE inclusion (%)`[which(Braun.Data$mutation == "WT mean")]

# calculate A
Braun.Data$A <- FromDeltaPsiToA(delta.psi = Braun.Data$`AE inclusion (%)` - WT.PSI,
                                starting.psi = WT.PSI)


# save the position of each mutation
Braun.Data$Position <- sapply(X = as.character(Braun.Data$mutation),
                              FUN = function(x){
                                if (x == "WT mean" | x == "WT SD") {
                                  answer <- 0
                                } else{
                                  x.vector <- strsplit(x,"")[[1]]
                                  dna <- which(x.vector %in% c("A", "T", "G", "C"))
                                  answer <- as.numeric(paste0(x.vector[-dna], collapse = ""))
                                }
                                answer
                              })

# what exon/intron is the mutation in?
Braun.Data$Region <- sapply(X = Braun.Data$Position,
                            FUN = function(x){
                              if (x == 0) {
                                region <- "WT"
                              } else if (x < 211) {
                                region <- "exon.1"
                              } else if (x < 298) {
                                region <- "intron.1"
                              } else if (x < 445) {
                                region <- "alt.exon"
                              } else if (x < 525) {
                                region <- "intron.2"
                              } else if (x < 691) {
                                region <- "exon.3"
                              } else {
                                region <- "intron.3"
                              }
                              region
                            })

save(Braun.Data,
     file = "Data/Braun_et_al.RData")
