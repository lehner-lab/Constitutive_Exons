# PSI & A values for _RON_ exon 11 data

Here I explain the code in [002\_ron\_dataset\_psi\_a.R](./002_ron_dataset_psi_a.R), where I calculate the inclusion level (+ splicing efficiency values) for variants in the _RON_ exon 11 mutagenesis data from [Braun et al., 2018](http://dx.doi.org/10.1038/s41467-018-05748-7).

Unless stated otherwise, all of the code in this document is written in R.


## 1. Load data

The first step is to load the dataset from Supplementary Table 3 in [Braun et al., 2018](http://dx.doi.org/10.1038/s41467-018-05748-7):

```r
# library to read excel files
library(readxl)

# load data
Braun.Data <- read_excel(path = "Data/Braun_Supplementary_Table_3.xlsx",
                         sheet = 5)
```


## 2. PSI and A values

The PSI values for variants in this library are already given (they are shown in the `AE inclusion (%)` column), so the only thing left to be calculated is the effect 'A' of each mutation on splicing efficiency. So I first of all defined a function that converts the phenotypic effect of a mutation on exon inclusion (ΔPSI) into its biophysical effect on splicing efficiency (A):

```r
# function that converts delta PSI values into delta splicing efficiency A
# (from equation 1 in figure 1)
FromDeltaPsiToA <- function(delta.psi, starting.psi){
  y <- delta.psi
  x <- starting.psi
  A <- (y + x - ((y*x + x^2)/(100)))/(x - ((y*x + x^2)/(100)))
  A
}
```
With this I can now calculate the effect 'A' of each mutation:

```r
# the delta PSI depends on the initial (wild-type) PSI
WT.PSI <- Braun.Data$`AE inclusion (%)`[which(Braun.Data$mutation == "WT mean")]

# calculate A
Braun.Data$A <- FromDeltaPsiToA(delta.psi = Braun.Data$`AE inclusion (%)` - WT.PSI,
                                starting.psi = WT.PSI)
```

## 3. Some extra information for later

I'll also save some extra information about the position of each mutation relative to the exon, which might come in handy later:

```r
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
```


## 4. Save

Save for use later.

```r
# save
save(Braun.Data,
     file = "Data/Braun_et_al.RData")
```
