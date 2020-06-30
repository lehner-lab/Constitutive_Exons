# PSI & A values for human _FAS_ exon 6 data

Here I explain the code in [002\_fas\_human\_dataset\_psi\_a.R](./002_fas_human_dataset_psi_a.R), where I calculate the inclusion level (+ splicing efficiency values) for variants in the _FAS_ exon 6 mutagenesis data from [Julien et al., 2018](http://dx.doi.org/10.1038/ncomms11558).

Unless stated otherwise, all of the code in this document is written in R.


## 1. Load data

The first step is to load the FAS exon 6 mutagenesis dataset as processed by in the section [001\_Processing\_Sequencing\_Data](../001_Processing_Sequencing_Data):

```r
# load data
load("Data/FAS_Human_fitness_replicates.RData")
```


## 2. From enrichment scores to PSI values

Note that DiMSum (the pipeline we used to process this dataset earlier) already calculates enrichment scores (and returns them in the column labelled 'fitness'). Therefore, we can easily calculate the PSI of a given variant using the following equation:

<p align="center">
  <img src="Equations/psi_var_1.gif">
</p>

and because the normalised ES of the wild-type sequence is, by definition, 1, that equation simplifies to:

<p align="center">
  <img src="Equations/psi_var_2.gif">
</p>

```r
# PSI of the wild type (determined in Julien et al 2016)
WT.PSI <- 49.1

# ES of the WT is 1 by definition
WT.ES <- 1

# calculate PSI (note that DiMSum returns log ES so need to exp
singles$PSI <- exp(singles$fitness) * WT.PSI
```


## 3. Corrected PSI estimates

Since DiMSum also returns an estimate of the error associated with a given enrichment score, we can use this in combination with the method described in [Li et al., 2019](http://dx.doi.org/10.1038/s41467-019-11735-3) to correct any PSI values that are estimated to be above 100%.


```r
# what is the maximum allowed ES value (ie which ES value = PSI 100%)?
Maximum.Normalised.ES <- 1/WT.PSI * 100

# split the allowed range of ES values into 1000 evenly distributed numbers
Range.of.ESs <- seq(0,Maximum.Normalised.ES, length.out = 1000)

# calculate corrected ES values
Corrected.ESs <- apply(X = singles[,c("fitness",
                                      "sigma")],
                       MARGIN = 1,
                       FUN = function(x){
                         # enrichment score
                         ES <- exp(x[1])
                         
                         # error associated with this enrichment score
                         Absolute.Sigma <- abs(x[2]*ES)
                         
                         # calculate the probability that the ES is each of the
                         # 1000 evenly distributed numbers, given the measured ES
                         # value and the error associated with that measurement
                         Vector.Densities <- sapply(X = Range.of.ESs,
                                                    FUN = function(y){
                                                      dnorm(ES,
                                                            mean = y,
                                                            sd = Absolute.Sigma)
                                                    })
                         
                         # the corrected enrichment score is given by the weighed
                         # average of these 1000 numbers, where the weights are
                         # equal to the probability of each number
                         Corrected.ES <- weighted.mean(x = Range.of.ESs,
                                                       w = Vector.Densities)
                         
                         # return corrected value
                         Corrected.ES
                       })

# add some columns to our data frame containing the corrected enrichment values
singles$Bounded.Merged.Normalised.ES <- Corrected.ESs
singles$Bounded.Log.Merged.Normalised.ES <- log(Corrected.ESs)

# calculate the PSI values derived from all the corrected enrichment scores 
singles$Bounded.PSI <- (100/Maximum.Normalised.ES)*singles$Bounded.Merged.Normalised.ES

# the distribution of mutant ES's should have a mode centred at the WT ES;
# if there is a small deviation, calculate this deviation 
Density.Object <- density(singles$Bounded.PSI)
WT.Estimated.PSI <- Density.Object$x[which(Density.Object$y == max(Density.Object$y))]
Deviation <- WT.Estimated.PSI - WT.PSI

# and correct the bounded PSI values accounting for this deviation
singles$Bounded.PSI <- singles$Bounded.PSI - Deviation

# finally, correct the PSI values of those mutants whose PSI value was estimated
# to be above 100%
singles$Corrected.PSI <- singles$PSI
singles$Corrected.PSI[which(singles$PSI > 100)] <- singles$Bounded.PSI[which(singles$PSI > 100)]
```


## 4. Effects on splicing efficiency

The effect of each mutation on splicing efficiency is given by equation 1 in figure 1 of the manuscript. I defined a function that converts the phenotypic effect of a mutation on exon inclusion (ΔPSI) into its biophysical effect on splicing efficiency (A):

```r
# function that converts delta PSI values into delta splicing efficiency A
# (from equation 1 in figure 1)
FromDeltaPsi2A <- function(starting.psi, delta.psi){
  s <- starting.psi
  d <- delta.psi
  
  (s^2 - 100*s + d*s - 100*d)/(s*(d + s - 100))
  
}
```
Calculate 'A' for all variants in the dataset:

```r
# calculate 'A'
singles$A <- FromDeltaPsi2A(starting.psi = 49.1,
                            delta.psi = singles$Corrected.PSI - 49.1)
```


## 5. Some extra information & save

Before saving the dataset, add a column to the data frame for the ID of each mutation.

```r
singles$ID <- apply(X = singles,
                    MARGIN = 1,
                    FUN = function(x){
                      paste(toupper(as.character(x[2])),
                            as.numeric(as.character(x[1])),
                            toupper(as.character(x[3])),
                            sep = "-"
                            )
                    })

```
Save:

```r
# save
save(singles,
     file = "Data/Fas_Human_Dataset.RData")
```

