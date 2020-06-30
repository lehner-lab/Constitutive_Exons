# get data
load("Data/WT1_ES.RData")

# psi for the 10 WTs
Ke.WT.PSIs <- c("A" = 7,
                "B" = 20,
                "C" = 65,
                "D" = 0.1,
                "E" = 3,
                "F" = 43,
                "G" = 4,
                "H" = 74,
                "I" = 53,
                "J" = 5)


# function that converts delta PSI values into delta splicing efficiency A
# (from equation 1 in figure 1)
FromDeltaPsiToA <- function(delta.psi, starting.psi){
  y <- delta.psi
  x <- starting.psi
  A <- (y + x - ((y*x + x^2)/(100)))/(x - ((y*x + x^2)/(100)))
  A
}


# loop through the 10 experiments/hexamers/exons in Ke et al 2018
for (i in 1:length(Singles.List)) {
  # the starting PSI is the PSI of the WT
  Singles.List[[i]]$Starting.PSI <- Ke.WT.PSIs[i]
  
  # the final PSI is the PSI after introducing the mutation
  Singles.List[[i]]$Final.PSI <- Singles.List[[i]]$Merged.Normalised.ES * Ke.WT.PSIs[i]
  
  # if we predicted anything to be above 100%, manually set to 100%
  idx <- which(Singles.List[[i]]$Final.PSI > 100)
  if (length(idx)>0){
    Singles.List[[i]]$Final.PSI[idx] <- 100
  }
  
  # if we predicted anything to be below 0%, manually set to 0%
  idx <- which(Singles.List[[i]]$Final.PSI < 0)
  if (length(idx)>0){
    Singles.List[[i]]$Final.PSI[idx] <- 0
  }
  
  # delta PSI
  Singles.List[[i]]$Delta.PSI <- Singles.List[[i]]$Final.PSI - Singles.List[[i]]$Starting.PSI
  
  # delta splicing efficiency
  Singles.List[[i]]$A <- FromDeltaPsiToA(delta.psi = Singles.List[[i]]$Delta.PSI,
                                         starting.psi = Singles.List[[i]]$Starting.PSI)
}

# distribution of A values
Ke.A.Distributions.Singles <- list("B" = Singles.List[[2]]$A,
                                   "C" = Singles.List[[3]]$A,
                                   "F" = Singles.List[[6]]$A,
                                   "H" = Singles.List[[8]]$A,
                                   "I" = Singles.List[[9]]$A)

# save
save(Ke.A.Distributions.Singles,
     file = "Data/Ke_A_Distributions.RData")


# save the whole data frames as well just in case
Singles.DFs <- list("B" = Singles.List[[2]],
                    "C" = Singles.List[[3]],
                    "F" = Singles.List[[6]],
                    "H" = Singles.List[[8]],
                    "I" = Singles.List[[9]])

# save
save(Singles.DFs,
     file = "Data/Ke_Singles_DataFrames.RData")
