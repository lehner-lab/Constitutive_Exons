
library(data.table)

DNA.Letters <- c("A", "T", "G", "C")

All.Hexamers <- apply(X = gtools::permutations(n = 4,
                                               r = 6,
                                               v = DNA.Letters,
                                               repeats.allowed = T),
                      MARGIN = 1,
                      FUN = function(x){
                        paste0(x, collapse='')
                      })





# all hexamers one Hamming distance away from each other hexamer





Hexamers1HammingDistanceAway <- function(this.hexamer, all.hexamers) {
  this.hexamer.vectorised <- strsplit(x = this.hexamer,
                                      split = "")[[1]]
  hamming.distances <- sapply(X = strsplit(all.hexamers,
                                           split = ""),
                              FUN = function(x){
                                sum(x != this.hexamer.vectorised)
                              })
  hamming.distance.1 <- which(hamming.distances == 1)
  all.hexamers[hamming.distance.1]
}





Hexamers.1.Hamming.Distance.Away <- sapply(X = All.Hexamers,
                                           FUN = function(x){
                                             Hexamers1HammingDistanceAway(this.hexamer = x,
                                                                          all.hexamers = All.Hexamers)
                                           })




# build a network table
Hexamer.Network <- data.table(Main = c(),
                              Neighbour = c())

for (each.column.name in colnames(Hexamers.1.Hamming.Distance.Away)) {
  print(each.column.name)
  temporary.data.table <- data.table(Main = each.column.name,
                                     Neighbour = Hexamers.1.Hamming.Distance.Away[,each.column.name])
  Hexamer.Network <- rbind(Hexamer.Network,
                           temporary.data.table)
}





# count how many (overlapping) ESEs we find, using
ESE <- readxl::read_excel(path = "238_Supplemental_Table_1.xls",
                          sheet = 2,
                          range = "A3:B1184",
                          col_names = c("Hexamer", "Score"))

ESS <- readxl::read_excel(path = "238_Supplemental_Table_1.xls",
                          sheet = 2,
                          range = "A1187:B2276",
                          col_names = c("Hexamer", "Score"))


ESE.Vector <- ESE$Score
names(ESE.Vector) <- ESE$Hexamer


ESS.Vector <- ESS$Score
names(ESS.Vector) <- ESS$Hexamer



save(Hexamer.Network,
     ESE,
     ESE.Vector,
     ESS,
     ESS.Vector,
     file = "239bis_Hexamer_Table.RData")
