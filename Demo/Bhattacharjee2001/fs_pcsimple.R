fs_pcsimple <- function(x,y, train_idx){
# function selects the parent children (pc) set of a variable (y); a wrapper on top of psSelect function from pcalg package
# input: x, each column is a variable, each row a observation. candidate variable set for pc
#        y, outcome of interest
#        train_idx, a vector sepcifying the set of observations used for training.
#################################################################################
# Sisi Ma, 11/02/2018 for AMIA Causal Feature Selection Workshop
# sisima@umn.edu
#################################################################################
  library(pcalg)
  tr_x<-x[train_idx,]
  tr_y<-y[train_idx]
  slct<-pcSelect(as.numeric(tr_y), as.matrix(tr_x), 0.001, corMethod = "standard",
    verbose = FALSE, directed = FALSE)
  features<-which(slct$G)
  return(features)
}
