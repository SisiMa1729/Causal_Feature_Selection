pcsimple_survival<-function (y, dm, alpha, verbose = FALSE, 
          directed = FALSE)
##################################################################################
# function output the Parents and Children set (PC) for a survival variable
# input: y, a survival object, Surv class (produced by the Surv function in survival package)
#        dm, a data frame containing variables that are candidate PC
#       alpha, threshold for determining conditional independency
# output: output$G contain the local graph structure of y (pcselect style)
##################################################################################
# Modified after pcselect function in pcalg by Markus Kalisch.
#
# The original pcselect procedure is modified minimally intentionally for illustrative purpose. 
# 
# This procedure is functional as it is, but not optimal. Some possible enhancements: 
# (1) Number of tests conducted can be reduced. e.g. outcome ~ coxph (X, Y) gives 
# if X is independent with outcome given Y and if Y is independent with outcome given X
# the current version run two cox model instead of one.
# (2) Alternative way of determining conditional independency can be constructed by comparing nested models 
#################################################################################
# Sisi Ma, 11/02/2018 for AMIA Causal Feature Selection Workshop
# sisima@umn.edu
#################################################################################
{
  stopifnot((n <- nrow(dm)) >= 1, (p <- ncol(dm)) >= 1)
  vNms <- colnames(dm)
  dm<-cbind(y,dm)
  # following commented lines prepares for fisher Z test, we do not need this, since 
  # we do not do fisher z for survival outcome
  #zMin <- c(0, rep.int(Inf, p)) 
  #C <- mcor(cbind(y, dm), method = corMethod)
  #cutoff <- qnorm(1 - alpha/2)
  n.edgetests <- numeric(1)
  G <- c(FALSE, rep.int(TRUE, p))
  seq_p <- seq_len(p + 1L)
  done <- FALSE
  ord <- 0
  while (!done && any(G)) {
    n.edgetests[ord + 1] <- 0
    done <- TRUE
    ind <- which(G)
    remEdges <- length(ind)
    if (verbose >= 1) 
      cat("Order=", ord, "; remaining edges:", remEdges, 
          "\n", sep = "")
    for (i in 1:remEdges) {
      if (verbose && (verbose >= 2 || i%%100 == 0)) 
        cat("|i=", i, "|iMax=", remEdges, "\n")
      y <- 1
      x <- ind[i]
      if (G[x]) {
        nbrsBool <- G
        nbrsBool[x] <- FALSE
        nbrs <- seq_p[nbrsBool]
        length_nbrs <- length(nbrs)
        if (length_nbrs >= ord) {
          if (length_nbrs > ord) 
            done <- FALSE
          S <- seq_len(ord)
          repeat {
            n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
            # the following commented lines conduct fisher Z test, we replace them
            # with line 66-68, the conditional independency test for survival
            #z <- zStat(x, y, nbrs[S], C, n)
            #if (abs(z) < zMin[x]) 
            # zMin[x] <- abs(z)  
            cdata<-dm[,c(y,x,nbrs[S])]
            cfit<-coxph(y~.,data=cdata)
            pval<-summary(cfit)$coef[1,5]

            if (verbose >= 2) 
              cat(paste("x:", vNms[x - 1], "y:", (ytmp <- round((p + 
                                                                   1)/2)), "S:"), c(ytmp, vNms)[nbrs[S]], 
                  paste("z:", z, "\n"))
            if (pval > alpha) {
              G[x] <- FALSE
              break
            }
            else {
              nextSet <- getNextSet(length_nbrs, ord, 
                                    S)
              if (nextSet$wasLast) 
                break
              S <- nextSet$nextSet
            }
          }
        }
      }
    }
    ord <- ord + 1
  }
  list(G = setNames(G[-1L], vNms))
}