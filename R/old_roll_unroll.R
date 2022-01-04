source("R/idx_to_vals_and_vals_to_idx.R")
source("R/normalize.R")

## Computes the joint table of list of tables
# ts: list of tables
joint_table <- function(ts) {
  N <- length(ts)  # Number of nodes
  
  # Node sizes (number of categories per node)
  ns <- numeric(N)
  for (i in 1:N) {
    if (is.vector(ts[[i]])) {
      ts[[i]] <- t(ts[[i]])  # Transform to row vector
    }
    ns[i] <- ncol(ts[[i]])
  }
  
  MM <- nrow(ts[[1]])  # Number of conditional combinations
  NN <- prod(ns)        # Node size of joint distribution
  
  P <- matrix(0, MM, NN)  # Initiate joint probability matrix
  for (i in 1:NN) {
    vals <- J(i, ns)  # Values of X_1^(t), X_2^(t), ..., X_N^(t)
    
    # Computes joint probability table
    temp <- 1
    for (j in 1:N) {
      t_j <- ts[[j]]  # Prior probability table for hidden node X_j
      temp <- temp * t_j[, vals[j]]
    }
    P[, i] <- temp    # Assigns P(vals | condition)
  }
  P[,, drop = TRUE]  # Transforms to vector if nrow(P) is 1.
}


## Rolls the marginal prior, transition and emission tables into a joint prior, joint
## transition and joint emission table
# qs:  list of marginal prior tables
# tms: list of marginal transition tables
# ems: list of marginal emission tables
DBN_to_HMM <- function(qs, tms, ems) {
  list(start = joint_table(qs), A = joint_table(tms), B = joint_table(ems))
}

marginal <- function(gamma, ns) {
  TT <- nrow(gamma)  # Number of slices
  NN <- ncol(gamma)  # Number of combinations of X_1^(t), X_2^(t), ..., X_H^(t) values
  N <- length(ns)    # Number of nodes
  X_margs <- list()  # List to store marginal posteriors
  
  # Store all combinations of X_1^(t), X_2^(t), ..., X_H^(t) values in matrix `vals` below
  #   1   1 ...   1
  #   1   1 ...   2
  #   1   1 ...   3
  #   :   :       :
  # h_1 h_2 ... h_H
  vals <- idx_to_vals(1:NN, ns)
  for (node in 1:N) {
    marg <- matrix(0, TT, ns[node])
    for (j in 1:ns[node]) {
      locs <- vals[, node] == j            # Locations where current node has value j
      marg[, j] <- rowSums(gamma[, locs])  # Marginalize
    }
    X_margs[[node]] <- marg
  }
  X_margs
}