source("R/idx_to_vals_and_vals_to_idx.R")
source("R/mk_cpt.R")

factor_product <- function(A, B) {
  # If any of A or B is an empty factor
  if (length(A$var) == 0) return(B)
  if (length(B$var) == 0) return(A)
  
  C <- list()  # Initiate joint factor C
  C$var <- union(A$var, B$var)
  
  mapA <- match(A$var, C$var)  # Indices in C$var where variables exist in A$var
  mapB <- match(B$var, C$var)  # Indices in C$var where variables exist in B$var
  C$ns <- numeric(length(C$var))
  C$ns[mapA] <- A$ns
  C$ns[mapB] <- B$ns
  
  K <- prod(C$ns)
  vals <- idx_to_vals(1:K, C$ns)
  idxA <- vals_to_idx(vals[, mapA, drop = FALSE], A$ns)
  idxB <- vals_to_idx(vals[, mapB, drop = FALSE], B$ns)
  C$val <- numeric(K)
  C$val <- A$val[idxA] * B$val[idxB]
  C
}

table_to_factor <- function(tables) {
  factors <- list()
  tab_names <- names(tables)
  for (i in 1:length(tables)) {
    # name <- tab_names[i]  # Name of current table
    # t <- tables[[name]]   # Current table
    # l <- list()           # Create list to store table as a factor
    # l$var <- with(t$var, c(prob_of, given))
    # l$ns <- with(t$ns, c(prob_of, given))
    # l$val <- c(t$prob)      # Vectorize emission/transition probability matrix
    # factors[[name]] <- l    # Store factor in factor list
    t <- tables[[i]]  # Current table
    l <- list()       # Create list to store table as a factor
    l$var <- with(t$var, c(prob_of, given))
    l$ns <- with(t$ns, c(prob_of, given))
    l$val <- c(t$prob)  # Vectorize emission/transition probability matrix
    factors[[i]] <- l   # Store factor in factor list
  }
  factors
}

factor_to_table <- function(fac, given = NULL) {
  vals <- with(fac, idx_to_vals(1:length(val), ns))
  are_given <- with(fac, var %in% given)  # Logical vector where given variables are TRUE
  ord_idx <- c(which(!are_given), which(are_given))  # Reorder indices s.t. current
                                                     # variables come before previous
  i_vals <- with(fac, order(vals_to_idx(vals[, ord_idx], ns[ord_idx])))  # Same as above
                                                                         # but with probs
  table <- list()
  table$var <- with(fac, list(prob_of = var[!are_given],
                              given   = var[are_given]))
  table$ns <- with(fac, list(prob_of = ns[!are_given],
                             given   = ns[are_given]))
  table$prob <- with(fac, matrix(val[i_vals], nrow = prod(ns[are_given]))[,, drop = TRUE])
  table
}

joint_CPT <- function(tables) {
  if (length(tables) == 1) {
    return(tables[[1]])
  }
  factors <- table_to_factor(tables)  # Turn probability matrices into factors
  joint_fac <- list(var = NULL, ns = NULL, val = NULL)
  for (fac in factors) {
    joint_fac <- factor_product(joint_fac, fac)
  }
  # Here we figure out what is given in the end after we join the tables
  prob_of <- c()
  given <- c()
  for (tab in tables) {
    prob_of <- union(prob_of, tab$var$prob_of)
    given <- union(given, tab$var$given)
  }
  given <- setdiff(given, prob_of)  # This is what is given in the joint table
  joint_tab <- factor_to_table(joint_fac, given)
  joint_tab
}

marginal <- function(gamma, ns) {
  if (is.vector(gamma)) gamma <- t(gamma)
  TT <- nrow(gamma)  # Number of slices
  K <- ncol(gamma)   # Number of combinations of X_1^(t), X_2^(t), ..., X_N^(t) values
  N <- length(ns)    # Number of nodes
  X_margs <- list()  # List to store marginal posteriors
  
  # Store all combinations of X_1^(t), X_2^(t), ..., X_N^(t) values in matrix `vals` below
  #   1   1 ...   1
  #   1   1 ...   2
  #   1   1 ...   3
  #   :   :       :
  # k_1 k_2 ... k_N
  vals <- idx_to_vals(1:K, ns)
  for (node in 1:N) {
    marg <- matrix(0, TT, ns[node])
    for (j in 1:ns[node]) {
      locs <- vals[, node] == j  # Locations where current node has value j
      marg[, j] <- rowSums(gamma[, locs, drop = FALSE])  # Marginalize
    }
    X_margs[[node]] <- marg
  }
  X_margs
}

cpt_from_joint_cpt <- function(joint_cpt, prob_of, given) {
  # joint_cpt
  #   var
  #     prob_of
  #     given
  #   ns
  #     prob_of
  #     given
  #   prob
  cond_idx <- match(given, joint_cpt$var$given)
  cond_ns <- joint_cpt$ns$given[cond_idx]
  cond_vals <- idx_to_vals(1:prod(cond_ns), cond_ns)
  cond_vals_full <- matrix(1, nrow(cond_vals), length(joint_cpt$var$given))
  cond_vals_full[, cond_idx] <- cond_vals  # Now we have the value combinations needed
  rows <- vals_to_idx(cond_vals_full, joint_cpt$ns$given)  # These are the rows we will use
  
  
  var_idx <- match(prob_of, joint_cpt$var$prob_of)
  var_ns <- joint_cpt$ns$prob_of[var_idx]
  var_vals <- with(joint_cpt$ns, idx_to_vals(1:prod(prob_of), prob_of))
  table <- matrix(0, length(rows), var_ns)
  for (i in 1:var_ns) {
    locs <- var_vals[, var_idx] == i
    table[, i] <- rowSums(joint_cpt$prob[rows, locs, drop = FALSE])
  }
  mk_cpt(table=table, prob_of=prob_of, given = given, given_ns = cond_ns)
}