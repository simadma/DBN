## Computes the node values from unique integer index
# n:  integer index
# ns: node sizes
## EXAMPLE
# n = 8, ns = (3, 4)
# returns (2, 4), because it is the eighth element in
# 1:     (1, 1),  2:     (1, 2),  3:     (1, 3),  4:     (1, ns[2]),
# 5:     (2, 1),  6:     (2, 2),  7:     (2, 3),  8:     (2, ns[2]),
# 9: (ns[1], 1), 10: (ns[1], 2), 11: (ns[1], 3), 12: (ns[1], ns[2])
J <- function(n, ns) {
  s <- length(ns)
  if (s == 1) return(n)
  NN <- prod(ns[2:s])
  x <- ((n - 1) %/% NN) + 1  # Index of node
  m <- ((n - 1) %% NN) + 1   # Remaining integer index
  c(x, J(m, ns[2:s]))
}

## Unnests values and returns a unique integer index
# vals: values of X_1^(t), X_2^(t), ... X_H^(t)
# ns:   node sizes
## EXAMPLE
# n = (2, 4), ns = (3, 4)
# returns 8, because it is the eighth element in
# 1: (    1, 1),  2: (    1, 2),  3: (    1, 3),  4: (    1, ns[2]),
# 5: (    2, 1),  6: (    2, 2),  7: (    2, 3),  8: (    2, ns[2]),
# 9: (ns[1], 1), 10: (ns[1], 2), 11: (ns[1], 3), 12: (ns[1], ns[2])
J_inv <- function(vals, ns) {
  s <- length(ns)  # same as length(vals)
  if (s == 1) return(vals)
  (vals[1] - 1)*prod(ns[2:s]) + J_inv(vals[2:s], ns[2:s])
}


## Returns a matrix where row i is J(idx[i], ns)
# Look at J()-function for details
idx_to_vals <- function(idx, ns) {
  N <- length(idx)
  vals <- matrix(0, nrow = N, ncol = length(ns))
  for (i in 1:N) {
    vals[i, ] <- J(idx[i], ns)
  }
  vals
}

## Returns a vector where element i is J_inv(vals[i, ], ns)
# Look at J_inv()-function for details
vals_to_idx <- function(vals, ns) {
  N <- nrow(vals)
  idx <- numeric(N)
  for (i in 1:N) {
    idx[i] <- J_inv(vals[i, ], ns)
  }
  idx
}

##########################################################################################
# Sanity check
# set.seed(1)
# ns <- sample(2:5, size = 4, replace = TRUE)
# vals_to_idx(idx_to_vals(1:prod(ns), ns), ns) - 1:prod(ns)
