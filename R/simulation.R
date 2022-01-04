## Simulates a Hidden Markov chain
# q:  start probability of states
# A:  transition probability matrix
# B:  emission probability matrix
# TT: number of steps
HMM_sim <- function(q, A, B, TT) {
  K <- nrow(B)  # X^{(t)} \in {1, 2, ..., K}
  L <- ncol(B)  # Y^{(t)} \in {1, 2, ..., L}
  
  X <- numeric(TT)
  Y <- numeric(TT)
  
  # First step
  X[1] <- sample(K, size = 1, prob = q)
  Y[1] <- sample(L, size = 1, prob = B[X[1], ])
  
  # Next steps
  if (TT > 1) for (t in 2:TT) {
    X[t] <- sample(K, size = 1, prob = A[X[t - 1], ])
    Y[t] <- sample(L, size = 1, prob = B[X[t], ])
  }
  data.frame(state = X, obs = Y)
}

## Simulates a (stationary) Markov chain
# q:  start probability of states
# P:  transition probability matrix
# N: number of steps
Markov_sim <- function(q, P, N) {
  K <- length(q)  # X_n \in {1, 2, ..., K}
  X <- numeric(N)
  
  # First step
  X[1] <- sample(K, size = 1, prob = q)
  
  # Next steps
  if (N > 1) for (t in 2:N) {
    X[t] <- sample(K, size = 1, prob = P[X[t - 1], ])
  }
  X
}