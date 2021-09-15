## Simulates a Hidden Markov chain
# q: start probability of states
# A: transition probability matrix
# B: emission probability matrix
# N: number of steps
HMM_sim <- function(q, A, B, N) {
  K <- nrow(B)  # X_t \in {1, 2, ..., K}
  L <- ncol(B)  # Y_t \in {1, 2, ..., L}
  
  X <- numeric(N)
  Y <- numeric(N)
  
  # First step
  X[1] <- sample(K, size = 1, prob = q)
  Y[1] <- sample(L, size = 1, prob = B[X[1], ])
  
  # Next steps
  if (N > 1) for (t in 2:N) {
    X[t] <- sample(K, size = 1, prob = A[X[t - 1], ])
    Y[t] <- sample(L, size = 1, prob = B[X[t], ])
  }
  data.frame(state = X, obs = Y)
}

## Computes the posterior marginals
# q: start probability of states
# A: transition probability matrix
# B: emission probability matrix
# y: observations
forwards_backwards <- function(q, A, B, y) {
  K <- nrow(B)    # X_t \in {1, 2, ..., K}
  L <- ncol(B)    # Y_t \in {1, 2, ..., L}
  N <- length(y)  # t = 1, 2, ..., N
  
  alpha <- matrix(0, nrow = N, ncol = K)  # filtered: alpha[t, i] = P(X_t = i | y_{1:t})
  const <- numeric(N)                     # const[t] = P(y_t | y_{1:t-1})
  xi <- array(0, dim = c(K, K, N - 1))    # xi[i, j, t-1] =
  #   P(X_{t-1} = i, X_t = j | y_{1:N})
  gamma <- matrix(0, nrow = N, ncol = K)  # smoothed: gamma[t, i] = P(X_t = i | y_{1:N})
  
  
  ## Forwards pass
  
  # First step
  t <- 1
  phi <- B[, y[t]] * q          # Proportional to P(X_1 | y_1)
  const[t] <- sum(phi)          # P(y_1)
  alpha[t, ] <- phi / const[t]  # Normalize
  
  # Next steps
  if (N > 1) for (t in 2:N) {
    phi <- B[, y[t]] * (t(A) %*% alpha[t - 1, ])  # Proportional to P(X_t | y_{1:t})
    const[t] <- sum(phi)                          # P(y_t | y_{1:t-1})
    alpha[t, ] <- phi / const[t]                  # Normalize
  }
  
  log_lik <- sum(log(const))  # log P(y_{1:N})
  
  ## Backwards pass
  
  # First step
  gamma[N, ] <- alpha[N, ]
  
  # Next steps
  if (N > 1) for (t in N:2) {
    ratio <- gamma[t, ] / alpha[t, ]              # Update factor
    r <- B[, y[t]] * ratio                        # Multiplied by P(y_t | X_t)
    prop_to_xi <- alpha[t - 1, ] * t(t(A) * r)    # Same as diag(alpha[t-1, ]) A diag(r)
    xi[,, t - 1] <- prop_to_xi / sum(prop_to_xi)  # Two-slice P(X_{t-1}, X_t | y_{1:N})
    gamma[t - 1, ] <- rowSums(xi[,, t - 1])       # Marginalize to P(X_{t-1} | y_{1:N})
  }
  
  list(filtered = alpha, smoothed = gamma, xi = xi, log_lik = log_lik)
}

##########################################################################################
## Example

K <- 2  # Number of categories for hidden variable
L <- 3  # Number of categories for observed variable

# States
rainy <- 1
sunny <- 2

# Define transition probability matrix
A <- matrix(0, K, K, dimnames = list(c("rainy", "sunny"), c("rainy", "sunny")))
A[rainy, rainy] <- 0.7
A[rainy, sunny] <- 0.3
A[sunny, rainy] <- 0.4
A[sunny, sunny] <- 0.6

# Symbols
walk <- 1
shop <- 2
clean <- 3

# Define emission probability matrix
B <- matrix(0, K, L, dimnames = list(c("rainy", "sunny"), c("walk", "shop", "clean")))
B[rainy, walk] <- 0.1
B[rainy, shop] <- 0.4
B[rainy, clean] <- 0.5
B[sunny, walk] <- 0.6
B[sunny, shop] <- 0.3
B[sunny, clean] <- 0.1

# Starting probabilities of the states
start <- c("rainy" = 0.6, "sunny" = 0.4)


set.seed(139)
chain <- HMM_sim(start, A, B, 10000)  # Simulate chain
inf_res <- forwards_backwards(start, A, B, y = chain$obs)


library(HMM)  # Compare result with HMM::posterior()
hmm <- initHMM(
  States        = c(rainy, sunny),
  Symbols       = c(walk, shop, clean),
  startProbs    = start,
  transProbs    = A,
  emissionProbs = B
)

post <- t(posterior(hmm, chain$obs))

head(post)              # posterior marginals
head(inf_res$smoothed)  # identical

norm(post - inf_res$smoothed, type = '2')  # ~zero