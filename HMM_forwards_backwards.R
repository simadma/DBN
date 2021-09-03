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
  gamma <- matrix(0, nrow = N, ncol = K)  # smoothed: gamma[t, i] = P(X_t = i | y_{1:N})
  
  ## Forwards pass
  
  # First step
  phi <- B[, y[1]] * q          # Proportional to P(X_1 | y_1)
  alpha[1, ] <- phi / sum(phi)  # Normalize
  
  # Next steps
  if (N > 1) for (t in 2:N) {
    q <- t(A) %*% alpha[t - 1, ]  # P(X_t | y_{1:(t - 1)})
    phi <- B[, y[t]] * q
    alpha[t, ] <- phi / sum(phi)
  }
  
  ## Backwards pass
  
  # First step
  gamma[N, ] <- alpha[N, ]
  
  # Next steps
  if (N > 1) for (t in (N - 1):1) {
    D <- t(  # D[j, i] = P(X_t = i | X_{t + 1} = j, y_{1:t})
           apply(
             X      = t(A) %*% diag(alpha[t, ]),
             MARGIN = 1, 
             FUN    = function(x) x / sum(x)  # Normalize each row
           )
         )
    gamma[t, ] <- t(D) %*% gamma[t + 1, ]
  }
  
  list(filtered = alpha, smoothed = gamma)
}

##########################################################################################
## Example

K <- 2  # Number of categories for hidden variable
L <- 3  # Number of categories for observed variable

# States
rainy <- 1
sunny <- 2

# Define transition probability matrix
A <- matrix(0, K, K)
A[rainy, rainy] <- 0.7
A[rainy, sunny] <- 0.3
A[sunny, rainy] <- 0.4
A[sunny, sunny] <- 0.6

# Symbols
walk <- 1
shop <- 2
clean <- 3

# Define emission probability matrix
B <- matrix(0, K, L)
B[rainy, walk] <- 0.1
B[rainy, shop] <- 0.4
B[rainy, clean] <- 0.5
B[sunny, walk] <- 0.6
B[sunny, shop] <- 0.3
B[sunny, clean] <- 0.1

# Starting probabilities of the states
start <- c(0.6, 0.4)


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

sum((post - inf_res$smoothed)^2)  # ~zero