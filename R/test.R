source("R/normalize.R")
source("R/HMM_forwards_backwards.R")
source("R/HMM_learning.R")
source("R/roll_unroll.R")
source("R/idx_to_vals_and_vals_to_idx.R")

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
N <- 4000
chain <- HMM_sim(start, A, B, N)  # Simulate chain
y <- chain$obs
y[sample(N, 1000)] <- NA
inf_res <- forwards_backwards(start, A, B, y)

library(HMM)  # Compare result with HMM::posterior()
hmm <- initHMM(
  States        = c(rainy, sunny),
  Symbols       = c(walk, shop, clean),
  startProbs    = start,
  transProbs    = A,
  emissionProbs = B
)

post <- t(posterior(hmm, chain$obs))
time <- 850:1000
plot(chain$state[time] - 1,
     type = 'l', main = "Probability of weather being sunny",
     xlab = "time", ylab = "P(X_t = sunny | obs)"
)
lines(inf_res$smoothed[time, 2], col = 'blue', lwd = 2)  # From our algorithm
lines(post[time, 2], col = 'red', lty = 2, lwd = 2)      # From HMM::posterior()
norm(post - inf_res$smoothed, type = '2')  # ~zero

#########################################################################################

set.seed(12)
hmm_new <- hmm
hmm_new$startProbs <- normalize(runif(K))
hmm_new$transProbs <- normalize(matrix(runif(K*K), K, K))
hmm_new$emissionProbs <- normalize(matrix(runif(K*L), K, L))
hmm_new
HMM_BW_res <- baumWelch(hmm_new, y)
HMM_BW_res$hmm
HMM_param <- with(HMM_BW_res$hmm, forwards_backwards(startProbs, transProbs, emissionProbs, y))

res <- HMM_learn(y, start = start, A = A, B = B)
param_adj <- with(res, forwards_backwards(start, A, B, y))

set.seed(13)
rng_res <- HMM_learn(y, K, L)
plot(rng_res$log_lik, type = 'l')
lines(res$log_lik[1:500], col = 'red')
rng_param_adj <- with(rng_res, forwards_backwards(start, A, B, y))


time <- 1:100 + 400
plot(chain$state[time] - 1,
     type = 'l', main = "Probability of weather being sunny",
     xlab = "time", ylab = "P(X_t = sunny | obs)"
)
lines(inf_res$smoothed[time, 2], col = 'blue', lwd = 1)  # From our algorithm
lines(param_adj$smoothed[time, 2], col = 'red', lwd = 1)  # With parameter estimation
lines(rng_param_adj$smoothed[time, 2], col = 'green', lwd = 1)  # From HMM::baumWelch()


##########################################################################################
#                                    Test DBN to HMM                                     #
#                                      Oct 5th 2021                                      #
##########################################################################################
set.seed(3)

H <- 4  # Number of hidden nodes
ns_H <- rep(2, H)  # Binary nodes

O <- 3  # Number of observed nodes
ns_O <- rep(2, O)  # Binary nodes

HH <- prod(ns_H)   # Number of conditional combinations (2^H)

qs <- list()
tms <- list()
for (i in 1:H) {
  ns <- ns_H[i]  # Node size
  qs[[i]] <- normalize(runif(ns))
  tms[[i]] <- normalize(matrix(runif(HH*ns), HH, ns))
}

ems <- list()
for (i in 1:O) {
  ns <- ns_O[i]
  ems[[i]] <- normalize(matrix(runif(HH*ns), HH, ns))
}

joint <- DBN_to_HMM(qs, tms, ems)
TT <- 50  # Number of slices
sim <- HMM_sim(joint$start, joint$A, joint$B, TT)
inf_res <- forwards_backwards(joint$start, joint$A, joint$B, sim$obs)

X_margs <- marginal(inf_res$smoothed, ns_H)

sim$state_vals <- matrix(0, TT, H)
for (t in 1:TT) {
  sim$state_vals[t, ] <- J_inv(sim$state[t], ns_H)
}

sim$state_vals

node <- 1
plot(sim$state_vals[, node] - 1, type = 'l')
lines(X_margs[[node]][, 2], col = 'blue')



##########################################################################################
forwards_backwardsTEST <- function(q, A, B, y) {
  K <- nrow(B)    # X_t \in {1, 2, ..., K}
  L <- ncol(B)    # Y_t \in {1, 2, ..., L}
  N <- length(y)  # t = 1, 2, ..., N
  
  alpha <- matrix(0, nrow = N, ncol = K)  # filtered: alpha[t, i] = P(X_t = i | y_{1:t})
  const <- numeric(N)                     # const[t] = P(y_t | y_{1:t-1})
  xi <- array(0, dim = c(K, K, N - 1))    # xi[i, j, t-1] =
  #   P(X_{t-1} = i, X_t = j | y_{1:N})
  gamma <- matrix(0, nrow = N, ncol = K)  # smoothed: gamma[t, i] = P(X_t = i | y_{1:N})
  
  o <- matrix(0, N, K)  # o[t, i] = P(Y^{(t)} = y^{(t)} | X^{(t)} = i)
  nas <- is.na(y)
  if (any(nas)) {  # If missing observations
    o[!nas, ] <- t(B[, y[!nas]])
    o[nas, ] <- 1/L
    # o[nas, ] <- rep(rowSums(B*B), each = sum(nas))  # o[t, i] E[P(Y^{(t)} | X^{(t)} = i)]
  } else {
    o <- t(B[, y])
  }
  
  
  ## Forwards pass
  t <- 1
  phi <- o[t, ] * q          # Proportional to P(X_1 | y_1)
  const[t] <- sum(phi)          # P(y_1)
  alpha[t, ] <- phi / const[t]  # Normalize
  if (N > 1) for (t in 2:N) {
    phi <- o[t, ] * (t(A) %*% alpha[t - 1, ])  # Proportional to P(X_t | y_{1:t})
    const[t] <- sum(phi)                          # P(y_t | y_{1:t-1})
    alpha[t, ] <- phi / const[t]                  # Normalize
  }
  
  log_lik <- sum(log(const))                      # Log likelihood = log P(y_{1:N})
  
  ## Backwards pass
  gamma[N, ] <- alpha[N, ]
  if (N > 1) for (t in N:2) {
    ratio <- gamma[t, ] / alpha[t, ]              # Update factor
    r <- o[t, ] * ratio                        # Multiplied by P(y_t | X_t)
    prop_to_xi <- alpha[t - 1, ] * t(t(A) * r)    # Same as diag(alpha[t-1, ]) A diag(r)
    xi[,, t - 1] <- prop_to_xi / sum(prop_to_xi)  # Two-slice P(X_{t-1}, X_t | y_{1:N})
    gamma[t - 1, ] <- rowSums(xi[,, t - 1])       # Marginalize to P(X_{t-1} | y_{1:N})
  }
  
  list(filtered = alpha, smoothed = gamma, xi = xi, log_lik = log_lik)
}

inf_res <- forwards_backwards(start, A, B, y)
inf_resTEST <- forwards_backwardsTEST(start, A, B, y)




time <- 1:50 + 50*79
plot(chain$state[time] - 1,
     type = 'l', main = "Probability of weather being sunny",
     xlab = "time", ylab = "P(X_t = sunny | obs)"
)
lines(inf_res$smoothed[time, 2], col = 'blue', lwd = 1)  # From our algorithm
lines(inf_resTEST$smoothed[time, 2], col = 'red', lwd = 1)  # From our algorithm







A_1 <- list(
  var = c("A", "B"),
  ns = c(2, 2),
  mat = matrix(c(0.2, 0.8, 0.4, 0.6), nrow = 2, byrow = TRUE)
)
A_1
