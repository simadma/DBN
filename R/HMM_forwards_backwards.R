source("R/normalize.R")

## Computes the posterior marginals
# q: start probability of states
# A: transition probability matrix
# B: emission probability matrix
# y: observations
forwards_backwards <- function(q, A, B, y, learning = FALSE) {
  K <- nrow(B)    # X_t \in {1, 2, ..., K}
  L <- ncol(B)    # Y_t \in {1, 2, ..., L}
  N <- length(y)  # t = 1, 2, ..., N
  
  const <- numeric(N)                     # const[t] = P(y_t | y_{1:t-1})
  o <- matrix(0, N, K)  # o[t, i] = P(Y^{(t)} = y^{(t)} | X^{(t)} = i)
  nas <- is.na(y)
  if (any(nas)) {  # If missing observations
    o[!nas, ] <- t(B[, y[!nas]])
    # For missing obs, we use the expectation, i.e. o[t, i] = E[P(Y^{(t)} | X^{(t)} = i)]
    o[nas, ] <- rep(rowSums(B*B), each = sum(nas))
  } else {
    o <- t(B[, y])
  }
  
  
  ## Forwards pass
  alpha <- matrix(0, nrow = N, ncol = K)  # filtered: alpha[t, i] = P(X_t = i | y_{1:t})
  t <- 1
  phi <- o[t, ] * q          # Proportional to P(X_1 | y_1)
  const[t] <- sum(phi)          # P(y_1)
  alpha[t, ] <- phi / const[t]  # Normalize
  if (N > 1) for (t in 2:N) {
    phi <- o[t, ] * (t(A) %*% alpha[t - 1, ])     # Proportional to P(X_t | y_{1:t})
    const[t] <- sum(phi)                          # P(y_t | y_{1:t-1})
    alpha[t, ] <- phi / const[t]                  # Normalize
  }
  
  log_lik <- sum(log(const))                      # Log likelihood = log P(y_{1:N})
  
  ## Backwards pass
  beta <- matrix(0, nrow = N, ncol = K)
  beta[N, ] <- 1
  if (N > 1) for (t in N:2) {
    # To ensure numerical stability, we normalize at each step
    # This is also the same as normalize(A %*% diag(o[t, ]) %*% beta[t, ])
    beta[t - 1, ] <- normalize(A %*% (o[t, ] * beta[t, ]))  
  }
  
  ## Posterior probability
  gamma <- normalize(alpha * beta)  # smoothed: gamma[t, i] = P(X_t = i | y_{1:N})
  
  # If we want to learn the parameters, we compute the two-slice distribution
  # P(X_{t-1} = i, X_t = j | y_{1:N}) for i,j = 1, ..., K
  if (learning) {
    xi <- array(0, dim = c(N - 1, K, K))    # xi[t-1, i, j] =
                                            #   P(X_{t-1} = i, X_t = j | y_{1:N})
    for (t in 2:N) {
      prop_to_xi <- A * (alpha[t - 1, ] %*% t(o[t, ] * beta[t, ]))
      xi[t - 1,, ] <- prop_to_xi / sum(prop_to_xi)
    }
    return(list(filtered = alpha, smoothed = gamma, xi = xi, log_lik = log_lik))
  }
  
  # gamma[N, ] <- alpha[N, ]
  # if (N > 1) for (t in N:2) {
  #   if (any(alpha[t,] == 0)) print("DIVISION BY ZERO!!!!!")
  #   ratio <- gamma[t, ] / alpha[t, ]              # Update factor
  #   r <- o[t, ] * ratio                           # Multiplied by P(y_t | X_t)
  #   prop_to_xi <- alpha[t - 1, ] * t(t(A) * r)    # Same as diag(alpha[t-1, ]) A diag(r)
  #   xi[,, t - 1] <- prop_to_xi / sum(prop_to_xi)  # Two-slice P(X_{t-1}, X_t | y_{1:N})
  #   gamma[t - 1, ] <- rowSums(xi[,, t - 1])       # Marginalize to P(X_{t-1} | y_{1:N})
  # }
  
  list(filtered = alpha, smoothed = gamma, log_lik = log_lik)
}