source("HMM_forwards_backwards.R")
source("normalize.R")

HMM_learn <- function(y, K, L,
               start   = normalize(runif(K)),
               A       = normalize(matrix(runif(K*K), K, K)),
               B       = normalize(matrix(runif(K*L), K, L)),
               maxiter = 500,
               epsilon = 1e-8
             ) {
  O <- model.matrix(~ factor(y) - 1)  # Convert obs to N by L matrix with L dummy variables
  
  count <- 0
  theta_prev <- 0
  theta <- c(start, A, B)  # Parameter of interest
  log_liks <- c()
  while (epsilon < norm(theta - theta_prev, type = '2') && count < maxiter) {
    theta_prev <- theta
    
    ## E-step
    result <- forwards_backwards(start, A, B, y, learning = TRUE)
    xi <- result$xi
    gamma <- result$smoothed
    log_liks <- c(log_liks, result$log_lik)
    ## M-step
    start <- gamma[1, ]  # P(X_1 = i), i = 1, ..., K
    
    A <- normalize(colSums(xi))     # Sum_{t=2:N}(xi[,, t-1]) row-wise normalized
    
    B <- normalize(t(gamma) %*% O)  # temp[i, j] = Sum_{t=1:N}(gamma_t(i) * I(y[t] = j))
                                    # row-wise normalized
    
    ## Update parameter
    theta <- c(start, A, B)
    
    count <- count + 1
  }
  if (count == maxiter) {
    print("Max iterations reached")
    print(sprintf("Norm: %.3e", norm(theta - theta_prev, type = '2')))
  }
  print(count)
  list(start = start, A = A, B = B, log_lik = log_liks)
}