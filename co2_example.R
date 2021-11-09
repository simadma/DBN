source("simulation.R")
source("roll_unroll.R")
source("HMM_forwards_backwards.R")
source("idx_to_vals_and_vals_to_idx.R")
source("mk_cpt.R")

# prob_leakage <- 0.025
# N <- 5
# tables <- create_CO2_DBN(prob_leakage, N)
# 
# q <- with(tables, joint_CPT(q_CPTs))
# A <- with(tables, joint_CPT(A_CPTs))
# B <- with(tables, joint_CPT(B_CPTs))
# 
# set.seed(28)
# TT <- 100  # Number of slices
# chain <- HMM_sim(q$prob, A$prob, B$prob, TT)
# # chain$state
# # chain$obs
# state_vals <- with(chain, idx_to_vals(state, rep(2, N))) - 1
# obs_vals <- with(chain, idx_to_vals(obs, rep(2, N))) - 1
# state_vals
# 
# y <- chain$obs
# y[sample(TT, size=80)] <- NA
# 
# inf_res <- forwards_backwards(q$prob, A$prob, B$prob, y)
# 
# margs <- marginal(inf_res$smoothed, A$ns$prob_of)
# 
# plot(margs[[1]][, 2], type = 'l', col=1, ylim=c(-0.1, 1.1))
# lines(state_vals[, 1], col=1, lty=2, lwd=2)
# lines(margs[[2]][, 2], col=2)
# lines(state_vals[, 2], col=2, lty=2, lwd=2)
# lines(margs[[3]][, 2], col=3)
# lines(state_vals[, 3], col=3, lty=2, lwd=2)
# lines(margs[[4]][, 2], col=4)
# lines(state_vals[, 4], col=4, lty=2, lwd=2)


# idx_to_vals(which(A$prob[1, ] > 0), ns = rep(2, N)) - 1
# A$prob




create_CO2_DBN <- function(N = 5, lambda, nu, kappa1, kappa2, thetaX, thetaP) {
  # Start probability tables
  q_CPTs <- list()
  q_CPTs[["C"]] <- mk_cpt(c(1/3, 1/3, 1/3), "C")
  q_CPTs[["X1"]] <- mk_cpt(c(0, 1), "X1")
  q_CPTs[["P1"]] <- mk_cpt(
    table    = matrix(c(1, 0, 0,
                        0, 0, 1), ncol = 3, byrow = T),
    prob_of  = "P1",
    given    = "X1",
    given_ns = 2
  )
  for (n in 2:N) {
    q_CPTs[[paste0("X", n)]] <- mk_cpt(      #   C    X_{n-1} P_{n-1}
      table    = matrix(c(         1,      0,  # 1       1       1
                                   1,      0,  #                 2
                                   1,      0,  #                 3
                          1 - lambda, lambda,  #         2       1
                                   0,      1,  #                 2
                                   0,      1,  #                 3
                                   1,      0,  # 2       1       1
                                   1,      0,  #                 2
                                   1,      0,  #                 3
                                   1,      0,  #         2       1
                          1 - lambda, lambda,  #                 2
                                   0,      1,  #                 3
                                   1,      0,  # 3       1       1
                                   1,      0,  #                 2
                                   1,      0,  #                 3
                                   1,      0,  #         2       1
                                   1,      0,  #                 2
                          1 - lambda, lambda), #                 3
                        ncol = 2, byrow = T),
      prob_of  = paste0("X", n),
      given    = c("C", paste0(c("X", "P"), n - 1)),
      given_ns = c(3, 2, 3)
    )
    if (n < N) {
      q_CPTs[[paste0("P", n)]] <- mk_cpt(                       #   P_{n-1}   X_n
        table    = matrix(c(     1,      0,                   0,  #   1        1
                                 1,      0,                   0,  #            2
                                 1,      0,                   0,  #   2        1
                            1 - nu,     nu,                   0,  #            2
                                 1,      0,                   0,  #   3        1
                            kappa1, kappa2, 1 - kappa1 - kappa2), #            2
                          ncol = 3, byrow = T),
        prob_of  = paste0("P", n),
        given    = paste0(c("P", "X"), (n - 1):n),
        given_ns = c(3, 2)
      )
    }
  }
  
  # Transition probability tables
  p <- 0.1
  A_CPTs <- list()
  A_CPTs[["C"]] <- mk_cpt(
    table    = diag(3),
    prob_of  = "C",
    given    = "C_prev",
    given_ns = 3
  )
  A_CPTs[["X1"]] <- mk_cpt(
    table    = diag(2),
    prob_of  = "X1",
    given    = "X1_prev",
    given_ns = 2
  )
  A_CPTs[["P1"]] <- mk_cpt(     #   P^(t-1)  X^(t)
    table    = matrix(c(1, 0, 0,  #   1       1
                        1, 0, 0,  #           2
                        0, 1, 0,  #   2       1
                        0, 1, 0,  #           2
                        0, 0, 1,  #   3       1
                        0, 0, 1), #           2
                      ncol = 3, byrow = T),
    prob_of  = "P1",
    given    = c("P1_prev", "X1"),
    given_ns = c(3, 2)
  )
  for (n in 2:N) {
    A_CPTs[[paste0("X", n)]] <- mk_cpt(       #  C^(t)  X_n^(t-1)  X_{n-1}^(t) P_{n-1}^(t)
      table    = matrix(c(         1,      0,  #   1        1          1           1
                                   1,      0,  #                                   2
                                   1,      0,  #                                   3
                          1 - lambda, lambda,  #                       2           1
                                   0,      1,  #                                   2
                                   0,      1,  #                                   3
                                   0,      1,  #            2          1           1
                                   0,      1,  #                                   2
                                   0,      1,  #                                   3
                                   0,      1,  #                       2           1
                                   0,      1,  #                                   2
                                   0,      1,  #                                   3
                                   1,      0,  #   2        1          1           1
                                   1,      0,  #                                   2
                                   1,      0,  #                                   3
                                   1,      0,  #                       2           1
                          1 - lambda, lambda,  #                                   2
                                    0,     1,  #                                   3
                                    0,     1,  #            2          1           1
                                    0,     1,  #                                   2
                                    0,     1,  #                                   3
                                    0,     1,  #                       2           1
                                    0,     1,  #                                   2
                                    0,     1,  #                                   3
                                    1,     0,  #   3        1          1           1
                                    1,     0,  #                                   2
                                    1,     0,  #                                   3
                                    1,     0,  #                       2           1
                                    1,     0,  #                                   2
                          1 - lambda, lambda,  #                                   3
                                    0,     1,  #            2          1           1
                                    0,     1,  #                                   2
                                    0,     1,  #                                   3
                                    0,     1,  #                       2           1
                                    0,     1,  #                                   2
                                    0,     1), #                                   3
                        ncol = 2, byrow = T),
      prob_of  = paste0("X", n),
      given    = c("C", paste0(c("X", "X", "P"), c(n, n - 1, n - 1), c("_prev", "", ""))),
      given_ns = c(3, 2, 2, 3)
    )
    if (n < N) {
      A_CPTs[[paste0("P", n)]] <- mk_cpt(                  # P_n^(t-1) P_{n-1}^(t) X_n^(t)
        table    = matrix(c(     1,      0,                   0,  # 1        1        1
                                 1,      0,                   0,  #                   2
                                 1,      0,                   0,  #          2        1
                            1 - nu,     nu,                   0,  #                   2
                                 1,      0,                   0,  #          3        1
                            kappa1, kappa2, 1 - kappa1 - kappa2,  #                   2
                                 0,      1,                   0,  # 2        1        1
                                 0,      1,                   0,  #                   2
                                 0,      1,                   0,  #          2        1
                                 0,      1,                   0,  #                   2
                                 0,      1,                   0,  #          3        1
                                 0, 1 - nu,                  nu,  #                   2
                                 0,      0,                   1,  # 3        1        1
                                 0,      0,                   1,  #                   2
                                 0,      0,                   1,  #          2        1
                                 0,      0,                   1,  #                   2
                                 0,      0,                   1,  #          3        1
                                 0,      0,                   1), #                   2
                          ncol = 3, byrow = T),
        prob_of  = paste0("P", n),
        given    = paste0(c("P", "P", "X"), c(n, n - 1, n), c("_prev", "", "")),
        given_ns = c(3, 3, 2)
      )
    }
  }
  
  # Emission probability tables
  B_CPTs <- list()
  for (n in 1:N) {
    B_CPTs[[paste0("Y", n)]] <- mk_cpt(
      table    = matrix(c(    thetaX, 1 - thetaX,
                          1 - thetaX,     thetaX),
                        ncol = 2, byrow = T),
      prob_of  = paste0("Y", n),
      given    = paste0("X", n),
      given_ns = 2
    )
    if (n < N) {
      B_CPTs[[paste0("Z", n)]] <- mk_cpt(
        table    = matrix(c(           thetaP, 0.75*(1 - thetaP), 0.25*(1 - thetaP),
                            0.50*(1 - thetaP),            thetaP, 0.50*(1 - thetaP),
                            0.25*(1 - thetaP), 0.75*(1 - thetaP),           thetaP),
                          ncol = 3, byrow = T),
        prob_of  = paste0("Z", n),
        given    = paste0("P", n),
        given_ns = 3
      )
    }
  }
  list(q_CPTs = q_CPTs, A_CPTs = A_CPTs, B_CPTs = B_CPTs)
}

condition_on_C <- function(tables, C = 1) {
  ns <- tables$C$ns$prob_of
  tables$C <- NULL
  for (idx in which(startsWith(names(tables), "X"))[-1]) {
    tables[[idx]]$var$given <- tables[[idx]]$var$given[-1]
    tables[[idx]]$ns$given <- tables[[idx]]$ns$given[-1]
    nrows <- nrow(tables[[idx]]$prob) / ns
    tables[[idx]]$prob <- tables[[idx]]$prob[(C - 1)*nrows + 1:nrows, ]
  }
  joint_CPT(tables)
}


N <- 4
lambda <- 0.10
nu <- 0.10
kappa1 <- 0.20
kappa2 <- 0.30
thetaX <- 0.70
thetaP <- 0.65
tables <- create_CO2_DBN(N, lambda, nu, kappa1, kappa2, thetaX, thetaP)

# Compute joint distributions given C
q_given_C <- list()
A_given_C <- list()
for (C in 1:tables$q_CPTs$C$ns$prob_of) {
  q_given_C[[C]] <- condition_on_C(tables$q_CPTs, C)
  A_given_C[[C]] <- condition_on_C(tables$A_CPTs, C)
}
B <- with(tables, joint_CPT(B_CPTs))


set.seed(21)
TT <- 20  # Number of slices
C <- sample(3, size = 1)  # Draw capillary threshold from {1, 2, 3} = {Low, Medium, High}
chain <- HMM_sim(q_given_C[[C]]$prob, A_given_C[[C]]$prob, B$prob, TT)

state_vals <- with(chain, idx_to_vals(state, B$ns$given))
colnames(state_vals) <- B$var$given
obs_vals <- with(chain, idx_to_vals(obs, B$ns$prob_of))
colnames(obs_vals) <- B$var$prob_of

state_vals[, c(1, 3, 5, 7)] - 1
obs_vals[, c(1, 3, 5, 7)] - 1

inf_res <- list()
for (C in 1:tables$q_CPTs$C$ns$prob_of) {
  inf_res[[C]] <- forwards_backwards(
    q = q_given_C[[C]]$prob,
    A = A_given_C[[C]]$prob,
    B = B$prob,
    y = chain$obs
  )
}
logLiks <- c(inf_res[[1]]$log_lik, inf_res[[2]]$log_lik, inf_res[[3]]$log_lik)
post_c <- normalize(exp(logLiks - max(logLiks))*tables$q_CPTs$C$prob)
post_c
gamma_mid <- post_c[1]*inf_res[[1]]$smoothed + 
  post_c[2]*inf_res[[2]]$smoothed +
  post_c[3]*inf_res[[3]]$smoothed
margs <- marginal(gamma_mid, B$ns$given)

plot(margs[[1]][, 2], type = 'l', col=1, ylim=c(-0.1, 1.1))
lines(state_vals[, 1] - 1, col=1, lty=2, lwd=2)
lines(margs[[3]][, 2], col=2)
lines(state_vals[, 3] - 1, col=2, lty=2, lwd=2)
lines(margs[[5]][, 2], col=3)
lines(state_vals[, 5] - 1, col=3, lty=2, lwd=2)
lines(margs[[7]][, 2], col=4)
lines(state_vals[, 7] - 1, col=4, lty=2, lwd=2)


# If we knew C:
margs <- marginal(inf_res[[3]]$smoothed, B$ns$given)

plot(margs[[1]][, 2], type = 'l', col=1, ylim=c(-0.1, 1.1))
lines(state_vals[, 1] - 1, col=1, lty=2, lwd=2)
lines(margs[[3]][, 2], col=2)
lines(state_vals[, 3] - 1, col=2, lty=2, lwd=2)
lines(margs[[5]][, 2], col=3)
lines(state_vals[, 5] - 1, col=3, lty=2, lwd=2)
lines(margs[[7]][, 2], col=4)
lines(state_vals[, 7] - 1, col=4, lty=2, lwd=2)