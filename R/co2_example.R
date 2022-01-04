source("R/simulation.R")
source("R/roll_unroll.R")
source("R/HMM_forwards_backwards.R")
source("R/idx_to_vals_and_vals_to_idx.R")
source("R/mk_cpt.R")
library(ggplot2)
library(dplyr)
library(tidyr)

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
                        1, 0, 0,  #   2       1
                        0, 1, 0,  #           2
                        1, 0, 0,  #   3       1
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
                                 1,      0,                   0,  # 2        1        1
                                 1,      0,                   0,  #                   2
                                 1,      0,                   0,  #          2        1
                                 0,      1,                   0,  #                   2
                                 1,      0,                   0,  #          3        1
                                 0, 1 - nu,                  nu,  #                   2
                                 1,      0,                   0,  # 3        1        1
                                 1,      0,                   0,  #                   2
                                 1,      0,                   0,  #          2        1
                                 1,      0,                   0,  #                   2
                                 1,      0,                   0,  #          3        1
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
lambda <- 0.02   #      P(X = present | P_below = C)
#              1 - nu = P(P = low | P_below = med)
nu <- 0.005      #      P(P = med | P_prev = low, P_below = med)
#                     = P(P = high | P_prev = med, P_below = high)
kappa1 <- 0.985  #      P(P = low | P_prev = low, P_below = high)
kappa2 <- 0.014  #      P(P = med | P_prev = low, P_below = high)
# 1 - kappa1 - kappa2 = P(P = high | P_below = high)
thetaX <- 0.80
thetaP <- 0.60
tables <- create_CO2_DBN(N, lambda, nu, kappa1, kappa2, thetaX, thetaP)

# Compute joint distributions given C
q_given_C <- list()
A_given_C <- list()
for (C in 1:tables$q_CPTs$C$ns$prob_of) {
  q_given_C[[C]] <- condition_on_C(tables$q_CPTs, C)
  A_given_C[[C]] <- condition_on_C(tables$A_CPTs, C)
}
B <- with(tables, joint_CPT(B_CPTs))

##########################################################################################
##########################################################################################



## SIMULATE
set.seed(2341 + 544)
TT <- 550  # Number of slices
C <- 2 #sample(3, size = 1)  # Draw capillary threshold from {1, 2, 3} = {Low, Medium, High}
chain <- HMM_sim(q_given_C[[C]]$prob, A_given_C[[C]]$prob, B$prob, TT)

state_vals <- with(chain, idx_to_vals(state, B$ns$given))
colnames(state_vals) <- B$var$given
obs_vals <- with(chain, idx_to_vals(obs, B$ns$prob_of))
colnames(obs_vals) <- B$var$prob_of


# Leakage
tail(state_vals[, 1:N * 2 - 1] - 1)
# obs_vals[, 1:N * 2 - 1] - 1

# Pressure
tail(state_vals[, 1:N * 2 - 2])
# obs_vals[, 1:N * 2 - 2]

y <- chain$obs
# miss <- -seq(1, TT, 30)  # observations every 30th day
miss <- -seq(1, TT, 90)  # observations every 90th day
y[miss] <- NA  
obs_vals[miss, ] <- NA

inf_res <- list()
for (C in 1:tables$q_CPTs$C$ns$prob_of) {
  inf_res[[C]] <- forwards_backwards(
    q = q_given_C[[C]]$prob,
    A = A_given_C[[C]]$prob,
    B = B$prob,
    y = y
  )
}
logLiks <- c(inf_res[[1]]$log_lik, inf_res[[2]]$log_lik, inf_res[[3]]$log_lik)
post_c <- normalize(exp(logLiks - max(logLiks))*tables$q_CPTs$C$prob)
post_c
gamma_mid <- post_c[1]*inf_res[[1]]$smoothed + 
  post_c[2]*inf_res[[2]]$smoothed +
  post_c[3]*inf_res[[3]]$smoothed
margs <- marginal(gamma_mid, B$ns$given)




## TIDY DATA
dfX <- NULL
dfP <- NULL
layer <- 1
for (i in 1:length(margs)) {
  m <- margs[[i]]
  if (i %% 2) {  # CO2 absence/presence
    dfX <- rbind(dfX, data.frame(
      days = 1:TT, absent = m[, 1], present = m[, 2], layer = as.character(layer),
      observed = factor(obs_vals[, i], levels = 1:2, labels = c("absent", "present")),
      actual = factor(state_vals[, i], levels = 1:2, labels = c("absent", "present"))
    ))
  } else {  # Pressure
    dfP <- rbind(dfP, data.frame(
      days = 1:TT, low = m[, 1], medium = m[, 2], high = m[, 3], layer = as.character(layer),
      observed = factor(obs_vals[, i], levels = 1:3, labels = c("low", "medium", "high")),
      actual = factor(state_vals[, i], levels = 1:3, labels = c("low", "medium", "high"))
    ))
    layer <- layer + 1
  }
}
dfXprob <- dfX %>% 
  select(days:layer) %>% 
  pivot_longer(c(absent, present), names_to="category", values_to="probability") %>% 
  mutate(category = factor(category, levels = c("present", "absent")))
dfPprob <- dfP %>% 
  select(days:layer) %>% 
  pivot_longer(c(low, medium, high), names_to="pressure", values_to="probability") %>% 
  mutate(pressure = factor(pressure, levels = c("high", "medium", "low")))

dfXvalues <- dfX %>% 
  select(days, layer, observed, actual) %>% 
  drop_na()
dfPvalues <- dfP %>% 
  select(days, layer, observed, actual) %>% 
  drop_na()

## PLOT
(p1 <- ggplot(dfXprob, aes(days, probability,
                           group = interaction(layer, category), color = layer, linetype = category)
  ) +
  geom_line() +
  theme_bw())

# w <- 12   # cm  instead of 14
# h <- 7    # cm  instead of 8
# ggsave("migProbN4allOBS.pdf",       # ALL OBSERVATIONS, C=3 AND SEED=2341
#   plot = p1, path = "Figures",
#   width = w, height = h, units = "cm"
# )
# ggsave("migProbN4missingOBS.pdf",   # OBSERVATION EVERY 30TH DAY, C=3 AND SEED=2341
#   plot = p1, path = "Figures",
#   width = w, height = h, units = "cm"
# )
# ggsave("migProbN4moremissOBS.pdf",  # OBSERVATION EVERY 90TH DAY, C=2 AND SEED=2341+544
#   plot = p1, path = "Figures",
#   width = w, height = h, units = "cm"
# )


(p2 <- ggplot(dfXvalues, aes(days, observed, group = interaction(1, layer), color = layer)) +
    geom_point(size = 1, position = position_jitter(height = .2)) +
    geom_step(aes(y=actual), direction = "mid", linetype=2) +
    labs(y = "") +
    theme_bw())

# ggsave("migValN4allOBS.pdf",
#   plot = p2, path = "Figures",
#   width = w, height = h, units = "cm"
# )
# ggsave("migValN4missingOBS.pdf",
#   plot = p2, path = "Figures",
#   width = w, height = h, units = "cm"
# )
# ggsave("migValN4moremissOBS.pdf",
#   plot = p2, path = "Figures",
#   width = w, height = h, units = "cm"
# )


(p3 <- ggplot(dfPprob, aes(days, probability,
                           group = interaction(layer, pressure), color = layer, linetype = pressure)
) +
    geom_line() +
    theme_bw())

# ggsave("pressureProbN4allOBS.pdf",
#   plot = p3, path = "Figures",
#   width = w, height = h, units = "cm"
# )
# ggsave("pressureProbN4missingOBS.pdf",
#   plot = p3, path = "Figures",
#   width = w, height = h, units = "cm"
# )
# ggsave("pressureProbN4moremissOBS.pdf",
#        plot = p3, path = "Figures",
#        width = w, height = h, units = "cm"
# )

(p4 <- ggplot(dfPvalues, aes(days, observed, group = interaction(1, layer), color = layer)) +
    geom_point(size = 1, position = position_jitter(height = .2)) +
    geom_step(aes(y=actual), direction = "mid", linetype=2) +
    labs(y = "") +
    theme_bw())

# ggsave("pressureValN4allOBS.pdf",
#        plot = p4, path = "Figures",
#        width = w, height = h, units = "cm"
# )
# ggsave("pressureValN4missingOBS.pdf",
#        plot = p4, path = "Figures",
#        width = w, height = h, units = "cm"
# )
# ggsave("pressureValN4moremissOBS.pdf",
#        plot = p4, path = "Figures",
#        width = w, height = h, units = "cm"
# )








# par(mfrow=c(1, 1))
# plot(margs[[1]][, 2], type = 'l', col=1, ylim=c(-0.1, 1.1), main="P(gas=present) vs true")
# lines(state_vals[, 1] - 1, col=1, lty=2, lwd=2)
# lines(margs[[3]][, 2], col=2)
# lines(state_vals[, 3] - 1, col=2, lty=2, lwd=2)
# lines(margs[[5]][, 2], col=3)
# lines(state_vals[, 5] - 1, col=3, lty=2, lwd=2)
# lines(margs[[7]][, 2], col=4)
# lines(state_vals[, 7] - 1, col=4, lty=2, lwd=2)
# lines(margs[[9]][, 2], col=5)
# lines(state_vals[, 9] - 1, col=5, lty=2, lwd=2)
# 
# 
# par(mfrow=c(2, 1))
# plot(margs[[2]][, 3], type = 'l', col=1, ylim=c(-0.1, 1.1), main="P(pressure=high)")
# lines(margs[[4]][, 3], col=2)
# lines(margs[[6]][, 3], col=3)
# lines(margs[[8]][, 3], col=4)
# plot(state_vals[, 2], type = 'l', col=1, lty=2, lwd=2, ylim=c(0.9, 3.1), main="True pressure")
# lines(state_vals[, 4], col=2, lty=2, lwd=2)
# lines(state_vals[, 6], col=3, lty=2, lwd=2)
# lines(state_vals[, 8], col=4, lty=2, lwd=2)
# 
# 
# 
# # If we knew C:
# margs <- marginal(inf_res[[3]]$smoothed, B$ns$given)
# 
# plot(margs[[1]][, 2], type = 'l', col=1, ylim=c(-0.1, 1.1))
# lines(state_vals[, 1] - 1, col=1, lty=2, lwd=2)
# lines(margs[[3]][, 2], col=2)
# lines(state_vals[, 3] - 1, col=2, lty=2, lwd=2)
# lines(margs[[5]][, 2], col=3)
# lines(state_vals[, 5] - 1, col=3, lty=2, lwd=2)
# lines(margs[[7]][, 2], col=4)
# lines(state_vals[, 7] - 1, col=4, lty=2, lwd=2)
