library(tidyverse)
source("R/normalize.R")
source("R/simulation.R")
source("R/HMM_forwards_backwards.R")
source("R/HMM_learning.R")


q <- c(`no rain` = 0.4, `rain` = 0.6)

A <- rbind(`no rain` = c(`no rain`=0.7, rain=0.3),
           rain      = c(0.2, 0.8))

# q <- normalize(eigen(t(A))$vectors[, 1])

B <- rbind(`no rain` = c(`sad`=0.1, `neutral`=0.3, `happy`=0.6),
           `rain`    = c(0.7, 0.2, 0.1))

set.seed(30)
TT <- 50  # Number of time slices
chain <- HMM_sim(q, A, B, TT)

(p1 <- data.frame(
    day     = 1:TT,
    weather = factor(chain$state, labels = colnames(A)),
    mood    = factor(chain$obs,   labels = colnames(B))
  ) %>%
  pivot_longer(c(weather, mood), names_to = "type") %>% 
  mutate(type = factor(type, levels = c("weather", "mood"))) %>% 
  ggplot(aes(day, value, group = 1)) +
  geom_step(direction = "mid") +
  geom_point(size = 1) +
  facet_wrap(~ type, nrow = 2, scales = "free_y") +
  labs(y = "") +
  theme_bw())

# w <- 14  # cm
# h <- 7   # cm
# ggsave("sim_mood_weather_example.pdf",
#   plot = p1, path = "Figures",
#   width = w, height = h, units = "cm"
# )

## INFERENCE

# Prediction step for h > 0
h <- 0
y <- chain$obs[1:(TT - h)]
inf_res <- forwards_backwards(q, A, B, y)

post <- rbind(inf_res$smoothed,  matrix(0, nrow = h, ncol = 2))
if (h >= 1) for (i in 1:h) {
  post[(TT - h) + i, ] <- t(A) %*% post[(TT - h) + i - 1, ]
}

colnames(post) <- colnames(A)

(p2 <- data.frame(
    day       = 1:TT,
    actual    = chain$state - 1,
    # estimated = inf_res$filtered[, 2]) %>% 
    estimated = post[, "rain"]
  ) %>%
  pivot_longer(c(actual, estimated), names_to = "rain", values_to = "probability") %>% 
  mutate(rain = factor(rain, levels = c("estimated", "actual"))) %>% 
  ggplot(aes(day, probability, linetype = rain)) +
  geom_step(direction = "mid") +
  theme_bw())

# w <- 14  # cm
# h <- 4   # cm
# ggsave("fwdbwd_mood_weather_example.pdf",
#   plot = p2, path = "Figures",
#   width = w, height = h, units = "cm"
# )

## MISSING OBSERVATIONS
y[2*(1:(TT/2))] <- NA  # Every 4th time slice is not observed
# y[1:(TT/2)] <- NA      # First half is not observed
inf_res <- forwards_backwards(q, A, B, y)

post <- inf_res$smoothed
colnames(post) <- colnames(A)

(p3 <- data.frame(
  day       = 1:TT,
  actual    = chain$state - 1,
  # estimated = inf_res$filtered[, 2]) %>% 
  estimated = post[, "rain"]
) %>%
    pivot_longer(c(actual, estimated), names_to = "rain", values_to = "probability") %>% 
    mutate(rain = factor(rain, levels = c("estimated", "actual"))) %>% 
    ggplot(aes(day, probability, linetype = rain)) +
    geom_step(direction = "mid") +
    theme_bw())

# w <- 14  # cm
# h <- 4   # cm
# ggsave("fwdbwd_mood_weather_missing_obs.pdf",
#   plot = p3, path = "Figures",
#   width = w, height = h, units = "cm"
# )


## LEARNING
y <- chain$obs
learn_res <- HMM_learn(y = y, start = q, A = A, B = B, epsilon = 1e-5)
# learn_res <- HMM_learn(y, 2, 3)
inf_res <- forwards_backwards(learn_res$start, learn_res$A, learn_res$B, y)
post <- inf_res$smoothed
colnames(post) <- colnames(A)
(p4 <- data.frame(
  day       = 1:TT,
  actual    = chain$state - 1,
  # estimated = inf_res$filtered[, 2]) %>% 
  estimated = post[, "rain"]
) %>%
    pivot_longer(c(actual, estimated), names_to = "rain", values_to = "probability") %>% 
    mutate(rain = factor(rain, levels = c("estimated", "actual"))) %>% 
    ggplot(aes(day, probability, linetype = rain)) +
    geom_step(direction = "mid") +
    theme_bw())


plot(learn_res$log_lik, type = "l")
learn_res$start

learn_res$A
A

learn_res$B
B



X <- matrix(c(-1, 1, -1, -1), 2)
X_inv <- solve(X)
Lam <- matrix(c(-1, 0, 0, 1), 2)

X%*%Lam%*%X_inv
