source("R/normalize.R")
source("R/roll_unroll.R")
source("R/idx_to_vals_and_vals_to_idx.R")
source("R/mk_cpt.R")

income <- mk_cpt(
  table   = c(0.2631578947368421, 0.4736842105263158, 0.2631578947368421),
  prob_of = "income"
)

ratio_of_debts_to_income <- mk_cpt(
  table   = c(0.3846153846153846, 0.6153846153846154),
  prob_of = "ratio_of_debts_to_income"
)


age <- mk_cpt(
  table   = c(0.14285714285714285, 0.5714285714285714, 0.2857142857142857),
  prob_of = "age"
)

assets <- mk_cpt(
  table    = rbind(
    c(0.6, 0.3, 0.09999999999999999),
    c(0.2, 0.6, 0.2),
    c(0.1, 0.3, 0.6)
  ),
  prob_of  = "assets",
  given    = "income",
  given_ns = 3
)

payment_history <- mk_cpt(
  table    = rbind(
    c(0.5, 0.3, 0.2),
    c(0.1, 0.4, 0.5),
    c(0.6, 0.3, 0.10000000000000009),
    c(0.2, 0.3, 0.5),
    c(0.9, 0.07, 0.03),
    c(0.3, 0.32, 0.38)
  ),
  prob_of  = "payment_history",
  given    = c("age", "ratio_of_debts_to_income"),
  given_ns = c(3, 2)
)

future_income <- mk_cpt(
  table    = rbind(
    c(0.95, 0.05),
    c(0.7, 0.3),
    c(0.5, 0.5),
    c(0.7, 0.3),
    c(0.5, 0.5),
    c(0.3, 0.7),
    c(0.5, 0.5),
    c(0.3, 0.7),
    c(0.05, 0.95)
  ),
  prob_of  = "future_income",
  given    = c("income", "assets"),
  given_ns = c(3, 3)
)

reliability <- mk_cpt(
  table    = rbind(
    c(0.7, 0.30000000000000004),
    c(0.8, 0.19999999999999996),
    c(0.9, 0.09999999999999998),
    c(0.3, 0.7),
    c(0.5, 0.5),
    c(0.6, 0.4),
    c(0.5, 0.5),
    c(0.5, 0.5),
    c(0.5, 0.5)
  ),
  prob_of  = "reliability",
  given    = c("payment_history", "age"),
  given_ns = c(3, 3)
)

credit_worthiness <- mk_cpt(
  table    = rbind(
    c(0.99, 0.01),
    c(0.85, 0.15),
    c(0.75, 0.25),
    c(0.5, 0.5),
    c(0.5, 0.5),
    c(0.25, 0.75),
    c(0.15, 0.85),
    c(0.01, 0.99)
  ),
  prob_of  = "credit_worthiness",
  given    = c("reliability", "future_income", "ratio_of_debts_to_income"),
  given_ns = c(2, 2, 2)
)

joint_prob <- joint_CPT(list(
  #income, ratio_of_debts_to_income, age,
  assets, payment_history,
  future_income, reliability,
  credit_worthiness
))


margs <- with(joint_prob, marginal(prob, ns$prob_of))
names(margs) <- joint_prob$var$prob_of


#######################################3
# TESTING
cpt_res <- cpt_from_joint_cpt(joint_prob, "payment_history", c("age", "ratio_of_debts_to_income"))
cpt_res
