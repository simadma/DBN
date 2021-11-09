mk_cpt <- function(table, prob_of, given = NULL, given_ns = NULL) {
  if (is.vector(table)) {
    table <- t(table)  # Vector becomes a matrix with one row
  }
  list(
    var = list(
      prob_of = prob_of,
      given   = given
    ),
    ns = list(
      prob_of = ncol(table),
      given   = given_ns
    ),
    prob = table
  )
}