normalize <- function(x) {
  if (is.vector(x) || ncol(x) == 1) {
    return(x / sum(x))
  } else if (is.matrix(x)) {
    return(x / rowSums(x))
  } else {
    stop("Datatype must be vector or matrix")
  }
}