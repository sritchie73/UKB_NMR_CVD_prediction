logit <- function(value) {
  log( value / (1 - value) )
}

invlogit <- function(value) {
  1 / (1 + exp(-value))
}
