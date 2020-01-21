#' @title Compute the log sum of exponentials.
#'
#' @description The reference for the approach is in the Stan Manual. The computation is for the log sum.
#'
#' @param x Data.
#' @return The log sum of x.
get_log_sum_exponential <- function(x) {
  maxVal <- max(x)
  log(maxVal) + log(sum(exp(log(x) - log(maxVal))))
}
