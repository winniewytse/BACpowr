#' Shape Parameters of a Beta Distribution
#'
#' \code{get_ab()} computes the shape parameters of a beta distribution from the mean
#' and standard deviation of a parameter.
#'
#' @param est Intraclass correlation estimate.
#' @param sd Uncertainty in the intraclass correlation estimate.
#' @return Shape paraemters for a beta distribution given the mean and standard deviation.
#' @export
#' @examples
#' get_ab(.1, .05)
#' get_ab(.5, .2)
get_ab <- function(est, sd){
  var <- sd^2
  a <- var
  b <- est^3 - est^2 + 7*est*var - 3*var
  c <- - 2*est^3 + est^2 + 16*est^2*var - 14*est*var + 3*var
  d <- 12*est^3*var - 16*est^2*var + 7*est*var - var

  alpha <- Re(RConics::cubic(c(a, b, c, d)))
  beta <- (alpha - 1 - est*alpha + 2*est)/est
  ab <- cbind(alpha, beta)
  ab <- ab[ab[, 1] > 1 & ab[, 2] > 1]
  if (identical(ab, numeric(0)))
    stop("The beta distribution is not meaningful with the given shape parameters
  because the standard deviation is too large. ")

  return(ab)
}
