#' Shape Parameters of a Beta Distribution
#'
#' \code{beta_ab()} computes the shape hyperparameters of a beta distribution
#' from the mode and standard deviation of a parameter (e.g., intraclass correlation).
#'
#' @param mode Mode of the parameter.
#' @param sd Standard deviation of the parameter.
#' @return Shape hyperparameters of a beta distribution.
#' @export
#' @examples
#' beta_ab(.1, .05)
#' beta_ab(.5, .2)
beta_ab <- function(mode, sd){
  var <- sd^2
  a <- var
  b <- mode^3 - mode^2 + 7*mode*var - 3*var
  c <- - 2*mode^3 + mode^2 + 16*mode^2*var - 14*mode*var + 3*var
  d <- 12*mode^3*var - 16*mode^2*var + 7*mode*var - var
  alpha <- Re(RConics::cubic(c(a, b, c, d)))
  beta <- (alpha - 1 - mode*alpha + 2*mode)/mode
  ab <- cbind(alpha, beta)
  ab <- ab[ab[, 1] > 1 & ab[, 2] > 1]
  if (identical(ab, numeric(0)))
    stop("With the given mode and standard deviation, the beta distribution is
         U-shaped, which may not serve as a reasonable prior.")

  return(ab)
}
