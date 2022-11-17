#' Shape and scale hyperparameters of a Gamma Distribution
#'
#' \code{gamma_ab()} computes the shape and scale hyperparameters of a gamma
#' distribution from the mode and standard deviation of a parameter (e.g.,
#' effect size heterogeneity).
#'
#' @param mode Mode of the parameter.
#' @param sd Standard deviation of the parameter.
#' @return Shape and scale hyperparameters of a gamma distribution.
#' @export
#' @examples
#' gamma_ab(.1, .05)
#' gamma_ab(.5, .2)

gamma_ab <- function(mode, sd) {
  bs <- c((mode + sqrt(mode^2 + 4 * sd^2)) / (2 * sd^2),
          (mode - sqrt(mode^2 + 4 * sd^2)) / (2 * sd^2))
  b <- bs[bs > 0]
  a <- b^2 * sd^2
  return(c(a, b))
}
