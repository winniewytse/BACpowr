prior_plot <- function(d_est, d_sd, rho_est, rho_sd, rsq2 = 0) {
  if (d_sd != 0) {
    pd <- ggplot2::ggplot(data.frame(d = c(d_est + 3 * d_sd, d_est - 3 * d_sd)),
                          ggplot2::aes(x = d)) +
      ggplot2::stat_function(fun = dnorm, n = 101,
                             args = list(mean = d_est, sd = d_sd)) +
      ggplot2::labs(x = expression(delta), y = "Density")
  }
  if (rho_sd != 0) {
    shapes <- get_ab(rho_est, rho_sd)
    prho <- ggplot2::ggplot(data.frame(rho = c(rho_est + 4 * rho_sd, 0)),
                          ggplot2::aes(x = rho)) +
      ggplot2::stat_function(fun = dbeta, n = 101,
                             args = list(shape1 = shapes[1], shape2 = shapes[2])) +
      ggplot2::labs(x = expression(rho), y = "Density")
  }
  if (d_sd != 0 & rho_sd == 0) {
    return(pd)
  } else if (d_sd == 0 & rho_sd != 0) {
    return(prho)
  } else if (d_sd != 0 & rho_sd != 0){
    return(gridExtra::grid.arrange(pd, prho, ncol = 2))
  } else {
    return("No prior distributions are specified. ")
  }
}
