plot_Jn <- function(J, n, delta, delta_sd, rho, rho_sd,
                    omega = NULL, omega_sd =  NULL, rsq1 = 0, rsq2 = 0,
                    K = 0, P = .5, power = .8, alpha = .05,
                    test = "two.sided", ep = NULL, al = NULL, smooth = 51) {
  ggplot2::theme_set(ggplot2::theme_bw())
  if (is.null(omega_sd)) {
    args <- list(
      delta = delta, delta_sd = delta_sd, rho = rho, rho_sd = rho_sd,
      rsq2 = rsq2, K = K, P = P, power = power, alpha = alpha, test = test
    )
    if (is.null(al)) {
      fn <- ep_crt2
    } else {
      fn <- al_crt2
    }
  } else {
    args <- list(
      delta = delta, delta_sd = delta_sd, rho = rho, rho_sd = rho_sd,
      omega = omega, omega_sd = omega_sd, rsq1 = rsq1, rsq2 = rsq2,
      K = K, P = P, power = power, alpha = alpha, test = test
    )
    if (is.null(al)) {
      fn <- ep_msrt2
    } else {
      fn <- al_msrt2
    }
  }

  if (is.null(omega_sd)) omega_sd <- 0
  if (is.null(al)) {
    if (delta_sd == 0 & rho_sd == 0 & omega_sd == 0) {
      y <- "Classical Power"
    } else {
      y <- "Expected Power"
    }
    yval <- power
  } else {
    y <- "Assurance Level"
    yval <- al
  }

  pJ <- ggplot2::ggplot(data.frame(J = c(J - J / 2, J + J / 5)), ggplot2::aes(x = J)) +
    ggplot2::stat_function(
      fun = Vectorize(fn, vectorize.args = "J"),
      args = append(args, list(n = n)),
      n = smooth) +
    ggplot2::geom_segment(x = J, xend = J, y = 0, yend = yval,
                          linetype = "dashed", col = "red") +
    ggplot2::geom_segment(x = 0, xend = J, y = yval, yend = yval,
                          linetype = "dashed", col = "red") +
    ggplot2::labs(x = "Number of Clusters (J)", y = y)
  pn <- ggplot2::ggplot(data.frame(n = c(1, n + n / 5)), ggplot2::aes(x = n)) +
    ggplot2::stat_function(
      fun = Vectorize(fn, vectorize.args = "n"),
      args = append(args, list(J = J)),
      n = smooth) +
    ggplot2::geom_segment(x = n, xend = n, y = 0, yend = yval,
                          linetype = "dashed", col = "red") +
    ggplot2::geom_segment(x = 0, xend = n, y = yval, yend = yval,
                          linetype = "dashed", col = "red") +
    ggplot2::labs(x = "Cluster Size (n)", y = y)

  return(list(J = pJ, n = pn))
}
