plot_Jn <- function(J, n, d_est, d_sd, rho_est, rho_sd,
                    rsq2 = 0, K = 0, P = .5, power = .8, alpha = .05,
                    test = "two.sided", ep = NULL, al = NULL) {

  if (is.null(al)) {
    p1 <- ggplot2::ggplot(data.frame(J = c(J - J / 2, J + J / 4)), ggplot2::aes(x = J)) +
      ggplot2::stat_function(
        fun = Vectorize(ep_crt2, vectorize.args = "J"),
        args = list(
          n = n, d_est = d_est, d_sd = d_sd, rho_est = rho_est, rho_sd = rho_sd,
          rsq2 = rsq2, K = K, P = P, power = power, alpha = alpha, test = test
        ),
        n = 31) +
      ggplot2::geom_segment(x = J, xend = J, y = 0, yend = power,
                            linetype = "dashed", col = "red") +
      ggplot2::geom_segment(x = 0, xend = J, y = power, yend = power,
                            linetype = "dashed", col = "red") +
      ggplot2::labs(x = "Number of Clusters (J)", y = "Expected Power")
    p2 <- ggplot2::ggplot(data.frame(n = c(1, n + n/3)), ggplot2::aes(x = n)) +
      ggplot2::stat_function(
        fun = Vectorize(ep_crt2, vectorize.args = "n"),
        args = list(
          J = J, d_est = d_est, d_sd = d_sd, rho_est = rho_est, rho_sd = rho_sd,
          rsq2 = rsq2, K = K, P = P, power = power, alpha = alpha, test = test
        ),
        n = 31) +
      ggplot2::geom_segment(x = n, xend = n, y = 0, yend = power,
                            linetype = "dashed", col = "red") +
      ggplot2::geom_segment(x = 0, xend = n, y = power, yend = power,
                            linetype = "dashed", col = "red") +
      ggplot2::labs(x = "Cluster Size (n)", y = "Expected Power")
  } else {
    p1 <- ggplot2::ggplot(data.frame(J = c(J - J / 2, J + J / 4)), ggplot2::aes(x = J)) +
      ggplot2::stat_function(
        fun = Vectorize(al_crt2, vectorize.args = "J"),
        args = list(
          n = n, d_est = d_est, d_sd = d_sd, rho_est = rho_est, rho_sd = rho_sd,
          rsq2 = rsq2, K = K, P = P, power = power, alpha = alpha, test = test
        ),
        n = 31) +
      ggplot2::geom_segment(x = J, xend = J, y = 0, yend = al,
                            linetype = "dashed", col = "red") +
      ggplot2::geom_segment(x = 0, xend = J, y = al, yend = al,
                            linetype = "dashed", col = "red") +
      ggplot2::labs(x = "Number of Clusters (J)", y = "Assurance Level")
    p2 <- ggplot2::ggplot(data.frame(n = c(1, n + n/3)), ggplot2::aes(x = n)) +
      ggplot2::stat_function(
        fun = Vectorize(al_crt2, vectorize.args = "n"),
        args = list(
          J = J, d_est = d_est, d_sd = d_sd, rho_est = rho_est, rho_sd = rho_sd,
          rsq2 = rsq2, K = K, P = P, power = power, alpha = alpha, test = test
        ),
        n = 31) +
      ggplot2::geom_segment(x = n, xend = n, y = 0, yend = al,
                            linetype = "dashed", col = "red") +
      ggplot2::geom_segment(x = 0, xend = n, y = al, yend = al,
                            linetype = "dashed", col = "red") +
      ggplot2::labs(x = "Cluster Size (n)", y = "Assurance Level")
  }
  return(list(p1, p2))
}
