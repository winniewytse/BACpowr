plot_n <- function(n, d_est, d_sd, power = .8, alpha = .05,
                   test = "two.sided", ep = NULL, al = NULL, smooth = 51) {
  ggplot2::theme_set(ggplot2::theme_bw())

  if (is.null(al)) {
    fn <- ep_2st
    if (d_sd == 0) {
      y <- "Classical Power"
    } else {
      y <- "Expected Power"
    }
    yval <- power
  } else {
    fn <- al_2st
    y <- "Assurance Level"
    yval <- al
  }

  pn <- ggplot2::ggplot(data.frame(n = c(3, n + n / 5)), ggplot2::aes(x = n)) +
    ggplot2::stat_function(
      fun = Vectorize(fn, vectorize.args = "n1"),
      args = list(d_est = d_est, d_sd = d_sd, n2 = NULL,
                  alpha = alpha, power = power, test = test),
      n = smooth) +
    ggplot2::geom_segment(x = n, xend = n, y = 0, yend = yval,
                          linetype = "dashed", col = "red") +
    ggplot2::geom_segment(x = 0, xend = n, y = yval, yend = yval,
                          linetype = "dashed", col = "red") +
    ggplot2::labs(x = "Group Size (n)", y = y)

  return(list(n = pn))
}
