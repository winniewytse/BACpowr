crtJn <- function(d_est, d_se, rho_est, rho_se, r2_est, r2_se, 
                  J = NULL, n = NULL, K = 0, power = .80, 
                  test = "two-tailed", plot = FALSE, 
                  abs.tol = 1e-10, x.tol = 1.5e-15, 
                  rel.tol = 1e-15, sing.tol = 1e-20){
  ggplot2::theme_set(ggplot2::theme_bw())
  lossJ <- function(J) {
    sum((ep_crt(J = J, n = n, d_est = d_est, d_se = d_se,
                    rho_est = rho_est, rho_se = rho_se, 
                    r2_est = r2_est, r2_se = r2_se, 
                    test = test) - power)^2)
  }
  lossn <- function(n) {
    sum((ep_crt(J = J, n = n, d_est = d_est, d_se = d_se,
                    rho_est = rho_est, rho_se = rho_se, 
                    r2_est = r2_est, r2_se = r2_se, 
                    test = test) - power)^2)
  }
  if (!is.null(n)) {
    J <- nlminb(start = 4, lossJ, lower = 1,
                control = list(abs.tol = abs.tol, x.tol = x.tol, 
                               rel.tol = rel.tol, sing.tol = sing.tol))$par
  } else if (!is.null(J)) {
    n <- nlminb(start = 0, lossn, lower = 1,
                control = list(abs.tol = abs.tol, x.tol = x.tol, 
                               rel.tol = rel.tol, sing.tol = sing.tol))$par
  } else {
    n <- 1e10
    J <- nlminb(start = c(4), lossJ, lower = c(1),
                control = list(abs.tol = abs.tol, x.tol = x.tol, 
                               rel.tol = rel.tol, sing.tol = sing.tol))$par
    rm(n)
    n <- nlminb(start = 0, lossn, lower = 1,
                control = list(abs.tol = abs.tol, x.tol = x.tol, 
                               rel.tol = rel.tol, sing.tol = sing.tol))$par
  }
  if (plot) {
    p1 <- ggplot2::ggplot(data.frame(J = c(4, J + J/3)), aes(x = J)) + 
      stat_function(fun = Vectorize(ep_crt, vectorize.args = "J"), 
                    args = list(n = n, r2_est = r2_est, r2_se = r2_se, 
                                d_est = d_est, d_se = d_se, 
                                rho_est = rho_est, rho_se = rho_se), 
                    n = 51) + 
      geom_segment(x = J, xend = J, y = 0, yend = power, 
                   linetype = "dashed", col = "red") +
      geom_segment(x = 0, xend = J, y = power, yend = power, 
                   linetype = "dashed", col = "red") +
      labs(x = "Number of Clusters (J)", y = "Generalized Power")
    p2 <- ggplot2::ggplot(data.frame(n = c(1, n + n/3)), aes(x = n)) + 
      stat_function(fun = Vectorize(ep_crt, vectorize.args = "n"), 
                    args = list(J = J, r2_est = r2_est, r2_se = r2_se, 
                                d_est = d_est, d_se = d_se, 
                                rho_est = rho_est, rho_se = rho_se), 
                    n = 51) + 
      geom_segment(x = n, xend = n, y = 0, yend = power, 
                   linetype = "dashed", col = "red") +
      geom_segment(x = 0, xend = n, y = power, yend = power, 
                   linetype = "dashed", col = "red") +
      labs(x = "Cluster Size (n)", y = "Generalized Power")
    p <- gridExtra::grid.arrange(p1, p2, nrow = 2)
    return(list(p, ceiling(cbind(J = J, n = n))))
  } else {
    return(ceiling(cbind(J = J, n = n)))
  }
}