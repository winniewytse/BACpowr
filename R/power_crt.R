power_crt <- function(J, n, rho, r2, delta, 
                      test = "two-tailed", K = 0) {
  df <- J - K - 2
  if (test == "two-tailed") {
    cv <- qt(.975, df)
    ncp <- delta * sqrt(J * n / 4 / (1 + (n * (1 - r2) - 1) * rho))
    pow <- pt(cv, df = df, ncp = ncp, lower.tail = FALSE) + 
      pt(-cv, df = df, ncp = ncp, lower.tail = TRUE)
  } else if (test == "one-tailed") {
    cv <- qt(.95, df)
    ncp <- delta * sqrt(J * n / 4 / (1 + (n * (1 - r2) - 1) * rho))
    pow <- pt(cv, df = df, ncp = ncp, lower.tail = FALSE)
  }
  return(pow)
}