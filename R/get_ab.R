get_ab <- function(est, se){
  var <- se^2
  a <- var 
  b <- est^3 - est^2 + 7*est*var - 3*var
  c <- - 2*est^3 + est^2 + 16*est^2*var - 14*est*var + 3*var
  d <- 12*est^3*var - 16*est^2*var + 7*est*var - var
  
  alpha <- Re(RConics::cubic(c(a, b, c, d)))
  beta <- (alpha - 1 - est*alpha + 2*est)/est
  ab <- cbind(alpha, beta)
  ab <- ab[ab[, 1] > 1 & ab[, 2] > 1]
  if (identical(a, numeric(0))) 
    stop("The beta distribution is not meaningful with the given shape parameters 
         because SE is too large. ")
  
  return(ab)
}