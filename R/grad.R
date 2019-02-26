grad <- function(Y, cov, alpha, beta, epsilon, theta){
  #internal function to compute gradient of main effects coefficients
  d <- dim(Y)
  n <- d[1]
  p <- d[2]
  m <- sum(!is.na(Y))
  Omega <- !is.na(Y)
  q <- ncol(cov)
  grad <- rep(0, n+p+q)
  X <- matrix(rep(alpha, p), nrow=n)
  X <- X + matrix(rep(beta, each=n), nrow=n)
  X <- X + matrix(cov%*% epsilon, nrow=n) + theta
  grad[1:n] <- rowSums(Omega*(-Y+exp(X)), na.rm=T)/m
  grad[(n+1):(n+p)] <- colSums(Omega*(-Y+exp(X)), na.rm=T)/m
  YY <- Y
  YY[is.na(Y)] <- 0
  mat <- exp(X)
  mat[is.na(Y)] <- 0
  grad[(n+p+1):(n+p+q)] <- t(cov)%*%(-c(YY)+c(mat))/m
  return(grad)
}
