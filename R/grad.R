grad <- function(Y, cov, mu, alpha, beta, epsilon, theta){
  #internal function to compute gradient of main effects coefficients
  d <- dim(Y)
  n <- d[1]
  p <- d[2]
  q <- ncol(cov)
  grad <- rep(0, n+p+q)
  X <- mu + matrix(as.matrix(cov)%*% epsilon, nrow=n) + theta
  X <- sweep(X, 1, alpha, "+")
  X <- sweep(X, 2, beta, "+")
  grad[1:n] <- rowSums(-Y+exp(X), na.rm=T)
  grad[(n+1):(n+p)] <- colSums(-Y+exp(X), na.rm=T)
  YY <- Y
  YY[is.na(Y)] <- 0
  mat <- exp(X)
  mat[is.na(Y)] <- 0
  grad[(n+p+1):(n+p+q)] <- t(cov)%*%(-c(YY)+c(mat))
  return(grad)
}
