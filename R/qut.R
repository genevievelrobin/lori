#' Automatic selection of nuclear norm regularization parameter
#' @param Y A matrix of counts (contingency table).
#' @param cov A (np)xK matrix of K covariates about rows and columns
#' @param lambda2 A positive number, the regularization parameter for covariates main effects
#' @param q A number between \code{0} and \code{1}. The quantile of the distribution of $lambda_{QUT}$ to take.
#' @param N An integer. The number of parametric bootstrap samples to draw.
#' @return the value of $lambda_{QUT}$ to use in LoRI.
#' @export
#' @examples
#' X = matrix(rnorm(30), 15)
#' Y = matrix(rpois(15, 1:15), 5)
#' lambda = qut(Y,X, 10, N=10)
qut <- function(Y,
                cov,
                lambda2 = 0,
                q = 0.95,
                N = 100,
                reff=T,
                ceff=T) {
  d <- dim(Y)
  n <- d[1]
  p <- d[2]
  nullest <-
    lori(Y, as.matrix(cov), lambda2 = lambda2, lambda1 = 1e5,
         trace.it = F, reff=reff, ceff=ceff)
  X0 <- nullest$X
  lambdas <- rep(0, N)
  m <- sum(!is.na(Y))
  for (i in 1:N)
  {
    Ysimul <- matrix(stats::rpois(n * p, exp(c(X0))), nrow = n)
    nullest2 <- lori(Ysimul, as.matrix(cov), lambda2 = lambda2, lambda1 = 1e5,
                     trace.it = F, reff=reff, ceff=ceff)
    X0simul <- nullest2$X
    dat <-
      sweep(Ysimul - exp(X0simul), 2, colMeans(Ysimul - exp(X0simul)))
    dat <-
      sweep(Ysimul - exp(X0simul), 1, rowMeans(Ysimul - exp(X0simul)))
    lambdas[i] <-
      (1 / m) * svd::propack.svd(dat, neig = 1, opts = list(maxiter = 1e5))$d
    cat('\r', i, "/", N)
  }
  return(stats::quantile(lambdas, q)[[1]])
}
