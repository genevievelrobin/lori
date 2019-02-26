#' cv.lori
#'
#' @param Y [matrix, data.frame] abundance table (nxp)
#' @param cov [matrix, data.frame] design matris (npxq)
#' @param N [integer] number of cross-validation folds
#' @param thresh [positive number] convergence threshold, default is 1e-5
#' @param maxit [integer] maximum number of iterations, default is 100
#' @param rank.max [integer] maximum rank of interaction matrix, default is 2
#' @param trace.it [boolean] whether information about convergence should be printed
#' @param parallel [boolean] whether the N-fold cross-validation should be parallelized, default value is TRUE
#' @param len [integer] the size of the grid
#' @param reff [boolean] whether row effects should be fitted, default value is TRUE
#' @param ceff [boolean] whether column effects should be fitted, default value is TRUE
#' @param prob [numeric in (0,1)] the proportion of entries to remove for cross-validation
#'
#' @return A list with the following elements
#' \item{lambda1}{regularization parameter estimated by cross-validation for nuclear norm penalty (interaction matrix)}
#' \item{lambda2}{regularization parameter estimated by cross-validation for l1 norm penalty (main effects)}
#' \item{errors}{a table containing the prediction errors for all pairs of parameters}

#' @export
#' @import data.table doParallel parallel corpcor foreach
#'
#' @examples
#' X <- matrix(rnorm(20), 10)
#' Y <- matrix(rpois(10, 1:10), 5)
#' res <- cv.lori(Y, X, N=2, len=2)
cv.lori <- function(Y,
                    cov = NULL,
                    N = 5,
                    thresh = 1e-6,
                    maxit = 100,
                    rank.max = 10,
                    trace.it = F,
                    parallel = F,
                    len = 20,
                    reff=T,
                    ceff=T,
                    prob=0.2) {
  Y <- as.matrix(Y)
  Y2 <- Y
  Y2[is.na(Y2)] <- 0
  m <- sum(!is.na(Y))
  d <- dim(Y)
  n <- d[1]
  p <- d[2]
  alpha <- rep(0, n)
  beta <- rep(0, p)
  epsilon <- rep(0, ncol(cov))
  theta <- matrix(0, n, p)
  gd_main <- grad(Y, cov, alpha, beta, epsilon, theta)
  lambda2.max <- max(abs(gd_main))
  #lambda2.max <-
  #  max(c(rowSums(Y, na.rm = T) / m, colSums(Y, na.rm = T) / m))
  lambda2.min <- max(1e-4,1e-3 * lambda2.max)
  alpmat <- matrix(rep(alpha, p), nrow=n)
  betmat <- matrix(rep(beta, each=n), nrow=n)
  epsmat <- matrix(cov%*% epsilon, nrow=n)
  grad_theta <- (-Y2 + exp(alpmat + betmat + epsmat + theta)) / m
  lambda1.max <- max(svd(grad_theta)$d)
  lambda1.min <- max(1e-3,1e-3*lambda1.max)
  grid.lambda1 <-
    exp(seq(log(lambda1.min), log(lambda1.max), length.out = len))
  grid.lambda2 <-
    exp(seq(log(lambda2.min), log(lambda2.max), length.out = len))
  grid <- as.matrix(data.table::CJ(grid.lambda1, grid.lambda2))
  grid <- grid[nrow(grid):1,]
  na_func <- function(x, prob = 0.2) {
    x <- as.matrix(x)
    omega <- !is.na(x)
    obs.idx <- which(omega)
    yp <- x
    yp[sample(obs.idx, round(prob * sum(omega)))] <- NA
    return(yp)
  }
  if (parallel) {
    nbco <- detectCores()
    cl <- makeCluster(nbco)
    registerDoParallel(cl)
    res.cv <-
      foreach(k = 1:N,
              .packages = c("lori", "corpcor", "parallel", "glmnet")) %dopar% {
                sapply(1:nrow(grid),
                       function(i) {
                         yy <- na_func(as.matrix(Y), prob = prob)
                         if (trace.it)
                           print(paste("lambda", i))
                         res <-
                           lori(as.matrix(yy),
                                cov,
                                lambda1 = grid[i, 1],
                                lambda2 = grid[i, 2],
                                reff=reff,
                                ceff=ceff)$imputed
                         return(sqrt(sum((res - Y) ^ 2, na.rm = T)))
                       })
              }
  } else{
    ylist <-
      lapply(1:N, function(k)
        na_func(as.matrix(Y), prob = prob))
    res.cv <-   lapply(1:N, function(k) {
        sapply(1:nrow(grid),
               function(i) {
                 if (trace.it)
                   cat("\r", "bootstrap sample ", k, "/",N," - ", round(100*i/len^2), "%")
                 res <-
                   lori(ylist[[k]], cov, lambda1 = grid[i, 1], lambda2 = grid[i, 2],
                        reff=reff, ceff=ceff)$imputed
                 return(sqrt(sum((res - Y) ^ 2, na.rm = T)))
               })
    })
  }
  res.cv <- colMeans(do.call(rbind, res.cv))
  l <- which.min(res.cv)
  lambda1 <- grid[l, 1]
  lambda2 <- grid[l, 2]
  dat <-
    data.frame(errors = res.cv,
               lambda1 = grid[, 1],
               lambda2 = grid[, 2])
  return(list(
    lambda1 = lambda1,
    lambda2 = lambda2,
    errors = dat
  ))
}
