#' The cv.lori method performs automatic selection of the
#' regularization parameters (lambda1 and lambda2) used in
#' the lori function. These parameters are selected by
#' cross-validation. The classical procedure is to apply
#' cv.lori to the data to select the regularization parameters,
#' and to then impute and analyze the data using the lori
#' function (or mi.lori for multiple imputation).
#'
#' @param Y [matrix, data.frame] abundance table (nxp)
#' @param cov [matrix, data.frame] design matris (npxq)
#' @param intercept [boolean] whether an intercept should be fitted, default value is FALSE
#' @param reff [boolean] whether row effects should be fitted, default value is TRUE
#' @param ceff [boolean] whether column effects should be fitted, default value is TRUE
#' @param rank.max [integer] maximum rank of interaction matrix, default is 2
#' @param N [integer] number of cross-validation folds
#' @param len [integer] the size of the grid
#' @param prob [numeric in (0,1)] the proportion of entries to remove for cross-validation
#' @param algo type of algorithm to use, either one of "mcgd" (mixed coordinate gradient descent, adapted to large dimensions) or "alt" (alternating minimization, adapted to small dimensions)
#' @param thresh [positive number] convergence threshold, default is 1e-5
#' @param maxit [integer] maximum number of iterations, default is 100
#' @param trace.it [boolean] whether information about convergence should be printed
#' @param parallel [boolean] whether computations should be performed in parallel on multiple cores
#'
#' @return A list with the following elements
#' \item{lambda1}{regularization parameter estimated by cross-validation for nuclear norm penalty (interaction matrix)}
#' \item{lambda2}{regularization parameter estimated by cross-validation for l1 norm penalty (main effects)}
#' \item{errors}{a table containing the prediction errors for all pairs of parameters}

#' @export
#' @import data.table parallel
#'
#' @examples
#' X <- matrix(rnorm(20), 10)
#' Y <- matrix(rpois(10, 1:10), 5)
#' res <- cv.lori(Y, X, N=2, len=2)
cv.lori <- function(Y,
                    cov = NULL,
                    intercept = T,
                    reff = T,
                    ceff = T,
                    rank.max = 5,
                    N = 5,
                    len = 20,
                    prob = 0.2,
                    algo = c("alt", "mcgd"),
                    thresh = 1e-5,
                    maxit = 10,
                    trace.it = F,
                    parallel=F) {
  Y <- as.matrix(Y)
  Y2 <- Y
  Y2[is.na(Y2)] <- 0
  m <- sum(!is.na(Y))
  d <- dim(Y)
  n <- d[1]
  p <- d[2]
  mu <- 0
  alpha <- rep(0, n)
  beta <- rep(0, p)
  epsilon <- rep(0, ncol(cov))
  alpmat <- matrix(rep(alpha, p), nrow = n)
  betmat <- matrix(rep(beta, each = n), nrow = n)
  epsmat <- matrix(as.matrix(cov) %*% epsilon, nrow = n)
  theta <- matrix(0, n, p)
  X <- alpmat + betmat + epsmat + theta
  gd_main <- grad(Y, as.matrix(cov), mu, alpha, beta, epsilon, theta)
  lambda2.max <- max(abs(gd_main))
  lambda2.min <- max(1e-4, 1e-3 * lambda2.max)
  grad_theta <- (-Y2 + exp(alpmat + betmat + epsmat + theta))
  lambda1.max <- max(svd(grad_theta)$d)
  lambda1.min <- max(1e-3, 1e-3 * lambda1.max)
  grid.lambda1 <-
    exp(seq(log(lambda1.min), log(lambda1.max), length.out = len))
  grid.lambda2 <-
    exp(seq(log(lambda2.min), log(lambda2.max), length.out = len))
  grid <- as.matrix(data.table::CJ(grid.lambda1, grid.lambda2))
  grid <- grid[nrow(grid):1, ]
  na_func <- function(x, prob=0.1){
    x <- as.matrix(x)
    yp <- x
    d <- dim(x)
    n <- d[1]
    p <- d[2]
    n_na <- round(prob*sum(!is.na(x)))
    idx <- sample(which(!is.na(x)), size=n_na, replace = FALSE)
    x[idx] <- NA
    return(x)
  }
  if (parallel) n_cores <- detectCores() else n_cores <- 1
  ylist <-
    mclapply(1:N, function(k)
      na_func(as.matrix(Y), prob = prob), mc.cores=n_cores)
  res.cv <-   mclapply(1:N, function(k) {
    sapply(1:nrow(grid),
           function(i) {
             if (trace.it)
               cat("\r",
                   "bootstrap sample ",
                   k,
                   "/",
                   N,
                   " - ",
                   round(100 * i / len ^ 2),
                   "%")
             res <-
               lori(
                 Y=ylist[[k]],
                 cov=cov,
                 lambda1 = grid[i, 1],
                 lambda2 = grid[i, 2],
                 intercept=intercept,
                 reff = reff,
                 ceff = ceff,
                 rank.max=rank.max,
                 algo=algo,
                 thresh=thresh,
                 maxit=maxit,
                 trace.it=trace.it
               )$imputed
             return(sqrt(sum((res - Y) ^ 2, na.rm = T)))
           })
  }, mc.cores = n_cores)
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
