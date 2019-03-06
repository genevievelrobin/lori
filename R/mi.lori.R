#' multiple imputation of count data using the lori model
#'
#' @param Y [matrix, data.frame] count table (nxp).
#' @param cov [matrix, data.frame] design matrix (np*q) in order row1xcol1,row2xcol2,..,rownxcol1,row1xcol2,row2xcol2,...,...,rownxcolp
#' @param lambda1 [positive number] the regularization parameter for the interaction matrix.
#' @param lambda2 [positive number] the regularization parameter for the covariate effects.
#' @param M [integer] the number of multiple imputations to perform
#' @param prob [positive number in (0,1)] the proportion of entries to remove in bootstrap
#' @param reff [boolean] whether row effects should be fitted, default value is TRUE
#' @param ceff [boolean] whether column effects should be fitted, default value is TRUE
#' @param rank.max [integer] maximum rank of interaction matrix (smaller than min(n-1,p-1))
#' @param algo type of algorithm to use, either one of "mcgd" (mixed coordinate gradient descent, adapted to large dimensions) or "alt" (alternating minimization, adapted to small dimensions)
#' @param thresh [positive number] convergence tolerance of algorithm, by default \code{1e-6}.
#' @param maxit [integer] maximum allowed number of iterations.
#' @param trace.it [boolean] whether convergence information should be printed
#'
#' @return
#' \item{mi.imputed}{a list of length M containing the imputed count tables}
#' \item{mi.alpha}{a (Mxn) matrix containing in rows the estimated row effects (one row corresponds to one single imputation)}
#' \item{mi.beta}{a (Mxp) matrix containing in rows the estimated column effects (one row corresponds to one single imputation)}
#' \item{mi.epsilon}{a (Mxq) matrix containing in rows the estimated effects of covariates (one row corresponds to one single imputation)}
#' \item{mi.theta}{a list of length M containing the estimated interaction matrices}
#' \item{mi.mu}{a list of length M containing the estimated Poisson means}
#' \item{mi.y}{list of bootstrapped count tables used fot multiple imputation}
#' \item{Y}{original incomplete count table}

#' @export
#'
#' @examples
#' X <- matrix(rnorm(50), 25)
#' Y <- matrix(rpois(25, 1:25), 5)
#' res <- mi.lori(Y, X, 10, 10, 2)
mi.lori <-   function(Y,
                      cov = NULL,
                      lambda1 = NULL,
                      lambda2 = NULL,
                      M = 25,
                      prob = 0.1,
                      reff = T,
                      ceff = T,
                      rank.max = 10,
                      algo = c("alt", "mcgd"),
                      thresh = 1e-5,
                      maxit = 1e3,
                      trace.it = F) {
  if (round(sqrt(M)) < sqrt(M))
    M <- round(sqrt(M)) + 1
  else
    M <- round(sqrt(M))
  ylist <- lapply(1:M, function(k)
    boot.lori(Y))
  reslist <- lapply(1:M, function(k) {
    cat('\r', round(100 * k / (2 * M)), "%", sep = "")
    return(lori(
      Y=ylist[[k]],
      cov=cov,
      lambda1=lambda1,
      lambda2=lambda2,
      reff=reff,
      ceff=ceff,
      rank.max=rank.max,
      algo=algo,
      thresh=thresh,
      maxit=maxit,
      trace.it=trace.it
    ))
  })
  d <- dim(Y)
  n <- d[1]
  p <- d[2]
  mi.imputed <- lapply(reslist, function(res)
    lapply(1:M, function(m)
      matrix(stats::rpois(
        n = n * p, lambda = c(res$imputed)
      ), nrow = n)))
  mi.imputed <- unlist(mi.imputed, recursive = F)
  reslist <- lapply(1:length(mi.imputed), function(m) {
    cat('\r', round(50 + 100 * m / (2 * M ^ 2)), "%", sep = "")
    return(
      lori(
        Y = mi.imputed[[m]],
        cov = cov,
        lambda1 = lambda1,
        lambda2=lambda2,
        reff=reff,
        ceff=ceff,
        rank.max=rank.max,
        algo=algo,
        thresh=thresh,
        maxit=maxit,
        trace.it=trace.it
      )
    )
  })
  mi.alpha <- t(sapply(reslist, function(res)
    res$alpha))
  colnames(mi.alpha) <- rownames(Y)
  mi.beta <- t(sapply(reslist, function(res)
    res$beta))
  colnames(mi.beta) <- colnames(Y)
  mi.epsilon <- t(sapply(reslist, function(res)
    res$epsilon))
  colnames(mi.epsilon) <- colnames(cov)
  mi.theta <- lapply(reslist, function(res)
    res$theta)
  return(
    list(
      mi.imputed = mi.imputed,
      mi.alpha = mi.alpha,
      mi.beta = mi.beta,
      mi.epsilon = mi.epsilon,
      mi.theta = mi.theta,
      mi.y = ylist,
      Y = Y
    )
  )
}

boot.lori <- function(Y) {
  Y <- as.matrix(Y)
  d <- dim(Y)
  n <- d[1]
  p <- d[2]
  m <- sum(!is.na(Y))
  new_counts <-
    stats::rmultinom(1, size = sum(Y, na.rm = T), prob = Y[!is.na(Y)] / m)
  Y[!is.na(Y)] <- new_counts
  return(Y)
}
