#' LORI method.
#'
#' @param Y [matrix, data.frame] count table (nxp).
#' @param cov [matrix, data.frame] design matrix (np*q) in order row1xcol1,row2xcol2,..,rownxcol1,row1xcol2,row2xcol2,...,...,rownxcolp
#' @param lambda1 [positive number] the regularization parameter for the interaction matrix.
#' @param lambda2 [positive number] the regularization parameter for the covariate effects.
#' @param reff [boolean] whether row effects should be fitted, default value is TRUE
#' @param ceff [boolean] whether column effects should be fitted, default value is TRUE
#' @param rank.max [integer] maximum rank of interaction matrix (smaller than min(n-1,p-1))
#' @param thresh [positive number] convergence tolerance of algorithm, by default \code{1e-6}.
#' @param maxit [integer] maximum allowed number of iterations.
#' @param trace.it [boolean] whether convergence information should be printed
#' @import glmnet stats
#' @export
#' @return A list with the following elements
#' \item{X}{nxp matrix of log of expected counts}
#' \item{alpha}{row effects}
#' \item{beta}{column effects}
#' \item{epsilon}{covariate effects}
#' \item{theta}{nxp matrix of row-column interactions}
#' \item{imputed}{nxp matrix of imputed counts}
#' \item{means}{nxp matrix of expected counts (exp(X))}
#' \item{cov}{npxK matrix of covariates}
#' @examples
#' \dontshow{
#' X <- matrix(rnorm(50), 25)
#' Y <- matrix(rpois(25, 1:25), 5)
#' res <- lori(Y, X, 10, 10)
#' }
lori <-
  function(Y,
           cov = NULL,
           lambda1 = NULL,
           lambda2 = NULL,
           reff = T,
           ceff = T,
           rank.max = 10,
           thresh = 1e-5,
           maxit = 1e3,
           trace.it = F)
  {
    Y <- as.matrix(Y)
    d <- dim(Y)
    n <- d[1]
    p <- d[2]
    if (!is.null(cov))
      q <- ncol(cov)
    else
      q <- 0
    rank.max = min(n - 1, p - 1, rank.max)
    #lambda2 <- rep(lambda2,n+p+q)
    #covsc <- cov[, apply(cov, 2, is.factor)==F]
    #sc <- mean(apply(covsc, 2, sd))
    #if(is.na(sc)) sc <- 1
    #if(is.null(sc)) sc <- 1
    #if(sc==0) sc <- 1
    #lambda2[1:n] <- lambda2[1:n]*sqrt(p*(1-1/n)^2/(n*p))/sc
    #lambda2[(n+1):(n+p)] <- lambda2[(n+1):(n+p)]*sqrt(n*(1-1/p)^2/(n*p))/sc
    res <- altmin(
      Y,
      lambda1,
      lambda2,
      cov = cov,
      thresh = thresh,
      trace.it = trace.it,
      rank.max = rank.max,
      reff = reff,
      ceff = ceff
    )
    X <- res$X
    alpha <- res$alpha
    beta <- res$beta
    epsilon <- res$epsilon
    theta <- res$theta
    #imputed <- matrix(rpois(n*p, lambda=c(exp(X))), nrow=n)
    imputed <- exp(X)
    imputed[!is.na(Y)] <- Y[!is.na(Y)]
    means <- exp(X)
    res2 <-
      structure(
        list(
          X = X,
          alpha = alpha,
          beta = beta,
          epsilon = epsilon,
          theta = theta,
          imputed = imputed,
          means = means,
          cov = cov,
          objective=res$objective
        )
      )
    return(res2)
  }


