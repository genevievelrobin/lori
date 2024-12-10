#' The lori method implements a method to analyze and impute
#' incomplete count tables. An important feature of the method
#' is that it can take into account main effects of rows and
#' columns, as well as effects of continuous or categorical
#' covariates, and interaction. The estimation procedure is
#' based on minimizing a Poisson loss penalized by a Lasso
#' type penalty (sparse vector of covariate effects) and a
#' nuclear norm penalty inducing a low-rank interaction matrix
#' (a few latent factors summarize the interactions).
#'
#'
#' @param Y [matrix, data.frame] count table (nxp).
#' @param cov [matrix, data.frame] design matrix (np*q) in order row1xcol1,row2xcol2,..,rownxcol1,row1xcol2,row2xcol2,...,...,rownxcolp
#' @param lambda1 [positive number] the regularization parameter for the interaction matrix.
#' @param lambda2 [positive number] the regularization parameter for the covariate effects.
#' @param intercept [boolean] whether an intercept should be fitted, default value is FALSE
#' @param reff [boolean] whether row effects should be fitted, default value is TRUE
#' @param ceff [boolean] whether column effects should be fitted, default value is TRUE
#' @param rank.max [integer] maximum rank of interaction matrix (smaller than min(n-1,p-1))
#' @param algo type of algorithm to use, either one of "mcgd" (mixed coordinate gradient descent, adapted to large dimensions) or "alt" (alternating minimization, adapted to small dimensions)
#' @param thresh [positive number] convergence tolerance of algorithm, by default \code{1e-6}.
#' @param maxit [integer] maximum allowed number of iterations.
#' @param trace.it [boolean] whether convergence information should be printed
#' @param parallel [boolean] whether computations should be performed in parallel on multiple cores
#' @import rARPACK svd
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
#' \item{nparam}{number of estimated parameters in the model}
#' \item{dfres}{residual degrees of freedom}
#' \item{chisq}{sum of squared deviations between observed and expected counts normalized by the expected value}
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
           intercept = T,
           reff = T,
           ceff = T,
           rank.max = 2,
           algo = c("alt","mcgd"),
           thresh = 1e-5,
           maxit = 100,
           trace.it = F,
           parallel=F)
  {
    Y <- as.matrix(Y)
    d <- dim(Y)
    n <- d[1]
    p <- d[2]
    if (!is.null(cov)){
      cov <- as.matrix(cov)
      q <- ncol(cov)
    }
    else
      q <- 0
    rank.max = min(n - 1, p - 1, rank.max)
    algo <- match.arg(algo,c("alt","mcgd"),several.ok=T)[1]
    if(min(n,p)<3) algo <- "alt"
    if(algo=="mcgd"){
      res <- mcgd(
        Y,
        lambda1,
        lambda2,
        cov = cov,
        thresh = thresh,
        trace.it = trace.it,
        rank.max = rank.max,
        intercept = intercept,
        reff = reff,
        ceff = ceff
      )
    } else{
      res <- altmin(
        Y,
        lambda1,
        lambda2,
        cov = cov,
        thresh = thresh,
        trace.it = trace.it,
        rank.max = rank.max,
        intercept = intercept,
        reff = reff,
        ceff = ceff,
        parallel=F
      )
    }
    X <- res$X
    mu <- res$mu
    alpha <- res$alpha
    beta <- res$beta
    epsilon <- res$epsilon
    theta <- res$theta
    imputed <- round(exp(X))
    imputed[!is.na(Y)] <- Y[!is.na(Y)]
    means <- exp(X)
    nparam <- sum(abs(alpha)>0) + sum(abs(beta)>0) + sum(abs(epsilon)>0)
    r <- min(sum(abs(svd(theta)$d)>0), rank.max)
    nparam <- nparam + r*(nrow(Y) + ncol(Y) - r)
    dfres <- nrow(Y)*ncol(Y) - nparam
    chisq <- sum((Y - means)^2/means, na.rm=T)
    res2 <-
      structure(
        list(
          X = X,
          mu=mu,
          alpha = alpha,
          beta = beta,
          epsilon = epsilon,
          theta = theta,
          imputed = imputed,
          means = means,
          cov = cov,
          objective=res$objective,
          nparam=nparam,
          dfres=dfres,
          chisq=chisq
        )
      )
    class(res2) <- c("lori", "list")
    return(res2)
  }
