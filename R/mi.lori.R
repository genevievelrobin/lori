#' Multiple imputation of count data using side information
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
                      M = 100,
                      prob=0.1,
                      reff = T,
                      ceff = T,
                      rank.max = 10,
                      thresh = 1e-5,
                      maxit = 1e3,
                      trace.it = F){

  ylist <- lapply(1:M, function(k) boot.lori(Y, prob))
  reslist <- lapply(1:M, function(k){
    cat('\r', round(100*k/M), "%", sep="")
    return(lori(ylist[[k]], cov, lambda1, lambda2, reff, ceff, rank.max,
                thresh, maxit, trace.it))
  })
  mi.imputed <- lapply(reslist, function(res) res$imputed)
  mi.means <- lapply(reslist, function(res) res$means)
  mi.alpha <- t(sapply(reslist, function(res) res$alpha))
  colnames(mi.alpha) <- rownames(Y)
  mi.beta <- t(sapply(reslist, function(res) res$beta))
  colnames(mi.beta) <- colnames(Y)
  mi.epsilon <- t(sapply(reslist, function(res) res$epsilon))
  colnames(mi.epsilon) <- colnames(cov)
  mi.theta <- lapply(reslist, function(res) res$theta)
  return(list(mi.imputed=mi.imputed, mi.means=mi.means,
              mi.alpha=mi.alpha, mi.beta=mi.beta,
              mi.epsilon=mi.epsilon, mi.theta=mi.theta,
              mi.y = ylist, Y=Y))
}

boot.lori <- function(Y, prob=0.1){
  Y <- as.matrix(Y)
  d <- dim(Y)
  n <- d[1]
  p <- d[2]
  n_na <- round(prob*sum(!is.na(Y)))
  n_na_sites <- round(n_na/n)
  n_na_years <- round(n_na/p)
  idx_sites <- sample(which(rowSums(!is.na(Y))>n_na_sites+1))
  idx_years <- sample(which(colSums(!is.na(Y))>n_na_years+1))
  yp <- Y[idx_sites,idx_years]
  yp[sample(which(!is.na(yp)), n_na)] <- NA
  Y[idx_sites,idx_years] <- yp
  Y[which(!is.na(Y))] <- stats::rpois(sum(!is.na(Y)), c(Y[which(!is.na(Y))]))
  return(Y)
}
