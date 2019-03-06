#' aggregate lori multiple imputation results
#'
#' @param res.mi a multiple imputation result from the function mi.lori
#'
#' @return
#' \item{pool.impute}{a list containing the pooled means (mean) and variance (var) of the imputed values}
#' \item{pool.alpha}{a list containing the pooled means (mean) and variance (var) of the row effects}
#' \item{pool.beta}{a list containing the pooled means (mean) and variance (var) of the column effects}
#' \item{pool.epsilon}{a list containing the pooled means (mean) and variance (var) of the covariate effects}
#' \item{pool.theta}{a list containing the pooled means (mean) and variance (var) of the interactions}
#' @export
#'
#' @examples
 #' X <- matrix(rnorm(50), 25)
#' Y <- matrix(rpois(25, 1:25), 5)
#' res <- mi.lori(Y, X, 10, 10, 2)
#' poolres <- pool.lori(res)
pool.lori <- function(res.mi){
  M <- length(res.mi$mi.imputed)
  impmean <- Reduce("+", res.mi$mi.imputed)/M
  var_between <- lapply(res.mi$mi.imputed, function(res) (res-impmean)^2)
  var_between <- (M+1)/(M*(M-1))*Reduce("+", var_between)
  resid <- lapply(1:M, function(m) (res.mi$mi.imputed[[m]]-res.mi$mi.means[[m]])^2/res.mi$mi.means[[m]])
  resid <- Reduce("+", resid)/M
  impmean2 <- impmean
  impmean2[!is.na(res.mi$Y)] <- res.mi$Y[!is.na(res.mi$Y)]
  var_imp <- resid+var_between
  var_imp[!is.na(res.mi$Y)] <- 0
  pool.impute <- list(mean=impmean2, var=var_imp)
  pool.alpha <- list(mean=colMeans(res.mi$mi.alpha),
                     var=apply(res.mi$mi.alpha, 2, var))
  pool.beta <- list(mean=colMeans(res.mi$mi.beta),
                     var=apply(res.mi$mi.beta, 2, var))
  pool.epsilon <- list(mean=colMeans(res.mi$mi.epsilon),
                     var=apply(res.mi$mi.epsilon, 2, var))
  meantheta <- Reduce("+", res.mi$mi.theta)/M
  pool.theta <- list(mean = meantheta,
                     var = Reduce("+", lapply(res.mi$mi.theta, function(t) (t-meantheta)^2))/(M-1))
  return(list(pool.impute=pool.impute, pool.alpha=pool.alpha,
              pool.beta=pool.beta, pool.epsilon=pool.epsilon,
              pool.theta=pool.theta))
}
