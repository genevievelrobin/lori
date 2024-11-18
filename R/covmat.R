#' covmat
#'
#' @param n number of rows
#' @param p number ofcolumns
#' @param R nxK1 matrix of row covariates
#' @param C nxK2 matrix of column covariates
#' @param E (n+p)xK3 matrix of row-column covariates
#' @param center boolean indicating whether the returned covariate matrix should be centered (for identifiability)
#'
#' @return the joint product of R and C column-binded with E, a (np)x(K1+K2+K3) matrix in order row1col1,row2col1,...,rowncol1, row1col2, row2col2,...,rowncolp
#' @export
#'
#' @examples
#' R <- matrix(rnorm(10), 5)
#' C <- matrix(rnorm(9), 3)
#' covs <- covmat(5,3,R,C)
covmat <- function(n, p, R=NULL, C=NULL, E=NULL, center=F) {
  if(center){
    if(!is.null(R)) R <- scale(R, scale=F)
    if(!is.null(C)) C <- scale(C, scale=F)
    if(!is.null(E)){
      Eprime <- E
      for(k in 1:ncol(E)){
        if(!is.factor(E[,k])){
          EE <- matrix(E[,k], nrow=n)
          EE <- scale(EE, scale=F)
          EE <- t(scale(t(EE), scale=F))
          Eprime[,k] <- c(EE)
        }
      }
      E <- Eprime
    }
  }
  if (is.null(R)) {
    covs <- covmatC(C, n)
  } else if (is.null(C)) {
    covs <- covmatR(R, p)
  } else {
    R <- as.matrix(R)
    C <- as.matrix(C)
    dR <- dim(R)
    dC <- dim(C)
    K1 <- dR[2]
    K2 <- dC[2]
    covs <-
      cbind(do.call(rbind, replicate(nrow(C), R, simplify = FALSE)),
            C[rep(seq_len(nrow(C)), each = nrow(R)), ])
  }
  if(!is.null(E)){
    covs <- cbind(covs, E)
  }
  return(covs)
}
covmatR <- function(R, p) {
  R <- as.matrix(R)
  dR <- dim(R)
  n <- dR[1]
  K1 <- dR[2]
  covs <- do.call(rbind, replicate(p, R, simplify = FALSE))
  return(covs)
}

covmatC <- function(C, n) {
  C <- as.matrix(C)
  dC <- dim(C)
  p <- dC[1]
  K2 <- dC[2]
  covs <- C[rep(seq_len(p), each = n), ]
  return(covs)
}




