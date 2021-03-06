% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.lori.R
\name{cv.lori}
\alias{cv.lori}
\title{selection of the regularization
parameters (lambda1 and lambda2) of the lori function by cross-validation}
\usage{
cv.lori(Y, cov = NULL, intercept = F, reff = T, ceff = T,
  rank.max = 10, N = 5, len = 20, prob = 0.2, algo = c("alt",
  "mcgd"), thresh = 1e-06, maxit = 100, trace.it = F, parallel = F)
}
\arguments{
\item{Y}{[matrix, data.frame] abundance table (nxp)}

\item{cov}{[matrix, data.frame] design matris (npxq)}

\item{intercept}{[boolean] whether an intercept should be fitted, default value is FALSE}

\item{reff}{[boolean] whether row effects should be fitted, default value is TRUE}

\item{ceff}{[boolean] whether column effects should be fitted, default value is TRUE}

\item{rank.max}{[integer] maximum rank of interaction matrix, default is 2}

\item{N}{[integer] number of cross-validation folds}

\item{len}{[integer] the size of the grid}

\item{prob}{[numeric in (0,1)] the proportion of entries to remove for cross-validation}

\item{algo}{type of algorithm to use, either one of "mcgd" (mixed coordinate gradient descent, adapted to large dimensions) or "alt" (alternating minimization, adapted to small dimensions)}

\item{thresh}{[positive number] convergence threshold, default is 1e-5}

\item{maxit}{[integer] maximum number of iterations, default is 100}

\item{trace.it}{[boolean] whether information about convergence should be printed}

\item{parallel}{[boolean] whether the N-fold cross-validation should be parallelized, default value is TRUE}
}
\value{
A list with the following elements
\item{lambda1}{regularization parameter estimated by cross-validation for nuclear norm penalty (interaction matrix)}
\item{lambda2}{regularization parameter estimated by cross-validation for l1 norm penalty (main effects)}
\item{errors}{a table containing the prediction errors for all pairs of parameters}
}
\description{
selection of the regularization
parameters (lambda1 and lambda2) of the lori function by cross-validation
}
\examples{
X <- matrix(rnorm(20), 10)
Y <- matrix(rpois(10, 1:10), 5)
res <- cv.lori(Y, X, N=2, len=2)
}
