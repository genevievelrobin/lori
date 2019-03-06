mcgd <- function(Y,
                 lambda1,
                 lambda2,
                 cov = NULL,
                 rank.max = 10,
                 thresh = 1e-6,
                 maxit = 1e3,
                 trace.it = F,
                 intercept = F,
                 reff = T,
                 ceff = T) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  q <- ncol(cov)
  rank.max <- min(n - 1, p - 1, rank.max)
  Omega <- 1 - 1 * is.na(Y)
  m <- sum(!is.na(Y))
  mu <- 0
  theta <- matrix(0, nrow = n, ncol = p)
  alpha <- rep(0, n)
  beta <- rep(0, p)
  R <- 0
  if (!is.null(cov)) {
    epsilon <- rep(0, ncol(cov))
  } else{
    epsilon <- 0
  }
  alpmat <- matrix(rep(alpha, p), nrow = n)
  betmat <- matrix(rep(beta, each = n), nrow = n)
  if(is.null(cov)){
    epsmat <- matrix(cov %*% epsilon, nrow = n)
  } else epsmat <- 0
  error = 1
  objective <- NULL
  count = 1
  Y2 <- Y
  Y2[is.na(Y)] <- 0
  param <- mu + alpmat + betmat + epsmat + theta
  Ytalp <- c(Y[!is.na(Y)])
  low_bound <- 0
  low_bound <-
    low_bound - sum(Y[!(Y == 0)] * (1 - log(Y[!(Y == 0)])), na.rm = T)
  obj0 <- low_bound
  obj0 <- obj0 + sum(-(Y * param) + exp(param), na.rm = T)
  U <- obj0 / lambda1
  while ((error > thresh) && (count < maxit)) {
    mu.tmp <- mu
    alpha.tmp <- alpha
    beta.tmp <- beta
    epsilon.tmp <- epsilon
    theta.tmp <- theta
    R.tmp <- R
    if (!intercept)
      mu <- 0
    else
      mu <- log(sum(Y, na.rm = T) / sum(Omega * exp(param - mu.tmp)))
    if (!reff)
      alpha <- rep(0, n)
    if (!ceff)
      beta <- rep(0, p)
    grad_alpha <- grad(Y, cov, mu, alpha, beta, epsilon, theta)
    if (!reff)
      grad[1:n] <- 0
    if (!ceff)
      grad[(n + 1):(n + p)] <- 0
    step <- 1
    flag <- TRUE
    alpmat <- matrix(rep(alpha, p), nrow = n)
    betmat <- matrix(rep(beta, each = n), nrow = n)
    epsmat <- matrix(cov %*% epsilon, nrow = n)
    armijo.step <- 1e-3
    ref.obj <-
      sum(
        -Y * (mu + alpmat + betmat + epsmat + theta) + exp(mu + alpmat + betmat +
                                                             epsmat + theta),
        na.rm = T
      ) / m + lambda1 * R + sum(lambda2 * c(abs(epsilon), abs(alpha), abs(beta)))
    while (flag) {
      step <- 0.5 * step
      mat <-
        abs(c(alpha.tmp, beta.tmp, epsilon.tmp) - step * grad_alpha) - lambda2 * step
      maineff <-
        sign(c(alpha.tmp, beta.tmp, epsilon.tmp) - step * grad_alpha) * pmax(as.matrix(mat), 0)
      alpha <- maineff[1:n]
      beta <- maineff[(n + 1):(n + p)]
      epsilon <- maineff[(n + p + 1):(n + p + q)]
      alpmat <- matrix(rep(alpha, p), nrow = n)
      betmat <- matrix(rep(beta, each = n), nrow = n)
      epsmat <- matrix(cov %*% epsilon, nrow = n)
      param <- mu + alpmat + betmat + epsmat + theta
      diff <-
        sum(
          -Y * (mu + alpmat + betmat + epsmat + theta) + exp(mu + alpmat + betmat +
                                                               epsmat + theta),
          na.rm = T
        ) / m + lambda1 * R + sum(lambda2 * c(abs(epsilon), abs(alpha), abs(beta))) - ref.obj
      flag <- diff > thresh * abs(ref.obj)
      step <- 0.5 * step
    }
    param <- mu + alpmat + betmat + epsmat + theta
    grad_theta <- -Omega * (Y2 - exp(param))
    svd_theta <- rARPACK::svds(grad_theta,
                               k = 1,
                               nu = 1,
                               nv = 1)
    D_t <- -svd_theta$u %*% t(svd_theta$v)
    step <- 2
    flag <- TRUE
    obj0 <- low_bound
    obj0 <- obj0 + sum(-(Y * param) + exp(param), na.rm = T)

    obj0 <-
      obj0 + sum(lambda2 * c(abs(epsilon), abs(alpha), abs(beta))) + lambda1 *
      R
    while ((flag)) {
      step <- 0.5 * step
      if (lambda1 >= -sum(D_t * grad_theta)) {
        R_hat <- 0
        theta_hat <- matrix(rep(0, n * p), nrow = n)
        if (norm(theta_hat - theta.tmp, type = "F") ^ 2 == 0) {
          eta <- 1
        } else {
          eta <- min(1, step)
        }
        theta <- theta.tmp + eta * (theta_hat - theta.tmp)
        R <- R.tmp + eta * (R_hat - R.tmp)
      } else{
        R_hat <- U
        theta_hat <- U * D_t
        if (norm(theta_hat - theta.tmp, type = "F") ^ 2 == 0) {
          beta <- 1
          theta <- theta.tmp + eta * (theta_hat - theta.tmp)
          R <- R.tmp + eta * (R_hat - R.tmp)
        } else {
          eta <- min(1, step)
          theta <- theta.tmp + eta * (theta_hat - theta.tmp)
          R <- R.tmp + eta * (R_hat - R.tmp)
        }
      }
      param <- mu + alpmat + betmat + epsmat + theta
      obj <- low_bound
      obj <- obj + sum(-(Y * param) + exp(param), na.rm = T)

      obj <-
        obj + sum(lambda2 * c(abs(epsilon), abs(alpha), abs(beta))) + lambda1 *
        R
      flag <- (obj > obj0 + thresh * abs(obj0))
      if (step <= 1e-10)
        flagflag <- TRUE
    }
    #
    param <- mu + alpmat + betmat + epsmat + theta
    obj <- low_bound
    obj <- obj + sum(-(Y * param) + exp(param), na.rm = T)
    obj <-
      obj + sum(lambda2 * c(abs(epsilon), abs(alpha), abs(beta))) + lambda1 *
      R
    objective <- c(objective, obj)
    U <- obj / lambda1
    if (count == 1) {
      error <- 1
    } else{
      error = abs(objective[count] - objective[count - 1]) / abs(objective[count])
    }
    count = count + 1
    # if (trace.it && (count %% 10 == 0)) {
    #   print(paste("iter ", count, ": error ", error, " - objective: ", objective[count-1]))
    # }
  }
  param <- mu + alpmat + betmat + epsmat + theta

  return(
    list(
      X = param,
      mu = mu,
      alpha = alpha,
      beta = beta,
      epsilon = epsilon,
      theta = theta,
      objective = unlist(objective),
      iter = count,
      convergence = count < maxit
    )
  )
}
