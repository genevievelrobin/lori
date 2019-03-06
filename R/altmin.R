altmin <-
  # internal function to do alternating minimization
  function(Y,
           lambda1,
           lambda2,
           cov = NULL,
           rank.max = 10,
           thresh = 1e-6,
           maxit = 1e3,
           trace.it = T,
           intercept = F,
           reff = T,
           ceff = T)
  {
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
    if (!is.null(cov)){
      epsilon <- rep(0, ncol(cov))
    } else{
      epsilon <- 0
    }
    alpmat <- matrix(rep(alpha, p), nrow = n)
    betmat <- matrix(rep(beta, each = n), nrow = n)
    epsmat <- matrix(cov %*% epsilon, nrow = n)
    error = 1
    objective <- NULL
    count = 1
    Y2 <- Y
    Y2[is.na(Y)] <- 0
    X <- mu + alpmat + betmat + epsmat + theta
    Ytalp <- c(Y[!is.na(Y)])
    while (error > thresh && count < maxit) {
      mu.tmp <- mu
      alpha.tmp <- alpha
      beta.tmp <- beta
      epsilon.tmp <- epsilon
      d.tmp <- svd(theta)$d
      if(!intercept) mu <- 0 else mu <- log(sum(Y, na.rm=T)/sum(Omega*exp(X-mu.tmp)))
      if(!reff) alpha <- rep(0,n)
      if(!ceff) beta <- rep(0,p)
      grad <- grad(Y, cov, mu, alpha, beta, epsilon, theta)
      if(!reff) grad[1:n] <- 0
      if(!ceff) grad[(n+1):(n+p)] <- 0
      step <- 1
      flag <- TRUE
      alpmat <- matrix(rep(alpha, p), nrow=n)
      betamat <- matrix(rep(beta, each=n), nrow=n)
      epsilonmat <- matrix(cov%*% epsilon, nrow=n)
      armijo.step <- 1e-3
      ref.obj <-
        sum(-Y * (mu+alpmat + betmat + epsmat + theta) + exp(mu+alpmat + betmat +
                                                             epsmat + theta),
            na.rm = T) / m + lambda1 * sum(d.tmp) + sum(lambda2 * c(abs(epsilon),abs(alpha), abs(beta)))
      while(flag){
        #print(step)
        maineff <- sign(c(alpha.tmp, beta.tmp, epsilon.tmp)-step*grad)*pmax(abs(c(alpha.tmp, beta.tmp, epsilon.tmp)-step*grad)-step*lambda2, 0)
        alpha <- maineff[1:n]
        beta <- maineff[(n+1):(n+p)]
        epsilon <- maineff[(n+p+1):(n+p+q)]
        alpmat <- matrix(rep(alpha, p), nrow=n)
        betmat <- matrix(rep(beta, each=n), nrow=n)
        epsmat <- matrix(cov%*% epsilon, nrow=n)
        diff <-
          sum(-Y * (mu + alpmat + betmat + epsmat + theta) + exp(mu + alpmat + betmat +
                                                                 epsmat + theta),
              na.rm = T) / m + lambda1 * sum(d.tmp) + sum(lambda2 * c(abs(epsilon),abs(alpha), abs(beta))) - ref.obj
        flag <- diff > thresh * abs(ref.obj)
        step <- 0.5*step
      }
      grad_theta <- (-Y2 + exp(mu + alpmat + betmat + epsmat + theta)) / m
      grad_theta <- sweep(grad_theta, 2, colMeans(grad_theta))
      grad_theta <- sweep(grad_theta, 1, rowMeans(grad_theta))
      flag <- TRUE
      step <- 1
      armijo.step <- 1e-2
      ref.obj <-
        sum(-Y * (mu + alpmat + betmat + epsmat + theta) + exp(mu + alpmat + betmat +
                                                               epsmat + theta),
            na.rm = T) / m + lambda1 * sum(d.tmp) + sum(lambda2 * c(abs(epsilon),abs(alpha), abs(beta)))
      while (flag) {
        #print(step)
        svd_theta <-
          svd(theta - step * grad_theta)
        d <- svd_theta$d[1:rank.max]
        u <- svd_theta$u[, 1:rank.max]
        v <- svd_theta$v[, 1:rank.max]
        d <- d[d > (lambda1 * step)] - lambda1 * step
        if (length(d) == 0) {
          theta2 <- matrix(rep(0, n * p), nrow = n)
        } else if (length(d) == 1) {
          u <- svd_theta$u[, 1:length(d)]
          v <- svd_theta$v[, 1:length(d)]
          theta2 <- d * u %*% t(v)
        } else {
          u <- svd_theta$u[, 1:length(d)]
          v <- svd_theta$v[, 1:length(d)]
          theta2 <- u %*% diag(d) %*% t(v)
        }
        diff <-
          sum(-Y * (mu + alpmat + betmat + epsmat + theta2) + exp(mu + alpmat + betmat +
                                                                  epsmat + theta2),
              na.rm = T) / m + lambda1 * sum(d) + sum(lambda2 * c(abs(epsilon),abs(alpha), abs(beta)))- ref.obj
        flag <- diff > thresh * abs(ref.obj)
        step <- 0.5 * step
      }
      theta <- theta2
      X <- mu + alpmat + betmat + epsmat + theta
      objective = c(objective,-sum(Y * X - Omega * exp(X), na.rm = T) / m + lambda1 *
                      sum(d) + sum(lambda2 * c(abs(epsilon),abs(alpha), abs(beta))))
      if(count==1){
        error <- 1
      } else{
        error = abs(objective[count] - objective[count - 1]) / abs(objective[count])
      }
      count = count + 1
      epsilon <- as.vector(epsilon)
      names(epsilon) <- colnames(cov)
      rownames(X) <- rownames(Y)
      colnames(X) <- colnames(Y)
      rownames(theta) <- rownames(Y)
      colnames(theta) <- colnames(Y)
    }
    return(structure(
      list(
        X = X,
        mu = mu,
        alpha = alpha,
        beta = beta,
        epsilon = epsilon,
        theta = theta,
        objective = unlist(objective),
        iter = count,
        convergence = count < maxit
      )
    ))
  }
