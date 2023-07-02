######## Likelihood function for share returns data ########


### multivariate HMM with contemporaneous conditional and longitudinal conditional independence assumption ###
# both likelihood functions
library(mvtnorm)

mllk1 <- function(theta.star, x) {
  Gamma <- diag(plogis(theta.star[1:2]))
  Gamma[1, 2] <- 1 - Gamma[1, 1]
  Gamma[2, 1] <- 1 - Gamma[2, 2]
  delta <- solve(t(diag(2) - Gamma + 1), c(1, 1))
  sigma.CI <- exp(theta.star[3:4])
  sigma.MS <- exp(theta.star[5:6])
  rho <- (exp(theta.star[7:8]) - 1) / (exp(theta.star[7:8]) + 1)
  allprobs <- matrix(1, dim(x)[1], 2)
  for (j in 1:2) {
    sigma <- matrix(0, 2, 2)
    sigma[1, 1] <- sigma.CI[j]^2
    sigma[1, 2] <- sigma[2, 1] <- sigma.CI[j] * sigma.MS[j] * rho[j]
    sigma[2, 2] <- sigma.MS[j]^2
    allprobs[, j] <- dmvnorm(x, mean = c(0, 0), sigma)
  }
  foo <- delta %*% diag(allprobs[1, ])
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  for (t in 2:dim(x)[1]) {
    foo <- phi %*% Gamma %*% diag(allprobs[t, ])
    l <- l + log(sum(foo))
    phi <- foo / sum(foo)
  }
  return(-l)
}

mllk2 <- function(theta.star, x) {
  Gamma <- diag(plogis(theta.star[1:2]))
  Gamma[1, 2] <- 1 - Gamma[1, 1]
  Gamma[2, 1] <- 1 - Gamma[2, 2]
  delta <- solve(t(diag(2) - Gamma + 1), c(1, 1))
  sigma.CI <- exp(theta.star[3:4])
  sigma.MS <- exp(theta.star[5:6])
  allprobs <- matrix(1, dim(x)[1], 2)
  for (j in 1:2) {
    allprobs[, j] <- dnorm(x$CI, mean = 0, sd = sigma.CI[j]) * dnorm(x$MS, mean = 0, sd = sigma.MS[j])
  }
  foo <- delta %*% diag(allprobs[1, ])
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  for (t in 2:dim(x)[1]) {
    foo <- phi %*% Gamma %*% diag(allprobs[t, ])
    l <- l + log(sum(foo))
    phi <- foo / sum(foo)
  }
  return(-l)
}




