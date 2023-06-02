##### HMM practicals
##### Session 5 - 26.05.2023
##### R code solution

### Exercise 1

## a) overview dataset, some plots
# load data:
elephant <- read.table("elephant.txt")

# a movement track
plot(elephant$x, elephant$y, type = "l", xlab = "x", ylab = "y", main = "Track")
# a time series of step lengths
plot(elephant$step, type = "h", xlab = "time", 
     ylab = "step length in metres", main = "Time series")
# serial correlation present! probably two or more different states
# a histogram of step lengths
hist(elephant$step, breaks = 60, 
     xlab = "step length in metres", main = "Histogram", prob = T)

# for a 2-state HMM with state-dependent gamma distributions of the step lengths,
# starting values could be:
mu0 <- c(10, 250)
sigma0 <- c(10, 200)
# plot the gamma distributions with the chosen starting values:
x <- seq(0, 600, by = 0.1)
lines(x, dgamma(x,
  shape = mu0[1]^2 / sigma0[1]^2,
  scale = sigma0[1]^2 / mu0[1]
),
col = "orange", lwd = 3
)
lines(x, dgamma(x,
  shape = mu0[2]^2 / sigma0[2]^2,
  scale = sigma0[2]^2 / mu0[2]
),
col = "green4", lwd = 3
)


## b) likelihood function accounting for missing observations
# given likelihood function:
llk <- function(theta, x) {
  Gamma <- diag(theta[1:2]) # build Gamma matrix
  Gamma[1, 2] <- 1 - Gamma[1, 1]
  Gamma[2, 1] <- 1 - Gamma[2, 2]
  delta <- rep(1 / 2, 2) # delta is not estimated here
  mu <- theta[3:4]
  sigma <- theta[5:6]
  allprobs <- cbind( # get the density values of the gamma distributions
    dgamma(x, shape = mu[1]^2 / sigma[1]^2, scale = sigma[1]^2 / mu[1]),
    dgamma(x, shape = mu[2]^2 / sigma[2]^2, scale = sigma[2]^2 / mu[2])
  )
  foo <- delta %*% diag(allprobs[1, ]) # alpha_1
  for (t in 2:length(x)) {
    foo <- foo %*% Gamma %*% diag(allprobs[t, ]) # alpha_t
  }
  return(log(sum(foo))) # log-likelihood
}
# starting values:
theta <- c(0.8, 0.8, mu0, sigma0)
# calculate the likelihood for the first 52 steps
llk(theta, x = elephant$step[1:52])

# now for the first 100 observations:
llk(theta, x = elephant$step[1:100])
# NA because of the missing step lengths at
which(is.na(elephant$step))

# adapt likelihood function to account for missing values.
llk <- function(theta, x) {
  Gamma <- diag(theta[1:2])
  Gamma[1, 2] <- 1 - Gamma[1, 1]
  Gamma[2, 1] <- 1 - Gamma[2, 2]
  delta <- rep(1 / 2, 2)
  mu <- theta[3:4]
  sigma <- theta[5:6]
  allprobs <- matrix(1, length(x), 2) # allprobs = 1 for every obs
  ind <- which(!is.na(x))
  allprobs[ind, ] <- cbind( # now only put density values for non-missing obs
    dgamma(x[ind], shape = mu[1]^2 / sigma[1]^2, scale = sigma[1]^2 / mu[1]),
    dgamma(x[ind], shape = mu[2]^2 / sigma[2]^2, scale = sigma[2]^2 / mu[2])
  )
  foo <- delta %*% diag(allprobs[1, ])
  for (t in 2:length(x)) {
    foo <- foo %*% Gamma %*% diag(allprobs[t, ])
  }
  return(log(sum(foo)))
}

# try again for the first 100 values
llk(theta, x = elephant$step[1:100])


## c) likelihood function with scaling to address numerical underflow
# get likelihood of whole dataset and guesses parameter values
llk(theta, x = elephant$step)
# numerical underflow! We need to use scaling techniques.

llk <- function(theta, x) {
  Gamma <- diag(theta[1:2])
  Gamma[1, 2] <- 1 - Gamma[1, 1]
  Gamma[2, 1] <- 1 - Gamma[2, 2]
  delta <- rep(1 / 2, 2)
  mu <- theta[3:4]
  sigma <- theta[5:6]
  allprobs <- matrix(1, length(x), 2)
  ind <- which(!is.na(x))
  allprobs[ind, ] <- cbind(
    dgamma(x[ind], shape = mu[1]^2 / sigma[1]^2, scale = sigma[1]^2 / mu[1]),
    dgamma(x[ind], shape = mu[2]^2 / sigma[2]^2, scale = sigma[2]^2 / mu[2])
  )
  foo <- delta %*% diag(allprobs[1, ]) # alpha_1
  l <- log(sum(foo)) # this will add up to the log-likelihood
  phi <- foo / sum(foo) # standardised alpha_1: phi_1
  for (t in 2:length(x)) {
    foo <- phi %*% Gamma %*% diag(allprobs[t, ]) # phi_t-1 * Gamma * allprobs
    l <- l + log(sum(foo)) # updated log-likelihood
    phi <- foo / sum(foo) # phi_t
  }
  return(l)
}

# try again:
llk(theta, x = elephant$step) # works!

## d) likelihood function dealing with parameter constraints
# and returning the negative log-likelihood:
mllk <- function(theta.star, x) { # argument is now called theta.star
  theta <- c(plogis(theta.star[1:2]), exp(theta.star[3:6])) # retransform to natural parameters
  Gamma <- diag(theta[1:2])
  Gamma[1, 2] <- 1 - Gamma[1, 1]
  Gamma[2, 1] <- 1 - Gamma[2, 2]
  delta <- rep(1 / 2, 2)
  mu <- theta[3:4]
  sigma <- theta[5:6]
  allprobs <- matrix(1, length(x), 2)
  ind <- which(!is.na(x))
  allprobs[ind, ] <- cbind(
    dgamma(x[ind], shape = mu[1]^2 / sigma[1]^2, scale = sigma[1]^2 / mu[1]),
    dgamma(x[ind], shape = mu[2]^2 / sigma[2]^2, scale = sigma[2]^2 / mu[2])
  )
  foo <- delta %*% diag(allprobs[1, ]) 
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  for (t in 2:length(x)) {
    foo <- phi %*% Gamma %*% diag(allprobs[t, ]) 
    l <- l + log(sum(foo)) 
    phi <- foo / sum(foo) 
  }
  return(-l) # return minus log-likelihood (mllk) in order to maximize log-likelihood with the minimizer nlm()
}

# double check whether it returns the same log-likelihood value as before, but negative:
theta.star <- c(qlogis(theta[1:2]), log(theta[3:6]))
mllk(theta.star, x = elephant$step) # works!

## e) run optimizer with random starting values

# transform the starting values for theta.star
theta.star <- c(qlogis(theta[1:2]), log(theta[3:6]))
mod <- nlm(mllk, theta.star, x = elephant$step, print.level = 2)


llks <- rep(NA, 10)
mods <- list()
set.seed(2806)
for (k in 1:10) {
  theta.star <- c(qlogis(runif(2, 0, 1)), rep(log(runif(2, 1, 400)), 2))
  mods[[k]] <- nlm(mllk, theta.star, x = elephant$step)
  llks[k] <- -mods[[k]]$minimum
  print(k)
}

llks # all but two are the same. -> probably a global maximum

## f) backtransform parameters
theta.star <- mods[[which(llks == max(llks))]]$estimate # choose best result (highest likelihood)
mu <- exp(theta.star[3:4])
sigma <- exp(theta.star[5:6])
Gamma <- diag(plogis(theta.star[1:2]))
Gamma[1, 2] <- 1 - Gamma[1, 1]
Gamma[2, 1] <- 1 - Gamma[2, 2]

# get stationary distribution of the resulting MC
delta <- solve(t(diag(2) - Gamma + 1), c(1, 1))

## g) plot results
hist(elephant$step, breaks = 60, col = "grey", 
     xlab = "step length in metres", main = "Histogram", prob = T)

x <- seq(0, 600, by = 0.1)
lines(x, dgamma(x,
  shape = mu[1]^2 / sigma[1]^2,
  scale = sigma[1]^2 / mu[1]
) * delta[1],
col = "orange", lwd = 3
)
lines(x, dgamma(x,
  shape = mu[2]^2 / sigma[2]^2,
  scale = sigma[2]^2 / mu[2]
) * delta[2],
col = "green4", lwd = 3
)
lines(x, dgamma(x,
  shape = mu[1]^2 / sigma[1]^2,
  scale = sigma[1]^2 / mu[1]
) * delta[1] +
  dgamma(x,
    shape = mu[2]^2 / sigma[2]^2,
    scale = sigma[2]^2 / mu[2]
  ) * delta[2],
lty = 2, lwd = 3
)
legend("topright", c("state 1", "state 2", "mixture"),
  col = c("orange", "green4", "black"),
  lty = c(1, 1, 2), lwd = c(3, 3, 1)
)
