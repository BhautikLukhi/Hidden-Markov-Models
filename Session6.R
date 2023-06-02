##### HMM practicals
##### Session 6 - 02.06.2023
##### R code solution

##### Exercise 1 #####
## a) data exploration

# daily share returns (percentage changes in the stock prices)
# for the Deutsche Bank share (2000-2016)
returns <- read.csv("http://www.rolandlangrock.com/Misc/deutschebank.csv")$x

plot(returns,
  type = "l", main = "share returns",
  ylab = "returns", xlab = "time"
)
summary(returns)
hist(returns, prob = T, breaks = 60)

# The histogram does not look like we have multiple states in the way we observed
# it for animal data.
# BUT: Thinking about the context, we can assume, that all returns
# fluctuate around the mean of zero.
# We assume that the standard deviations of the state-dependent distributions
# vary with respect to the underlying state.

# a 3 state Gaussian HMM could be a good fit:


## b) fit HMM, plot results
# likelihood code from slide 145
mllk <- function(theta.star, x, N) {
  mu <- theta.star[1:N]
  sigma <- exp(theta.star[N + 1:N])
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[(2 * N + 1):(length(theta.star))])
  Gamma <- Gamma / rowSums(Gamma)
  delta <- solve(t(diag(N) - Gamma + 1), rep(1, N))
  allprobs <- matrix(1, length(x), N)
  ind <- which(!is.na(x))
  for (j in 1:N) {
    allprobs[ind, j] <- dnorm(x[ind], mean = mu[j], sd = sigma[j])
  }
  foo <- delta %*% diag(allprobs[1, ])
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  for (t in 2:length(x)) {
    foo <- phi %*% Gamma %*% diag(allprobs[t, ])
    l <- l + log(sum(foo))
    phi <- foo / sum(foo)
  }
  return(-l)
}

N <- 3
# starting values for the nondiagonal of the tpm:
tpm.non.diag <- rep(0.05, N*(N-1))
# starting values for the means of the normal distributions all zero
mu0 <- rep(0, N)
# starting values for the variances of the normal distributions:
sigma0 <- c(0.01, 0.1, 1)
theta <- c(mu0, sigma0, tpm.non.diag)

# transform starting values 
# no transformation necessary for the mean
# log transform sigma to ensure positive values
theta.star <- c(theta[1:3], log(theta[4:6]), qlogis(theta[7:12]))

# run estimation / optimiser
mod <- nlm(mllk, theta.star, x = returns, N = N, print.level = 2, iterlim = 500)
theta.star.mle <- mod$estimate

# transform back the results (similar to first lines of code in mllk function):
mu <- theta.star.mle[1:N]
sigma <- exp(theta.star.mle[N + 1:N])
Gamma <- diag(N)
Gamma[!Gamma] <- exp(theta.star.mle[(2 * N + 1):(length(theta.star.mle))])
Gamma <- Gamma / rowSums(Gamma)
delta <- solve(t(diag(N) - Gamma + 1), rep(1, N))

# look at ml estimates
Gamma
options(scipen = 99) # to get small numbers printed as decimals
round(Gamma, 4) # rounding just for a quick look on the result
delta
mu
sigma

# plot the histogram & weighted state-dep. distributions
# plot the histogram
hist(returns,
  probability = TRUE, main = "histogram of log-returns",
  xlab = "time index", breaks = 60, col = "light grey"
)

# add the weighted densities of the state-dependent distributions
z <- seq(-0.2, 0.2, by = 0.001)
lines(z, delta[1] * dnorm(z, mu[1], sigma[1]), col = "red", lwd = 2) # for state 1
lines(z, delta[2] * dnorm(z, mu[2], sigma[2]), col = "blue", lwd = 2) # for state 2
lines(z, delta[3] * dnorm(z, mu[3], sigma[3]), col = "green", lwd = 2) # for state 3

# add the overall density (sum of the state-dependent densities)
lines(z, delta[1] * dnorm(z, mu[1], sigma[1])
  + delta[2] * dnorm(z, mu[2], sigma[2])
  + delta[3] * dnorm(z, mu[3], sigma[3]), lwd = 2, lty = 2)


legend("topright",
  legend = c("state 1", "state 2", "state 3", "overall"),
  lwd = 2, col = c("red", "blue", "green", "black"), lty = c(1,1,1,2)
)
