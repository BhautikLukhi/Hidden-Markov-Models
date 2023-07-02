##### HMM practicals
##### Session 9 - 30.06.2023
##### R code solution

### Exercise 1

## a)

data <- read.table("elephant_full.txt", header = T)
head(data)

# exclude missing values
data.new <- data[2:1521, ]

# get an overview of the dataset
par(mfrow = c(2, 1), mar = c(1, 4, 1, 1))
plot(data.new$step[1:750], type = "h", xlab = "", ylab = "step length in metres", main = "", xaxt = "n")
plot(data.new$angle[1:750], type = "h", xlab = "", ylab = "turning angle in radians", main = "", xaxt = "n")
par(mar = c(5.1, 4.1, 4.1, 2.1))
hist(data.new$step, breaks = 60, col = "grey", xlab = "step length in metres", main = "Histogram", prob = T)
hist(data.new$angle, breaks = 60, col = "grey", xlab = "step length in metres", main = "Histogram", prob = T)
summary(data.new$step)
summary(data.new$angle)


## b)

# function from lecture slide 140: (for univariate/ step length HMM with N states)
mllk <- function(theta.star, x, N) {
  mu <- exp(theta.star[1:N])
  sigma <- exp(theta.star[N + 1:N])
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[(2 * N + 1):(length(theta.star))])
  Gamma <- Gamma / rowSums(Gamma)
  delta <- solve(t(diag(N) - Gamma + 1), rep(1, N))
  allprobs <- matrix(1, length(x), N)
  ind <- which(!is.na(x))
  for (j in 1:N) {
    allprobs[ind, j] <- dgamma(x[ind], shape = mu[j]^2 / sigma[j]^2, 
                               scale = sigma[j]^2 / mu[j])
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

# rewrite it to include turning angles as well:
library(CircStats)

L <- function(theta.star, x, N) {
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N - 1) * N)])
  Gamma <- Gamma / rowSums(Gamma)
  delta <- solve(t(diag(N) - Gamma + 1), rep(1, N))
  mu.step <- exp(theta.star[(N - 1) * N + 1:N])
  sigma <- exp(theta.star[(N - 1) * N + (N + 1):(2 * N)])
  kappa <- exp(theta.star[(N - 1) * N + 2 * N + 1:N])
  mu.turn <- rep(0, N)
  allprobs <- matrix(1, dim(x)[1], N)
  for (j in 1:N) {
    allprobs[, j] <- dgamma(x$step, shape = mu.step[j]^2 / sigma[j]^2, 
                            scale = sigma[j]^2 / mu.step[j]) *
      dvm(x$angle, mu.turn[j], kappa[j])
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


N <- 3
mu.step0 <- seq(10, 300, length = N)
sigma0 <- seq(10, 200, length = N)
kappa0 <- seq(0.1, 8, length = N)
theta.star <- c(rep(-2, (N - 1) * N), log(mu.step0), log(sigma0), log(kappa0))
mod <- nlm(L, theta.star, x = data.new, N = N, print.level = 1)

## c)

# transform estimates to natural parameters
Gamma <- diag(N)
Gamma[!Gamma] <- exp(mod$estimate[1:((N - 1) * N)])
Gamma <- Gamma / rowSums(Gamma)
delta <- solve(t(diag(N) - Gamma + 1), rep(1, N))
mu.step <- exp(mod$estimate[(N - 1) * N + 1:N])
sigma <- exp(mod$estimate[(N - 1) * N + (N + 1):(2 * N)])
kappa <- exp(mod$estimate[(N - 1) * N + 2 * N + 1:N])


## d)

# plot state-dependent distributions
par(mfrow = c(2, 1))

hist(data.new$step, probability = TRUE, breaks = 60, col = "light grey", xlab = "step length in metres", main = "State-dependent distributions")
z <- seq(0, 1200, by = 0.1)
lines(z, delta[1] * dgamma(z, shape = mu.step[1]^2 / sigma[1]^2, scale = sigma[1]^2 / mu.step[1]), col = "red", lwd = 2)
lines(z, delta[2] * dgamma(z, shape = mu.step[2]^2 / sigma[2]^2, scale = sigma[2]^2 / mu.step[2]), col = "blue", lwd = 2)
lines(z, delta[3] * dgamma(z, shape = mu.step[3]^2 / sigma[3]^2, scale = sigma[3]^2 / mu.step[3]), col = "green", lwd = 2)

hist(data.new$angle, probability = TRUE, breaks = 40, col = "light grey", xlab = "turning angle in radians", main = "")
z <- seq(-pi, pi, by = 0.01)
lines(z, delta[1] * dvm(z, 0, kappa[1]), col = "red", lwd = 2)
lines(z, delta[2] * dvm(z, 0, kappa[2]), col = "blue", lwd = 2)
lines(z, delta[3] * dvm(z, 0, kappa[3]), col = "green", lwd = 2)




##### Exercise 2 #####

## a) read and plot data

data <- read.table("share_returns.txt", header = T)
head(data)

# get an overview of the dataset
par(mfrow = c(2, 1), mar = c(1, 4, 1, 1))
plot(data$CI, type = "h", ylab = "Citygroup returns", xlab = "", xaxt = "n", main = "", ylim = c(-0.4, 0.4))
plot(data$MS, type = "h", ylab = "Morgan Stanley returns", xlab = "", xaxt = "n", main = "", ylim = c(-0.4, 0.4))


## c) fit model with contemporaneous conditional independence assumption
source("mllk_shareReturns.R")

sigma.MS <- sigma.CI <- c(0.01, 0.04)
theta.star <- c(qlogis(c(0.99, 0.99)), log(sigma.CI), log(sigma.MS))

mod.cci <- nlm(mllk2, theta.star, x = data, print.level = 1, stepmax = 10)

# natural parameters
Gamma.cci <- diag(plogis(mod.cci$estimate[1:2]))
Gamma.cci[1, 2] <- 1 - Gamma.cci[1, 1]
Gamma.cci[2, 1] <- 1 - Gamma.cci[2, 2]
delta.cci <- solve(t(diag(2) - Gamma.cci + 1), c(1, 1))
sigma.CI.cci <- exp(mod.cci$estimate[3:4])
sigma.MS.cci <- exp(mod.cci$estimate[5:6])


## d) fit model with longitudinal conditional independence assumption

sigma.MS <- sigma.CI <- c(0.01, 0.04)
rho <- c(0.6, 0.6)
trho <- log((1 + rho) / (1 - rho))
theta.star <- c(qlogis(c(0.99, 0.99)), log(sigma.CI), log(sigma.MS), trho)

mod.lci <- nlm(mllk1, theta.star, x = data, print.level = 1, stepmax = 10)

# natural parameters
Gamma.lci <- diag(plogis(mod.lci$estimate[1:2]))
Gamma.lci[1, 2] <- 1 - Gamma.lci[1, 1]
Gamma.lci[2, 1] <- 1 - Gamma.lci[2, 2]
delta.lci <- solve(t(diag(2) - Gamma.lci + 1), c(1, 1))
sigma.CI.lci <- exp(mod.lci$estimate[3:4])
sigma.MS.lci <- exp(mod.lci$estimate[5:6])
rho.lci <- (exp(mod.lci$estimate[7:8]) - 1) / (exp(mod.lci$estimate[7:8]) + 1)

# compare estimated parameters from both models
Gamma.cci
Gamma.lci
delta.cci
delta.lci
sigma.CI.cci
sigma.CI.lci
sigma.MS.cci
sigma.MS.lci


## e)

# number of parameters for both models: 2 sigmas for CI, 2 sigmas for MS, 2 tpm entries
# for longitudinal conditional independence, additionally 2 rhos for covariance between CI an MS
# AIC
aic.cci <- 2 * mod.cci$minimum + 2 * 6
aic.cci
aic.lci <- 2 * mod.lci$minimum + 2 * 8
aic.lci

# BIC
T <- dim(data)[1]
bic.cci <- 2 * mod.cci$minimum + log(T) * 6
bic.cci
bic.lci <- 2 * mod.lci$minimum + log(T) * 8
bic.lci
