##### HMM practicals
##### Session 4 - 12.05.2023
##### R code solution

### Exercise 1 - Likelihood evaluation, Poisson HMM

llk <- function(theta, x) {
  Gamma <- matrix(c(0.9, 0.1, 0.2, 0.8), byrow = TRUE, nrow = 2) # fixed tpm
  delta <- c(1, 0) # fixed initial distribution
  lambda <- theta[1:2] # vector containing parameters of the state-dependent distribution
  allprobs <- cbind(dpois(x, lambda[1]), dpois(x, lambda[2])) # densities of both state-dependent distributions for all states and time points (P matrix)
  foo <- delta %*% diag(allprobs[1, ]) # forward variable at time 1
  for (t in 2:length(x)) {
    foo <- foo %*% Gamma %*% diag(allprobs[t, ]) # calculate & overwrite forward variables up until alpha_t recusively
  }
  return(log(sum(foo))) # return log likelihood
}

## a)
# load and plot observations
load("simulatedPoisson.RData")
plot(x, type = "l") # plot time series
points(x, pch = 19)

## b)
# guessed values upon visual inspection:
theta <- c(1, 5)
llk(theta, x)

# try differing lambda1 first:
llk(c(1.1, 5), x) # worse than 1
llk(c(1.2, 5), x) # even worse
# larger than 1 seems to make it worse, let's try smaller
llk(c(0.9, 5), x) # bit better
llk(c(0.8, 5), x) # worse than 0.9
# now differ lambda2
llk(c(0.9, 4.9), x) # worse than 5
llk(c(0.9, 4.8), x) # even worse
# smaller than 5 seems to make it worse, let's try larger
llk(c(0.9, 5.1), x) # better
llk(c(0.9, 5.2), x) # a bit worse
theta <- c(0.9, 5.1) # let's stick with these values for now

## c)
# use code from exercise sheet
lambda1 <- seq(0, 2, length = 50)
lambda2 <- seq(3, 8, length = 50)
z <- matrix(NA, 50, 50)
for (i in 1:50) {
  for (j in 1:50) {
    z[i, j] <- llk(c(lambda1[i], lambda2[j]), x)
  }
}
contour(lambda1, lambda2, z, nlevels = 200)
abline(v = theta[1], lty = 2)
abline(h = theta[2], lty = 2)

# the parameters chosen in b) seem to maximize this likelihood
# so we'll continue with the same theta

## d)
# state probabilities for last observation:
# copy code from likelihood function used to calculate forward variables
Gamma <- matrix(c(0.9, 0.1, 0.2, 0.8), byrow = TRUE, nrow = 2) # fixed tpm
delta <- c(1, 0) # fixed initial distribution
lambda <- theta[1:2] # vector containing parameters of the state-dependent distribution
allprobs <- cbind(dpois(x, lambda[1]), dpois(x, lambda[2])) # densities of both state-dependent distributions for all states and time points (P matrix)
foo <- delta %*% diag(allprobs[1, ]) # forward variable at time 1
for (t in 2:length(x)) {
  foo <- foo %*% Gamma %*% diag(allprobs[t, ]) # calculate & overwrite forward variables up until alpha_t recusively
}
# divide each entry in alpha_30 by sum of both entries:
foo / sum(foo) # most likely, last observation is generated by second state

## e)
# plot data again
plot(x, type = "l") # plot time series
points(x, pch = 19)
points(15, x[15], pch = 19, col = "red")
# because of the high diagonal entries in the tpm (0.8 and 0.9),
# it is unlikely to switch states a lot, so staying in state 1 is more likely
# but also, a count of 3 is more likely under second Poisson distribution

# calculate state probabilities at time 15 given data x1, ..., x15
foo <- delta %*% diag(allprobs[1, ]) # forward variable at time 1
for (t in 2:15) {
  foo <- foo %*% Gamma %*% diag(allprobs[t, ]) # calculate & overwrite forward variables up until alpha_t recusively
}
foo / sum(foo)
# given x1-x15 and the parameter values, state 1 is 75% likely to have generated x15
# but: what happens after time 15 (smaller count values) also makes it more likely that
# the first state has generated x15 which is not considered here (will be covered later in the semestre)

## f)
# initialize matrix with probabilities of each observation belonging to either state 1 or state 2
stateProbs <- matrix(NA, nrow = length(x), ncol = 2)

# same code as in c), only that within the for-loop, 
# we now store the vector of probabilities for each observation belonging to state 1 or 2
foo <- delta %*% diag(allprobs[1, ])
stateProbs[1, ] <- foo / sum(foo)
for (t in 2:length(x)) {
  foo <- foo %*% Gamma %*% diag(allprobs[t, ])
  stateProbs[t, ] <- foo / sum(foo) # calculate & store probabilities
}

# column 1 contains the probabilities of belonging to state 1, 
# column 2 the probabilities of belonging to state 2
stateProbs
# thus for each row, the columns sum to 1
rowSums(stateProbs)

# plot state probabilities
plot(stateProbs[, 1], pch = 19, 
     main = "probability of observations to belong to state 1", 
     xlab = "observation", ylab = "probability")
points(stateProbs[, 2], pch = 19, col = "red")

# check which state has most likely generated each observation, 
# given all previous observations, to get most likely state sequence
states <- rep(1, length(x)) # initialize states vector with all states = 1
# all observations with a probability > 0.5 to belong to state 2 are assigned state 2 in our vector
states[which(stateProbs[, 2] > 0.5)] <- 2 
# plot time series, colour coded according to most probable state
plot(x, type = "l", main = "colour coded observations") 
points(x, col = states, pch = 19)


#### Exercise 2 - Gaussian HMM

## a)
# likelihood function for normal distributions
L <- function(theta, x) {
  Gamma <- diag(theta[1:2])
  Gamma[1, 2] <- 1 - Gamma[1, 1]
  Gamma[2, 1] <- 1 - Gamma[2, 2]
  delta <- solve(t(diag(2) - Gamma + 1), c(1, 1))
  mu <- theta[3:4]
  sigma <- theta[5:6]
  allprobs <- cbind(
    dnorm(x, mu[1], sigma[1]),
    dnorm(x, mu[2], sigma[2])
  )
  foo <- delta %*% diag(allprobs[1, ])
  for (t in 2:length(x)) {
    foo <- foo %*% Gamma %*% diag(allprobs[t, ])
  }
  return(log(sum(foo)))
}

## b)
# simulate HMM:
simNormHMM <- function(T, N, Gamma, delta, mu, sigma) {
  s <- numeric(T)
  s[1] <- sample(1:N, 1, prob = delta)
  for (t in 2:T) {
    s[t] <- sample(1:N, 1, prob = Gamma[s[t - 1], ])
  }
  x <- rnorm(T, mean = mu[s], sd = sigma[s])
  return(list(s = s, x = x))
}

N <- 2
Gamma <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = N, byrow = TRUE)
delta <- c(1, 0)
mu <- c(0.9, 6)
sigma <- c(1, 1)
data <- simNormHMM(30, N, Gamma, delta, mu, sigma)

L(c(0.9, 0.8, mu, sigma), data$x)

## c)
# use code from exercise sheet
mu1 <- seq(0, 2, length = 50)
mu2 <- seq(5, 7, length = 50)
z <- matrix(NA, 50, 50)
for (i in 1:50) {
  for (j in 1:50) {
    z[i, j] <- L(c(0.9, 0.8, mu1[i], mu2[j], sigma), data$x)
  }
}
contour(mu1, mu2, z, nlevels = 200)
abline(v = mu[1], lty = 2)
abline(h = mu[2], lty = 2)

# true parameter values might not always be the values maximizing the likelihood for each sample
