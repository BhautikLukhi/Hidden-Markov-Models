##### HMM practicals
##### Session 3 - 05.05.2023
##### R code solution

### Exercise 1 - Markov chains

## a)
# function for simulating a Markov chain
simMC <- function(T, N, Gamma, delta) {
  s <- numeric(T)
  s[1] <- sample(1:N, 1, prob = delta)
  for (t in 2:T) {
    s[t] <- sample(1:N, 1, prob = Gamma[s[t - 1], ])
  }
  return(s)
}

## b)
# T = 300, N = 4, Gamma and delta of own choice:
Gamma <- matrix(0.05, nrow = 4, ncol = 4)
diag(Gamma) <- 0.85
delta <- c(1,0,0,0)
s <- simMC(300, 4, Gamma, delta)

## c)
# function to calculate stationary distribution 
getDelta <- function(Gamma, N) {
  delta <- solve(t(diag(N) - Gamma + 1), rep(1, N)) 
  return(delta)
}

# or you could do it without having to define N as input
getDelta2 <- function(Gamma) {
  N <- nrow(Gamma)
  delta <- solve(t(diag(N) - Gamma + 1), rep(1, N))
  return(delta)
}

## d)
# get stationary for given tpm
Gamma <- matrix(0.05, nrow = 3, ncol = 3)
diag(Gamma) <- 0.9
getDelta(Gamma, 3)

#### Exercise 2 - HMMs

## a)
# simulate 3-state-HMM
# parameters:
T = 300
N = 3
delta <- getDelta(Gamma, N)
mu <- c(-2, 0, 2)
sigma <- rep(0.7, 3)

set.seed(123)
# simulate 3-state Markov chain
s <- simMC(T, N, Gamma, delta)
# then simulate observations
x <- rnorm(T, mean = mu[s], sd = sigma[s])

# alternatively, simulate both in the same for-loop
# initialize vectors: 
s <- numeric(T)
x <- numeric(T)
set.seed(123)
s[1] <- sample(1:3, 1, prob = delta)
x[1] <- rnorm(1, mean = mu[s[1]], sd = sigma[s[1]])
for (t in 2:T) {
  # generate the state sequence 
  s[t] <- sample(1:3, 1, prob = Gamma[s[t - 1], ])
  # use generated state sequence to generate the state-dependent observations:
  x[t] <- rnorm(1, mean = mu[s[t]], sd = sigma[s[t]])
}


## b)
# plot the observation sequence from a)
plot(x, main = "simulated data", xlab = "time t", pch = 19)
# plot observations in different colours depending on their underlying state
plot(x, main = "simulated data", xlab = "time t", pch = 19, col = s)

# we can also plot the observations as a line
plot(x, main = "simulated data", xlab = "time t", type = "l")
# and then add points for each state in a different colour
points(x, col = s, pch = 19)
# or choose the colours ourselves:
points(which(s == 1), x[s == 1], col = "orange", pch = 19)
points(which(s == 2), x[s == 2], col = "blue", pch = 19)
points(which(s == 3), x[s == 3], col = "red", pch = 19)


## d*)
# function to simulate HMM with Gaussian distributions:
simNormHMM <- function(T, N, Gamma, delta, mu, sigma) {
  s <- simMC(T, N, Gamma, delta)
  x <- rnorm(T, mean = mu[s], sd = sigma[s])
  return(list(s = s, x = x))
}

## c)
# things to try:
# smaller diagonal entries in t.p.m.
# smaller or larger difference between state-dependent means
# smaller or larger state-dependent variances

## d*)
# 4-state HMM
N <- 4
Gamma <- matrix(0.05, 4, 4)
diag(Gamma) <- 0.85
Gamma
delta <- c(0.2, 0.3, 0.25, 0.25)
mu <- c(-2, 0, 2, 4)
sigma <- rep(0.5, 4)

set.seed(123)
simHMM4 <- simNormHMM(T, N, Gamma, delta, mu, sigma)
plot(simHMM4$x, main = "simulated data", xlab = "time t", type = "l")
points(simHMM4$x, pch = 19, col = simHMM4$s)
