##### HMM practicals
##### Session 2 - 21.04.2023
##### R code solution

### Exercise 1 - independent mixture models

## a)
# read in the data set
data <- faithful

# overview data set:
?faithful
head(data)

## b)
# histogram
hist(data$eruptions,
  probability = TRUE, breaks = 20,
  main = "Eruption durations of the Old Faithful",
  xlab = "durations in minutes"
)

## c)
# function for calculating negative log-likelihood:
mllk <- function(theta, x) {
  mu <- theta[1:2]
  sigma <- theta[3:4]
  pi <- theta[5]
  logl <- sum(log(pi * dnorm(x, mu[1], sigma[1]) +
    (1 - pi) * dnorm(x, mu[2], sigma[2])))
  return(-logl)
}

# positive log-likelihood for given parameter set:
-mllk(c(2, 4.5, 1, 1, 0.5), data$eruptions)

## d)
# use optimiser:
mod <- nlminb(c(2, 4.5, 1, 1, 0.5), mllk,
  x = data$eruptions,
  lower = c(0, 0, 0, 0, 0), upper = c(Inf, Inf, Inf, Inf, 1)
)
# maximum likelihood estimates/ parameter values:
mod$par

## e)
# extract parameters
mu <- mod$par[1:2]
sigma <- mod$par[3:4]
pi <- c(mod$par[5], 1 - mod$par[5])

# to add the densities, you need to create a sequence for which you want to plot the normal densities
# (the x-coordinates for the plot)
z <- seq(1.5, 5.5, by = 0.01)

lines(z, pi[1] * dnorm(z, mean = mu[1], sd = sigma[1]),
  lwd = 3, col = "orange"
)
lines(z, pi[2] * dnorm(z, mean = mu[2], sd = sigma[2]),
  lwd = 3, col = "lightgreen"
)
# and joint density:
lines(z, pi[1] * dnorm(z, mean = mu[1], sd = sigma[1]) +
  pi[2] * dnorm(z, mean = mu[2], sd = sigma[2]),
lwd = 2, lty = 2
)
legend("top", c("density 1", "density 2", "mixture"),
  col = c("orange", "lightgreen", "black"),
  lty = c(1, 1, 2), lwd = c(3, 3, 2)
)

## f)
# simulate data from mixture model
# using the vectors mu and sigma from part d)

# sample size:
n <- 272

# generate sequence of which component is active
# get component 1 or 2 with probability pi1 and (1-pi1), respectively
s <- sample(1:2, n, replace = TRUE, prob = c(pi[1], pi[2]))

# generate realizations
# generate n eruption durations using a sequence of parameters of the corresponding mixture components
x <- rnorm(n, mean = mu[s], sd = sigma[s])

# # # or alternatively:
# # generate a number of how many observations should be simulated from component 1
# n1 <- rbinom(1, n , pi[1])
# # generate data from first mixture component
# x1 <- rnorm(n1, mean = mu[1], sd = sigma[1])
# # generate data from second mixture component
# x2 <- rnorm(n-n1, mean = mu[2], sd = sigma[2])
# x <- c(x1, x2) # combine data

## g)
# compare the simulated data to the real eruptions data
par(mfrow = c(2, 1))
hist(data$eruptions, probability = TRUE, breaks = 20, 
     main = "histogram of eruptions data", xlab = "eruptions", xlim = c(1, 5.5))
hist(x, probability = TRUE, breaks = 20, 
     main = "histogram of simulated data", xlim = c(1, 5.5))

# or compare it to multiple simulations of n=272:
par(mfrow = c(3, 3))
hist(data$eruptions, probability = TRUE, breaks = 20, col.main = "blue", 
     main = "histogram of eruptions data", xlab = "eruptions", xlim = c(1, 5.5))
for (i in 1:8){
  s <- sample(1:2, n, replace = TRUE, prob = c(pi[1], pi[2]))
  x <- rnorm(n, mean = mu[s], sd = sigma[s])
  hist(x, probability = TRUE, breaks = 20, 
       main = "histogram of simulated data", xlim = c(1, 5.5))
}

### Exercise 2 - Markov chains

## a)
# create the transition probability matrix:
N = 3
Gamma <- matrix(c(0.6, 0.3, 0.1,
                  0.2, 0.7, 0.1,
                  0.2, 0.3, 0.5), nrow = N, byrow = TRUE)

## b)
Gamma[1,] # probabilities of weather states tomorrow if today is sunny
# (60% sunny, 30% cloudy, 10% rainy)
Gamma %*% Gamma # Markov chain two steps from now
(Gamma %*% Gamma)[1,] # probabilities of weather states the day after tomorrow if today is sunny

# define a 7-step-tpm:
Gamma7 <- Gamma # initialize the 7-step tpm
for (i in 1:6) Gamma7 <- Gamma7 %*% Gamma
Gamma7

## c)
delta <- solve(t(diag(N) - Gamma + 1), rep(1, N))
delta # each row of the 7-step-tpm is already quite close to this distribution

## d) 
T = 1000 # length of Markov chain
s = numeric(T) # initialize vector
# state at t=1:
s[1] = 1 # sunny day

# states at t = 2,.., T
for(t in 2:T){
  s[t] = sample(x = 1:3, size = 1, prob = Gamma[s[t-1],])
}
# proportions of states in the simulated Markov chain:
table(s)/T
# similar to stationary distribution (even more so with larger T)
