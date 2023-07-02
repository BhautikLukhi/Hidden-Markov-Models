##### HMM practicals
##### Session 8 - 23.06.2023
##### R code solution

### Exercise 1

## load the models from Practical Session 7:
# install.packages("moveHMM")
library(moveHMM)
# load data:
rawdata <- read.table("elephant_rawdata.txt", header = T)
data <- prepData(rawdata, type = "UTM")
N <- 3
stepMean <- seq(10, 300, length = N) # mean of gamma distribution (step lengths)
stepSD <- seq(10, 100, length = N) # SD of gamma distribution (step lengths)
stepPar <- c(stepMean, stepSD)
mod3 <- fitHMM(data,
  nbStates = N, stepPar0 = stepPar, verbose = 1,
  angleDist = "none", stationary = T
)

# define some colors that I will use 
colors <- c("#E69F00", "#56B4E9", "#009E73")

## a) model checking: marginal distribution
#
# 3 states
hist(data$step, prob = T, breaks = 60)
# extract parameters from fitted model
delta <- mod3$mle$delta
mean <- mod3$mle$stepPar[1, ]
sd <- mod3$mle$stepPar[2, ]
# calculate shape and scale
shape <- mean^2 / sd^2
scale <- sd^2 / mean
lines(dgamma(1:800, shape = shape[1], scale = scale[1]) * delta[1], col = colors[1], lwd = 2)
lines(dgamma(1:800, shape = shape[2], scale = scale[2]) * delta[2], col = colors[2], lwd = 2)
lines(dgamma(1:800, shape = shape[3], scale = scale[3]) * delta[3], col = colors[3], lwd = 2)
lines(dgamma(1:800, shape = shape[1], scale = scale[1]) * delta[1] +
  dgamma(1:800, shape = shape[2], scale = scale[2]) * delta[2] +
  dgamma(1:800, shape = shape[3], scale = scale[3]) * delta[3], col = "black", lwd = 2)

## b) simulation based checks
# copy and adapt simulation function from Practical Session 3
simGammaHMM <- function(T, N, Gamma, delta, shape, scale) {
  s <- numeric(T)
  s[1] <- sample(1:N, 1, prob = delta)
  for (t in 2:T) {
    s[t] <- sample(1:N, 1, prob = Gamma[s[t - 1], ])
  }
  x <- rgamma(T, shape = shape[s], scale = scale[s])
  return(list(s = s, x = x))
}
# extract Gamma matrix from fitted model
Gamma <- mod3$mle$gamma
T <- dim(data)[1]
hmm_simulated <- simGammaHMM(T, N, Gamma, delta, shape, scale)

# compare histograms
par(mfrow = c(2, 1))
hist(data$step, main = "data", xlim = c(0, 1000))
hist(hmm_simulated$x, main = "simulation", xlim = c(0, 1000))
# compare time series:
plot(data$step[1:300], main = "data", ylim = c(0, 800), type = "l")
plot(hmm_simulated$x[1:300], main = "simulation", ylim = c(0, 800), type = "l")
# compare autocorrelation
acf(data$step, na.action = na.pass, ylim = c(0, 1), main = "data")
acf(hmm_simulated$x, na.action = na.pass, ylim = c(0, 1), main = "simulation")
# compare summary statistics
summary(data$step)
summary(hmm_simulated$x)

### Exercise 2 - State decoding

## b)
plotStates(mod3)

## c)
sp3 <- stateProbs(mod3)
round(sp3[1:100, 1:3], 3)

## d)
vitstates3 <- moveHMM::viterbi(mod3)
vitstates3[1:100]

## e)
plot(data$step, type = "h", xlab = "time", ylab = "step length in metres", 
     main = "Time series without state decoding")
plot(data$step, type = "h", xlab = "time", ylab = "step length in metres", 
     main = "Time series with 3 decoded states", col = colors[vitstates3])


## f)
# define state sequence based on local decoding
spStates <- apply(sp3, 1, which.max)
# index of the observations for which the local decoding differs from the global decoding
which(spStates != vitstates3)


par(mfrow = c(1, 1))
plot(vitstates3, pch = 16, col = colors[spStates], 
     ylab = "states according to viterbi", main = "Decoded time series")
legend("topright", inset = c(0.1, 0.1), title = "states acc. to state prob.s", 
       legend = c("State 1", "State 2", "State 3"), pch = 16, col = colors)


par(mfrow = c(2, 1))
plot(data$step[1200:1300], type = "h", 
     xlab = "time", ylab = "step length in metres", 
     main = "Global decoding", 
     col = colors[vitstates3[1200:1300]], lwd = 2)
plot(data$step[1200:1300], type = "h", 
     xlab = "time", ylab = "step length in metres", 
     main = "Local decoding", 
     col = colors[spStates[1200:1300]], lwd = 2)


### Bonus exercise
# Viterbi code for general N-state case

viterbi <- function(x, mu, sigma, gamma, delta, N) {
  n <- length(x)
  allprobs <- matrix(1, n, N)
  ind <- which(!is.na(x))
  for (j in 1:N) {
    allprobs[ind, j] <- dgamma(x[ind], shape = mu[j]^2 / sigma[j]^2, scale = sigma[j]^2 / mu[j])
  }
  xi <- matrix(0, n, N)
  foo <- delta * allprobs[1, ]
  xi[1, ] <- foo / sum(foo)
  for (t in 2:n) {
    foo <- apply(xi[t - 1, ] * gamma, 2, max) * allprobs[t, ]
    xi[t, ] <- foo / sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ])
  for (t in (n - 1):1) {
    iv[t] <- which.max(gamma[, iv[t + 1]] * xi[t, ])
  }
  iv
}

decStates <- viterbi(data$step, mean, sd, Gamma, delta, N = 3)
plot(decStates, col = colors[decStates], pch = 16)
