##### HMM practicals
##### Session 1 - 14.04.2023
##### R code solution

### Exercise 1
# set workspace:
setwd("P:/HMM_folder") # change according to your folder

## a)
#read in dataset
data <- read.table("muskox.txt", header = TRUE)

## b)
# get an overview
head(data)
summary(data)

# only one muskox in the dataset
unique(data$ID)

# dimension of dataset: 1440 observations, 14 variables
dim(data)

# calculation of mean and standard deviation for step length
mean(data$step) # doesn't work because of the missing values
# set na.rm=TRUE to ignore the NAs for the calculation
mean(data$step, na.rm = T)
sd(data$step, na.rm = T)

# how many NAs in step length?
ind <- which(is.na(data$step))
length(ind)
# or
sum(is.na(data$step))

## c)
# euclidean distance (step length) between the first two observations
sqrt((data$xutm[2] - data$xutm[1])^2 + (data$yutm[2] - data$yutm[1])^2)
# compare calculated step length to the one in the data set
data$step[1]

## d)
# plot step length time series
plot(data$step, type = "h", 
     xlab = "time index", ylab = "step length", 
     main = "step length of muskox")
# autocorrelation present in the data:
acf(data$step, na.action = na.pass)

## e)
# histograms for step and angle
hist(data$step, probability = TRUE, breaks = 20,
     main = "histogram for step length", 
     xlab = "step length", ylab = "density")
hist(data$angle, probability = TRUE, breaks = 20,
     main = "histogram for turning angle", 
     xlab = "turning angle", ylab = "density")

## f)
# plot the track of the animal
plot(data$x, data$y, type = "l", 
     main = "movement of muskox M07", xlab = "x", ylab = "y")

# add green points to the plot to mark the densely vegetated areas
ind_green <- which(data$habitat == "Dense_vegetated")
points(data$x[ind_green], data$y[ind_green], col = "green", pch = 20)

# you could continue with the other habitats:
ind_orange <- which(data$habitat == "Sparse_vegetated")
points(data$x[ind_orange], data$y[ind_orange], col = "darkorange", pch = 20)

ind_brown <- which(data$habitat == "Bareground")
points(data$x[ind_brown], data$y[ind_brown], col = "brown", pch = 20)

ind_ice <- which(data$habitat == "Water_snow_ice")
points(data$x[ind_ice], data$y[ind_ice], col = "blue", pch = 20)

## g*)
library(tidyverse)
strptime(data$time, "%H:%M:%S") %>% hour() %>% diff() %>% table()

### Exercise 2
## a) (functions)
mysum <- function(a, b) {
  c <- a + b
  return(c)
}
mysum(3, 5)

# write function for step length:
calc.step <- function(x1, x2, y1, y2) {
  d <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  return(d)
}
# e.g.:
calc.step(x1 = 2, y1 = 10, x2 = 2, y2 = 15)
# test and compare for first step length:
calc.step(x1 = data$xutm[1], y1 = data$yutm[1], x2 = data$xutm[2], y2 = data$yutm[2])
data$step[1]

## b) (for-loops)
# e.g. summing all values between 1 and 100:
result <- 0
for (i in 1:100) {
  result <- result + i
}
result
# e.g. filling a vector with the squares of all values between 1 to 10:
v <- numeric(10)
for (i in 1:10) {
  v[i] <- i^2
}
v

# calculate first 100 step lengths using function calc.step
step.length <- numeric(100)
for (i in 1:100) {
  step.length[i] <- calc.step(x1 = data$xutm[i], x2 = data$xutm[i + 1], y1 = data$yutm[i], y2 = data$yutm[i + 1])
}
# compare with variable:
step.length
cbind(step.length, data$step[1:100])


## c) (the normal distribution)

# plot the density of the normal distribution
x <- seq(-4, 4, 0.01)
plot(x, dnorm(x), type = "l")

dnorm(0.5) # value of the density function at x=0.5
points(x = 0.5, y = dnorm(0.5), col = "red", pch = 20, cex = 1.5)

pnorm(0.5) # proportion of distribution left of x=0.5
polygon(c(x[x <= 0.5], 0.5),
  c(dnorm(x)[x <= 0.5], 0),
  col = "#FF00004C"
)

qnorm(0.5) # x value for which 0.5(half) of the distribution is left of it: quantile
points(x = qnorm(0.5), y = dnorm(qnorm(0.5)), col = "blue", pch = 20, cex = 1.5)
polygon(c(x[x <= qnorm(0.5)], qnorm(0.5)),
  c(dnorm(x)[x <= qnorm(0.5)], 0),
  col = "#0000FF4C"
)

rnorm(1) # 1 random value from the standard normal distribution
set.seed(123)
rnorm(1) # another random value from the standard normal distribution
set.seed(123)
rnorm(1) # the same value as before because of the seed

