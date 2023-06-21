##### HMM practicals
##### Session 7 - 16.06.2023
##### R code solution

### Exercise 1

## a) data preparation and exploration
# load package (install beforehand if you haven't already)
# install.packages("moveHMM")
library(moveHMM)

# load data:
rawdata <- read.table("elephant_rawdata.txt", header = T)
head(rawdata)

data <- prepData(rawdata, type = "UTM")
head(data)

plot(data)

## b) fit models with 2,3,4,5 states
# create empty list
mods <- list()
 for (N in 2:5){
  stepMean <- seq(10, 300, length = N) # mean of gamma distribution (step lengths)
  stepSD <- seq(10, 100, length = N) # SD of gamma distribution (step lengths)
  stepPar <- c(stepMean, stepSD)
  mods[[N]] <- fitHMM(data, nbStates = N, stepPar0 = stepPar, verbose = 1, 
                      angleDist = "none", stationary = T)

 }

## c) 10 sets of random starting values each

mods2 <- mods3 <- mods4 <- mods5 <-list()
mllks2 <- mllks3 <- mllks4 <- mllks5 <- rep(NA, 10)

N = 2  
for (i in 1:10){
    stepMean <- runif(N, 1, 400)
    stepSD <- runif(N, 1, 400)
    stepPar <- c(stepMean, stepSD)
    mods2[[i]] <- fitHMM(data, nbStates = N, stepPar0 = stepPar, 
                        angleDist = "none", stationary = T)
    mllks2[i] <- mods2[[i]]$mod$minimum
    print(i)
}

N = 3
for (i in 1:10){
  stepMean <- runif(N, 1, 400)
  stepSD <- runif(N, 1, 400)
  stepPar <- c(stepMean, stepSD)
  mods3[[i]] <- fitHMM(data, nbStates = N, stepPar0 = stepPar, 
                       angleDist = "none", stationary = T)
  mllks3[i] <- mods3[[i]]$mod$minimum
  print(i)
}

N = 4  
for (i in 1:10){
  stepMean <- runif(N, 1, 400)
  stepSD <- runif(N, 1, 400)
  stepPar <- c(stepMean, stepSD)
  mods4[[i]] <- fitHMM(data, nbStates = N, stepPar0 = stepPar, 
                       angleDist = "none", stationary = T)
  mllks4[i] <- mods4[[i]]$mod$minimum
  print(i)
}

N = 5  
for (i in 1:10){
  stepMean <- runif(N, 1, 400)
  stepSD <- runif(N, 1, 400)
  stepPar <- c(stepMean, stepSD)
  mods5[[i]] <- fitHMM(data, nbStates = N, stepPar0 = stepPar, 
                       angleDist = "none", stationary = T)
  mllks5[i] <- mods5[[i]]$mod$minimum
  print(i)
}

# compare minimum mllk values of the 10 iterations with the ones we found before
mllks2; min(mllks2); mods[[2]]$mod$minimum
mllks3; min(mllks3); mods[[3]]$mod$minimum
mllks4; min(mllks4); mods[[4]]$mod$minimum
mllks5; min(mllks5); mods[[5]]$mod$minimum
# all good, we can continue with the models we had before

## d) model selection

# look at models, compare them visually
# with moveHMM, you can plot the results:
plot(mods[[2]])
plot(mods[[3]])
plot(mods[[4]])
plot(mods[[5]])


# calculate AIC/BIC and identify lowest.
# AIC
aic <- rep(NA, 5)

# calculate AIC for all using a for-loop
for (i in 2:5) {
  aic[i] <- 2*mods[[i]]$mod$minimum + 2*length(mods[[i]]$mod$estimate)
}
# you can double check your calculations by comparing with the AIC value given by moveHMM
AIC(mods[[2]]); aic[2]
# find minimum
which.min(aic)


# BIC
bic <- rep(NA, 5)
T <- dim(data)[1]

# calculate BIC for all using a for-loop
for (i in 2:5) {
  bic[i] <- 2*mods[[i]]$mod$minimum + log(T)*length(mods[[i]]$mod$estimate)
}

# find minimum
which.min(bic)

# let's continue with the 3- and 4-state model

## e) residual analysis
mod3 <- mods[[3]]
mod4 <- mods[[4]]

# either use the plotting residual function
plotPR(mod3)
plotPR(mod4)

# or calculate the residuals with the function and plot manually
pr3 <- pseudoRes(mod3)
pr4 <- pseudoRes(mod4)
par(mfrow=c(2,1))
plot(pr3$stepRes,ylab="pseudo-residuals",main="3 states")
plot(pr4$stepRes,ylab="pseudo-residuals",main="4 states")
# indication of daily patterns
acf(pr3$stepRes,lag.max=100,na.action=na.pass,main="3 states")
acf(pr4$stepRes,lag.max=100,na.action=na.pass,main="4 states")


