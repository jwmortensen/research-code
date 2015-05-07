library(foreach)
library(doParallel)
source("fitMaternGP.R")
## Brings in a list called temp.data that has NAs that need to be estimated
load("./RData/TemperatureData_withMissing.RData")

# Fix nu at 3.5
nu <- 3.5
numCores <- 24
registerDoParallel(cores=numCores)
time <- system.time(temp.data.nomiss <- foreach(i=1:length(temp.data)) %dopar% {
  for (j in 3:ncol(temp.data[[i]])) {
    where.NA <- is.na(temp.data[[i]][,j])
    unobs.dat <- temp.data[[i]][where.NA,]
    obs.dat <- temp.data[[i]][!where.NA,]
    locs <- cbind(obs.dat$LON, obs.dat$LAT)
    pred.locs <- cbind(unobs.dat$LON, unobs.dat$LAT)
    
    # Estimate parameters using the fit.Matern.GP function
    params <- fit.Matern.GP(obs.dat[,j], matrix(rep(1, length(obs.dat[,j])), ncol=1), locs, nu)
    s2 <- params$sigma2
    phi <- params$ph
    tau2 <- params$tau2
    mu <- params$beta.hat
    
    D <- rdist(rbind(pred.locs, locs))
    K <- nrow(pred.locs)
    N <- nrow(locs)
    
    V <- s2*Matern(D,alpha=phi,nu=nu) + tau2*diag(nrow(D))
    predictions <- mu + V[1:K, K+(1:N)]%*%solve(V[K+(1:N), K+(1:N)])%*%(obs.dat[,j]-mu)
    temp.data[[i]][where.NA,j] <- predictions
  }
  temp.data[[i]]
})

# Took 8.5 hours to do 576 days with 6 variables on 24 cores
names(temp.data.nomiss) <- names(temp.data)
save(temp.data.nomiss, file="./RData/TempDataNoMiss.RData")
