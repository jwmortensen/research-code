library(LatticeKrig)
library(FNN)
library(MASS)
library(utils)
library(MCMCpack)

source("AMCMCUpdate.R")
load("./RData/Spatial911PtPtrn.RData")
load("./RData/TempDataNoMiss.RData")

# Variables used throughout
num.pred.locs <- nrow(pp.grid)

# Create matrix to allow for plotting of values
x.grid <- matrix(hrldas.grid[,2], nrow=125)
y.grid <- matrix(hrldas.grid[,1], nrow=125)

# Break 1428 grid points into B blocks in order to speed up updating
num.blocks <- 51  # Chose this because it factors evenly into 1428
b.pts <- pp.grid[seq.int(from=1, to=num.pred.locs, length=num.blocks),]

close.pts <- vector("list", length=num.blocks)
close.pts.index <- vector("list", length=num.blocks)
pp.copy <- pp.grid

for (i in 1:num.blocks) {
  nn <- get.knnx(pp.copy, b.pts[i,], num.pred.locs/num.blocks)
  close.pts[[i]] <- pp.copy[nn$nn.index,]
  close.pts.index[[i]] <- which(paste(pp.grid$Latitude, pp.grid$Longitude)
    %in% paste(close.pts[[i]]$Latitude, close.pts[[i]]$Longitude))
  pp.copy <- pp.copy[-nn$nn.index,]
}

# Get nearest prediction locations for each observed call
call.locs <- calls[,1:2]
nn <- get.knnx(pp.grid, call.locs, 1)
nn.index <- nn$nn.index

# Get dates where a 911 call occurred
split.dates <- strsplit(as.character(strptime(calls$DOY, format="%j")), "-")
months <- numeric(length(split.dates))
days <- numeric(length(split.dates))
for (i in 1:length(split.dates)) {
  months[i] <- split.dates[[i]][2]
  days[i] <- split.dates[[i]][3]
}
dates <- as.Date(paste(calls$Year, months, days), "%Y %m %d")
unique.dates <- unique(dates)
date.ind <- cbind(unique.dates, 1:length(unique.dates))

## Merge to get lambda indices. Note that there are 1384 unique out of 1389 total so 
## there is little to be gained from combining them.
colnames(date.ind) <- c("dates", "date.index")
loc.date <- cbind(nn.index, dates)
colnames(loc.date) <- c("location.index","dates")
merged <- merge(date.ind, loc.date)

# Calculate distance between elements for use in Gaussian process
spaceDist <- rdist(pp.grid)
timeDist <- rdist(0:3)

# Fix nu
nu <- 3.5

# Generate matern matrices for varying values of phi
lambda.phi <- 500
lambda.matern <- Matern(spaceDist, alpha=lambda.phi, nu=nu) 
lambda.inverse.matern <- solve(lambda.matern)

# Fix this because it wasn't able to estimate it from the data.
beta.phi <- 5
beta.matern <- Matern(timeDist, alpha=beta.phi, nu=nu)
beta.inverse.matern <- solve(beta.matern)

# Adaptive MCMC stuff
amcmc <- vector("list", length=num.blocks)
for (i in 1:num.blocks) {
  amcmc[[i]]$mn <- matrix(0, ncol=1, nrow=length(close.pts.index[[i]]))
  amcmc[[i]]$var <- matrix(0, ncol=length(close.pts.index[[i]]), 
                           nrow=length(close.pts.index[[i]]))
}

# Adaptive MCMC stuff for beta
beta.amcmc <- vector("list")
num.lags <- 4
beta.amcmc$mn <- matrix(0, ncol=1, nrow=num.lags)
beta.amcmc$var <- matrix(0, ncol=num.lags, nrow=num.lags)
amcmc.it <- 100
beta.amcmc.it <- 1000

# Utility functions
FacToNum <- function(x) {
  as.numeric(as.character(x))
}

# Metropolis within Gibbs sampler to estimate lambda and beta coefficients
# Gibbs sampler to estimate sigma2 and lambda
MHGibbs <- function(ndraws, thin.factor, lambda.var.start, lambda.var.a, 
                    lambda.var.b, beta.var.start, beta.var.a, beta.var.b) {
  # define variables to use throughout mh.gibbs
  vars <- c("HI_MAX","HI_MIN","T2MAX","T2MIN","SW_MIN","SW_MAX")
  postfix <- c(".0",".1",".2",".3")
  # for now, this only evaluates for HI_MAX
  lagged.vars <- paste(vars[1], postfix, sep="")
  H <- lapply(temp.data.nomiss, function(x) { as.matrix(x[,lagged.vars]) })
  pb <- txtProgressBar(min=0, max=ndraws*thin.factor, style=3)
  n.draw <- ifelse(thin.factor == 1, 1, 0)
  
  # initialize containers to hold drawsta
  lambda.star <- matrix(NA, nrow=ndraws, ncol=num.pred.locs)
  lambda.var <- numeric(ndraws)
  beta <- matrix(NA, nrow=ndraws, ncol=num.lags)
  beta.var <- numeric(ndraws)
  
  # create initial proposal and initialize variables 
  lambda.star[1, ] <- rep(1/num.pred.locs, num.pred.locs)
  lambda.var[1] <- lambda.var.start
  beta[1, ] <- c(-0.15, 0.15, 0.08, 0)
  beta.var[1] <- beta.var.start
  temp.lstar <- lambda.star[1, ]
  temp.lvar <- lambda.var[1]
  temp.beta <- beta[1,]
  temp.bvar <- beta.var[1]
  
  
  # Initialize functions for use in M-H
  CalcLogLambda <- function(H, lambda.star, beta) {
    lambda.star + H%*%beta - 
      log(sum(exp(lambda.star + H%*%beta)))
  }
  
  LogLike <- function(lambda.star, beta) {
    log.lambda <- lapply(H, CalcLogLambda, 
                           lambda.star=lambda.star, beta=beta)
    GetLogLambda <- function(x) {
      log.lambda[[x["date.index"]]][x["location.index"]]
    }
    obs.log.lambda <- apply(merged, 1, GetLogLambda)
    sum(obs.log.lambda)
  }
  
  LogLambdaPrior <- function(lambda.star, sig2) {
    -0.5*(t(lambda.star)%*%lambda.inverse.matern%*%lambda.star)/sig2
  }
  
  LogBetaPrior <- function(beta, sig2) {
    -0.5*(t(beta)%*%beta.inverse.matern%*%beta)/sig2
  }
  
  # delta is conjugate so we just draw it up front
  delta <- rgamma(ndraws, shape=nrow(calls)+0.001, rate=1.001)
  
  # fit age with a dirichlet prior, which is conjugate
  age.alpha <- 2
  count.ages <- data.frame(table(calls$Age))
  names(count.ages) <- c("age", "count")
  count.ages$age <- FacToNum(count.ages$age)
  count.ages$count <- FacToNum(count.ages$count)
  n.ages <- merge(data.frame(age=0:100), count.ages, all.x=TRUE)
  n.ages$count[is.na(n.ages$count)] <- 0
  age.draws <- rdirichlet(ndraws, n.ages$count + age.alpha)
  
  # fit gender with a dirichlet prior
  gender.alpha <- 2
  n.gender <- data.frame(table(calls$Gender))
  names(n.gender) <- c("gender", "count")
  gender.draws <- rdirichlet(ndraws, n.gender$count + gender.alpha)

  # fit race with a dirichlet prior
  race.alpha <- 2
  n.race <- data.frame(table(calls$Eth))
  names(n.race) <- c("race", "count")
  race.draws <- rdirichlet(ndraws, n.race$count + race.alpha)
  
  for (i in 2:(ndraws*thin.factor)) {
    update <- ifelse(i %% thin.factor == 0, TRUE, FALSE)
    if (update) n.draw <- n.draw + 1

    # Get draws for lambda.var using the complete conditional
    lambda.a <- lambda.var.a + num.pred.locs/2
    lambda.b <- lambda.var.b + (1/2)*t(temp.lstar)%*%
      lambda.inverse.matern%*%(temp.lstar) 
    temp.lvar <- 1/rgamma(1, shape=lambda.a, rate=lambda.b)
    if (update) lambda.var[n.draw] <- temp.lvar

    # Here we implement metropolis hastings to get draws for lambda
    new.lstar <- temp.lstar
    for (j in 1:num.blocks) {
      # If enough iterations have been completed, begin using adaptive MCMC techniques
      prop.var.const <- 1e-8
      if (i < amcmc.it) {
        prop.var <- prop.var.const*diag(num.pred.locs/num.blocks)
      } else {
       prop.var <- (2.4^2/(num.pred.locs/num.blocks))*
         (prop.var.const*diag(num.pred.locs/num.blocks)+amcmc[[j]]$var)
      }
      
      prop.lstar <- mvrnorm(1, new.lstar[close.pts.index[[j]]], prop.var) 
      prop.lstar.vec <- new.lstar
      prop.lstar.vec[close.pts.index[[j]]] <- prop.lstar

      log.MH <- LogLike(prop.lstar.vec, temp.beta) - LogLike(new.lstar, temp.beta) 
      log.MH <- log.MH + LogLambdaPrior(prop.lstar.vec, temp.lvar) - 
        LogLambdaPrior(new.lstar, temp.lvar)
      
      if (log(runif(1)) < log.MH) {
        new.lstar[close.pts.index[[j]]] <- prop.lstar
      }
      new.amcmc <- AMCMC.update(new.lstar[close.pts.index[[j]]], 
                                amcmc[[j]]$mn, amcmc[[j]]$var, i-1)
      amcmc[[j]] <- new.amcmc
    }
    temp.lstar <- new.lstar
    if (update) lambda.star[n.draw, ] <- temp.lstar
    
    # Get draws for beta.var using the complete conditional
    beta.a <- beta.var.a + num.lags/2
    beta.b <- beta.var.b + (1/2)*t(temp.beta)%*%
      beta.inverse.matern%*%(temp.beta)
    temp.bvar <- 1/rgamma(1, shape=beta.a, rate=beta.b)
    if (update) beta.var[n.draw] <- temp.bvar
    
    # Metropolis hastings to get draws for beta
    prop.var.const <- 5e-6
    if (i < beta.amcmc.it) {
      prop.var <- prop.var.const*diag(num.lags)
    } else {
      prop.var <-  (2.4^2/num.lags)*
        (prop.var.const*diag(num.lags)+beta.amcmc$var)
    }
    prop.beta <- mvrnorm(1, temp.beta, prop.var)
    log.MH <- LogLike(temp.lstar, prop.beta) - LogLike(temp.lstar, temp.beta)
    log.MH <- log.MH + LogBetaPrior(prop.beta, temp.bvar) -
      LogBetaPrior(temp.beta, temp.bvar)
    temp.beta <- if(log(runif(1)) < log.MH) prop.beta else temp.beta
    if (update) beta[n.draw, ] <- temp.beta
    
    new.amcmc <- AMCMC.update(temp.beta, 
                              beta.amcmc$mn, beta.amcmc$var, i-1)
    beta.amcmc <- new.amcmc

    setTxtProgressBar(pb, i)
  }
 
  return(list(delta=delta, lambda.star=lambda.star, beta=beta, race=race.draws, gender=gender.draws, age=age.draws))
}

Rprof()
time <- system.time(draws.new <- MHGibbs(5, 1, 0.01, 0.01, 0.01, 1, 0.01, 0.01))
Rprof(NULL)
# save(draws, file="./RData/MHDrawsHI_MAX.RData")

