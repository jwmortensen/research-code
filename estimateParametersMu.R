library(parallel)
library(LatticeKrig)
library(FNN)
library(MASS)
library(utils)
library(MCMCpack)
library(Rcpp)
library(RcppArmadillo)

load("./RData/Spatial911PtPtrn.RData")
load("./RData/AllTempData.RData")
load("./RData/InterceptData.RData")
sourceCpp("./estimateParametersConstTemp.cpp")
death <- read.csv("./RData/cleanDeathData.csv", header=T)

# Variables used throughout
num.pred.locs <- nrow(pp.grid)

# Functions used throughout
FacToNum <- function(x) {
  as.numeric(as.character(x))
}

# Break 1428 grid points into B blocks in order to speed up updating
num.blocks <- 21  # Chose this because it factors evenly into 1428
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
call.nn <- get.knnx(pp.grid, call.locs, 1)
call.index <- call.nn$nn.index
Nk.911 <- numeric(num.pred.locs)
for (i in call.index) Nk.911[i] <- Nk.911[i] + 1

# Get nearest prediction locations for each death
death.locs <- data.frame(lat=death$D_LAT, lon=death$D_LONG)
death.nn <- get.knnx(pp.grid, death.locs, 1)
death.index <- death.nn$nn.index
Nk.death <- numeric(num.pred.locs)
for (i in death.index) Nk.death[i] <- Nk.death[i] + 1

Nk.mu <- Nk.911 + Nk.death

# Calculate distance between elements for use in Gaussian process
space.dist <- rdist(pp.grid)

# Fix nu
nu <- 3.5

# Generate matern matrices for varying values of phi
lambda.phi <- 500
lambda.matern <- Matern(space.dist, alpha=lambda.phi, nu=nu)
chol.matern <- chol(lambda.matern)
lambda.inverse.matern <- solve(lambda.matern)

# Adaptive MCMC stuff
amcmc.911 <- vector("list", length=num.blocks)
amcmc.death <- vector("list", length=num.blocks)
amcmc.mu <- vector("list", length=num.blocks)
for (i in 1:num.blocks) {
  amcmc.911[[i]]$mn <- matrix(0, ncol=1, nrow=length(close.pts.index[[i]]))
  amcmc.911[[i]]$var <- matrix(0, ncol=length(close.pts.index[[i]]), 
                           nrow=length(close.pts.index[[i]]))
  amcmc.death[[i]]$mn <- matrix(0, ncol=1, nrow=length(close.pts.index[[i]]))
  amcmc.death[[i]]$var <- matrix(0, ncol=length(close.pts.index[[i]]), 
                               nrow=length(close.pts.index[[i]]))
  amcmc.mu[[i]]$mn <- matrix(0, ncol=1, nrow=length(close.pts.index[[i]]))
  amcmc.mu[[i]]$var <- matrix(0, ncol=length(close.pts.index[[i]]), 
                                 nrow=length(close.pts.index[[i]]))
}
amcmc.it <- 100

D <- function(lambda.star) {
  -2 * LogLike(lambda.star, Nk)
}

DIC <- function(draws) {
  lambda.star <- apply(draws$lambda.star, 2, mean)
  2 * mean(draws$d.vals) - D(lambda.star)
}

MHGibbs <- function(ndraws, thin.factor, init.lambda, init.beta, temp.index) {
  # define variables to use throughout mh.gibbs
  temp.var <- switch(temp.index,
                     "1" = "HI_MAX",
                     "2" = "HI_MIN",
                     "3" = "T2MAX",
                     "4" = "T2MIN",
                     "5" = "SW_MIN",
                     "6" = "SW_MAX")
  intercept <- as.matrix(cbind(1, intercept.df[c(temp.var, "Population", "PercentOver65")]))
  
  pb <- txtProgressBar(min=0, max=ndraws*thin.factor, style=3)
  n <- ifelse(thin.factor == 1, 1, 0)
  
  # initialize containers to hold drawsta
  lstar.mu <- lstar.911 <- lstar.death <- lstar.avg <- matrix(NA, nrow=ndraws, ncol=num.pred.locs)
  lvar.mu <- lvar.911 <- lvar.death <- numeric(ndraws)
  beta <- matrix(NA, nrow=ndraws, ncol=ncol(intercept))
  d.vals <- numeric(ndraws)
  #   bvar <- numeric(ndraws)
  
  # create initial proposal and initialize variables 
  lstar.mu[1, ] <- temp.lstar.mu <- init.lambda
  lstar.911[1, ] <- temp.lstar.911 <- init.lambda #- log(1/num.pred.locs)
  lstar.death[1, ] <- temp.lstar.death <- init.lambda #- log(1/num.pred.locs)
  beta[1, ] <- temp.beta <- init.beta
  #   d.vals[1] <- D(temp.lstar)
  temp.lvar.mu <- temp.lvar.911 <- temp.lvar.death <- 50
#   temp.lvar.mu <- 5000
  lambda.var.a <- 2
  lambda.var.b <- 1
  
  
  # delta is conjugate so we just draw it up front
  delta.911 <- rgamma(ndraws, shape=nrow(calls)+0.001, rate=1.001)
  delta.death <- rgamma(ndraws, shape=nrow(death)+0.001, rate=1.001)
  
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
    # Get draws for lvar.911 using the complete conditional
        lambda.a <- lambda.var.a + num.pred.locs/2
        lambda.b <- lambda.var.b + (1/2)*t(temp.lstar.911)%*%
          lambda.inverse.matern%*%(temp.lstar.911) 
        temp.lvar.911 <- 1/rgamma(1, shape=lambda.a, rate=lambda.b)
    
    # Here we implement metropolis hastings to get draws for lstar.911
    new.lstar.911 <- temp.lstar.911
    for (j in 1:num.blocks) {
      # If enough iterations have been completed, begin using adaptive MCMC techniques
      prop.var.const.911 <- 1e-3
      if (i < amcmc.it) {
        prop.var <- prop.var.const.911*diag(num.pred.locs/num.blocks)
      } else {
        prop.var <- (2.4^2/(num.pred.locs/num.blocks))*
          (prop.var.const.911*diag(num.pred.locs/num.blocks)+amcmc.911[[j]]$var)
      }
      
      prop.lstar.911 <- t(mvrnormC(1, new.lstar.911[close.pts.index[[j]]], prop.var))
      prop.lstar.vec.911 <- new.lstar.911
      prop.lstar.vec.911[close.pts.index[[j]]] <- prop.lstar.911
      
      log.MH <- LogLike(prop.lstar.vec.911, Nk.911) - LogLike(new.lstar.911, Nk.911) 
      log.MH <- log.MH + LogLambdaPrior(prop.lstar.vec.911, temp.lvar.911, lambda.inverse.matern) - 
        LogLambdaPrior(new.lstar.911, temp.lvar.911, lambda.inverse.matern)
      
      if (log(runif(1)) < log.MH) {
        new.lstar.911[close.pts.index[[j]]] <- prop.lstar.911
      }
      new.amcmc <- amcmcUpdate(new.lstar.911[close.pts.index[[j]]], 
                               amcmc.911[[j]]$mn, amcmc.911[[j]]$var, i-1)
      amcmc.911[[j]] <- new.amcmc
    }
    temp.lstar.911 <- new.lstar.911
    
    # Get draws for lvar.death using the complete conditional
        lambda.a <- lambda.var.a + num.pred.locs/2
        lambda.b <- lambda.var.b + (1/2)*t(temp.lstar.death)%*%
          lambda.inverse.matern%*%(temp.lstar.death) 
        temp.lvar.death <- 1/rgamma(1, shape=lambda.a, rate=lambda.b)
    
    # Here we implement metropolis hastings to get draws for lstar.death
    new.lstar.death <- temp.lstar.death
    for (j in 1:num.blocks) {
      # If enough iterations have been completed, begin using adaptive MCMC techniques
      prop.var.const <- 2e-4
      if (i < amcmc.it) {
        prop.var <- prop.var.const.death*diag(num.pred.locs/num.blocks)
      } else {
        prop.var <- (2.4^2/(num.pred.locs/num.blocks))*
          (prop.var.const.death*diag(num.pred.locs/num.blocks)+amcmc.death[[j]]$var)
      }
      
      prop.lstar.death <- t(mvrnormC(1, new.lstar.death[close.pts.index[[j]]], prop.var)) 
      prop.lstar.vec.death <- new.lstar.death
      prop.lstar.vec.death[close.pts.index[[j]]] <- prop.lstar.death
      
      log.MH <- LogLike(prop.lstar.vec.death, Nk.death) - LogLike(new.lstar.death, Nk.death) 
      log.MH <- log.MH + LogLambdaPrior(prop.lstar.vec.death, temp.lvar.death, lambda.inverse.matern) - 
        LogLambdaPrior(new.lstar.death, temp.lvar.death, lambda.inverse.matern)
      
      if (log(runif(1)) < log.MH) {
        new.lstar.death[close.pts.index[[j]]] <- prop.lstar.death
      }
      new.amcmc <- amcmcUpdate(new.lstar.death[close.pts.index[[j]]], 
                               amcmc.death[[j]]$mn, amcmc.death[[j]]$var, i-1)
      amcmc.death[[j]] <- new.amcmc
    }
    temp.lstar.death <- new.lstar.death
    
    # Get draws for lvar.mu using the complete conditional
        lambda.a <- lambda.var.a + num.pred.locs/2
        lambda.b <- lambda.var.b + (1/2)*t(temp.lstar.mu - intercept%*%temp.beta)%*%
          lambda.inverse.matern%*%(temp.lstar.mu - intercept%*%temp.beta) 
        temp.lvar.mu <- 1/rgamma(1, shape=lambda.a, rate=lambda.b)
    
    # Here we implement metropolis hastings to get draws for lstar.death
    new.lstar.mu <- temp.lstar.mu
    for (j in 1:num.blocks) {
      # If enough iterations have been completed, begin using adaptive MCMC techniques
      prop.var.const.mu <- 2e-4
      if (i < amcmc.it) {
        prop.var <- prop.var.const.mu*diag(num.pred.locs/num.blocks)
      } else {
        prop.var <- (2.4^2/(num.pred.locs/num.blocks))*
          (prop.var.const.mu*diag(num.pred.locs/num.blocks)+amcmc.mu[[j]]$var)
      }
      
      prop.lstar.mu <- t(mvrnormC(1, new.lstar.mu[close.pts.index[[j]]], prop.var)) 
      prop.lstar.vec.mu <- new.lstar.mu
      prop.lstar.vec.mu[close.pts.index[[j]]] <- prop.lstar.mu
      
      log.MH <- LogLike(prop.lstar.vec.mu, Nk.mu) - LogLike(new.lstar.mu, Nk.mu) 
      log.MH <- log.MH + LogLambdaMuPrior(prop.lstar.vec.mu, temp.lvar.mu, lambda.inverse.matern, intercept, temp.beta) - 
        LogLambdaMuPrior(new.lstar.mu, temp.lvar.mu, lambda.inverse.matern, intercept, temp.beta)
      
      if (log(runif(1)) < log.MH) {
        new.lstar.mu[close.pts.index[[j]]] <- prop.lstar.mu
      }
      new.amcmc <- amcmcUpdate(new.lstar.mu[close.pts.index[[j]]], 
                               amcmc.mu[[j]]$mn, amcmc.mu[[j]]$var, i-1)
      amcmc.mu[[j]] <- new.amcmc
    }
    temp.lstar.mu <- new.lstar.mu
    
    # Complete conditional to get draws for beta
    s2 <- 50
    beta.var <- solve(1/s2 * diag(ncol(intercept)) + t(intercept)%*%lambda.inverse.matern%*%intercept)
    beta.mean <- (1/temp.lvar.mu)*beta.var%*%t(intercept)%*%lambda.inverse.matern%*%temp.lstar.mu
    temp.beta <- t(mvrnormC(1, beta.mean, beta.var))
    
    if (i %% thin.factor == 0) {
      n <- n + 1
      lvar.911[n] <- temp.lvar.911
      lstar.911[n, ] <- temp.lstar.911
      lvar.death[n] <- temp.lvar.death
      lstar.death[n, ] <- temp.lstar.death
      lvar.mu[n] <- temp.lvar.mu
      lstar.mu[n, ] <- temp.lstar.mu
      beta[n, ] <- temp.beta
      d.vals[n] <- D(temp.lstar.mu)
    }
    
    setTxtProgressBar(pb, i)
  }
  
  return(list(delta.911 = delta.911, 
              lstar.911 = lstar.911, 
              lvar.911 = lvar.911, 
              delta.death = delta.death,
              lstar.death = lstar.death, 
              lvar.death = lvar.death, 
              lstar.mu = lstar.mu, 
              lvar.mu = lvar.mu, 
              beta = beta, 
              race=race.draws, 
              gender=gender.draws, 
              age=age.draws, 
              d.vals = d.vals))
}

init.beta <- c(0,0,0,0)
init.lambda <- rep(log(1/num.pred.locs), num.pred.locs)




