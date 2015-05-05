library(parallel)
library(LatticeKrig)
library(FNN)

source("AMCMCUpdate.R")
load("./RData/Spatial911PtPtrn.RData")
load("./RData/TempDataNoMiss.RData")

# Create matrix to allow for plotting of values
x.grid <- matrix(hrldas.grid[,2], nrow=125)
y.grid <- matrix(hrldas.grid[,1], nrow=125)

# Get number of prediction locations to use throughout
K <- nrow(pp.grid)

# Break 1428 grid points into B blocks in order to speed up updating
num.blocks <- 51 # Chose this because it factors evenly into 1428
b.pts <- pp.grid[seq.int(from=1, to=K, length=num.blocks),]

close.pts <- vector("list", length=num.blocks)
close.pts.index <- vector("list", length=num.blocks)
pp.copy <- pp.grid

for (i in 1:num.blocks) {
  nn <- get.knnx(pp.copy, b.pts[i,], K/num.blocks)
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
SPACE.DIST <- rdist(pp.grid)
TIME.DIST <- rdist(0:3)

# Fix nu
nu <- 3.5

# Generate matern matrices for varying values of phi
lambda.phi.seq <- seq(5, 500, length=16)
lambda.matern.matrices <- mclapply(lambda.phi.seq, function(phi) { Matern(SPACE.DIST, alpha=phi, nu=nu)}, mc.cores=16) 
lambda.inverse.matern <- mclapply(lambda.matern.matrices, solve, mc.cores=16)

beta.phi.seq <- seq(0,16,length=16)
beta.matern.matrices <- mclapply(beta.phi.seq, function(phi) { Matern(TIME.DIST, alpha=phi, nu=nu) }, mc.cores=16)
beta.inverse.matern <- mclapply(beta.matern.matrices, solve, mc.cores=16)

# Adaptive MCMC stuff
amcmc - vector("list", length=num.blocks)
for (i in 1:num.blocks) {
  amcmc[[i]]$mn <- matrix(0, ncol=1, nrow=length(close.pts.index[[i]]))
  amcmc[[i]]$var <- matrix(0, ncol=length(close.pts.index[[i]]), nrow=length(close.pts.index[[i]]))
}
amcmc.it <- 100

# Metropolis within Gibbs sampler to estimate lambda and beta coefficients
# Gibbs sampler to estimate sigma2 and lambda
mh.gibbs <- function(ndraws, lambda.var.start, lambda.var.a, lambda.var.b, beta.var.start, beta.var.a, beta.var.b) {
  # define variables to use throughout mh.gibbs
  num.lags <- 4
  vars <- c("HI_MAX","HI_MIN","T2MAX","T2MIN","SW_MIN","SW_MAX")
  postfix <- c(".0",".1",".2",".3")
  # for now, this only evaluates for HI_MAX
  lagged.vars <- paste(vars[1], postfix, sep="")
  
  # initialize containers to hold draws
  lambda.star <- matrix(NA, nrow=ndraws, ncol=K)
  lambda.var <- numeric(ndraws)
  beta <- matrix(NA, nrow=ndraws, ncol=num.lags)
  beta.var <- numeric(ndraws)
  lambda.phi <- numeric(ndraws)
  beta.phi <- numeric(ndraws)
  
  # create initial proposal and initialize variables 
  lambda.star[1,] <- rep(1/K, K)
  lambda.var[1] <- lambda.var.start
  lambda.phi.ind <- sample(1:length(lambda.phi.seq), 1)
  lambda.phi[1] <- lambda.phi.seq[lambda.phi.ind]
  beta[1,] <- rep(1,num.lags)
  beta.var[1] <- beta.var.start
  beta.phi.ind <- sample(1:length(beta.phi.seq), 1)
  beta.phi[1] <- beta.phi.seq[beta.phi.ind]
  
  # Initialize functions for use in M-H
  calc.log.lambda <- function(H, lambda.star, beta, vars) {
    lambda.star + as.matrix(H[,vars])%*%beta - log(sum(exp(lambda.star + as.matrix(H[,vars])%*%beta)))
  }
  
  get.log.lambda <- function(x) {
    log.lambda[[x["date.index"]]][x["location.index"]]
  }
  
  log.like <- function(lambda.star, beta) {
    log.lambda <- mclapply(temp.data.nomiss, calc.log.lambda, lambda.star=lambda.star, beta=beta, vars=lagged.vars, mc.cores=1)
    obs.log.lambda <- apply(merged, 1, get.log.lambda)
    sum(obs.log.lambda)
  }
  
  log.lambda.prior <- function(lambda.star,sig2,matern.ind){
    -0.5*(t(lambda.star)%*%lambda.inverse.matern[[matern.ind]]%*%lambda.star)/sig2
  }
  
  log.beta.prior <- function(beta, sig2, matern.ind) {
    -0.5*(t(beta)%*%beta.inverse.matern[[matern.ind]]%*%beta)/sig2
  }
  
  lambda.phi.post.dens <- function(lambda.star, sig2, matern.ind) {
    -log(length(lambda.phi.seq)) + log.lambda.prior(lambda.star, sig2, matern.ind)
  }
  
  beta.phi.post.dens <- function(beta, sig2, matern.ind) {
    -log(length(beta.phi.seq)) + log.beta.prior(beta, sig2, matern.ind)
  }
  
  # delta is conjugate so we just draw it up front
  delta <- rgamma(ndraws, shape=nrow(calls)+0.001, rate=1.001)
  
  for (i in 2:ndraws) {
    
    # Get draws for lambda.var using the complete conditional
    lambda.a <- lambda.var.a + K/2
    lambda.b <- lambda.var.b + (1/2)*t(lambda.star[i-1,])%*%lambda.inverse.matern[[lambda.phi.ind]]%*%(lambda.star[i-1,]) 
    lambda.var[i] <- 1/rgamma(1, shape=lambda.a, rate=lambda.b)
    
    # Get draws for beta.var using the complete conditional
    beta.a <- beta.var.a + num.lags/2
    beta.b <- beta.var.b + (1/2)*t(beta[i-1,])%*%beta.inverse.matern[[beta.phi.ind]]%*%(beta[i-1,])
    beta.var[i] <- 1/rgamma(1, shape=beta.a, rate=beta.b)
    
    # Here we implement metropolis hastings to get draws for lambda
    lambda.star[i,] <- lambda.star[i-1,]
    for (j in 1:num.blocks) {
      # If enough iterations have been completed, begin using adaptive MCMC techniques
      if (i < amcmc.it) {
        prop.var <- 0.01*diag(K/num.blocks)
      }
      else {
        prop.var <- (2.4^2/(K/num.blocks))*(0.01*diag(K/num.blocks)+amcmc[[j]]$var)
      }
      
      prop.lstar <- mvrnorm(1, lambda.star[i,close.pts.index[[j]]], prop.var)
      prop.lstar.vec <- lambda.star[i,]
      prop.lstar.vec[close.pts.index[[j]]] <- prop.lstar
      prop.lvec <- exp(prop.lstar.vec)/sum(exp(prop.lstar.vec))
      cur.lvec <- exp(lambda.star[i,])/sum(exp(lambda.star[i,]))
      
      log.MH <- log.like(prop.lvec) - log.like(cur.lvec) 
      log.MH <- log.MH + lambda.prior(prop.lstar.vec,sigma2[i]) - lambda.prior(lambda.star[i,],sigma2[i])
      
      if (log(runif(1)) < log.MH) {
        lambda.star[i,close.pts.index[[j]]] <- prop.lstar
      }
      new.amcmc <- AMCMC.update(lambda.star[i, close.pts.index[[j]]], amcmc[[j]]$mn, amcmc[[j]]$var, i-1)
      amcmc[[j]] <- new.amcmc
    }
  }
  lambda <- exp(lambda.star) / apply(exp(lambda.star), 1, sum)  
  return(list(delta=delta, sigma2=sigma2, lambda=lambda))
}



lambda.star <- rep(1/1428, 1428)
beta <- rep(1,4)



log.lambda <- mclapply(temp.data.nomiss, calc.log.lambda, lambda.star=lambda.star, beta=beta, vars=lagged.vars, mc.cores=1)
