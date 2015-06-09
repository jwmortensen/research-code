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
sourceCpp("./estimateParametersConstTemp.cpp")

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
nn <- get.knnx(pp.grid, call.locs, 1)
nn.index <- nn$nn.index
Nk <- rep(0, num.pred.locs)
for (i in nn.index) Nk[i] <- Nk[i] + 1


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
merged.mat <- as.matrix(merged[,2:3])

# Calculate distance between elements for use in Gaussian process
space.dist <- rdist(pp.grid)
time.dist <- rdist(0:3)

# Fix nu
nu <- 3.5

# Generate matern matrices for varying values of phi
lambda.phi <- 500
lambda.matern <- Matern(space.dist, alpha=lambda.phi, nu=nu) 
lambda.inverse.matern <- solve(lambda.matern)

# Adaptive MCMC stuff
amcmc <- vector("list", length=num.blocks)
for (i in 1:num.blocks) {
  amcmc[[i]]$mn <- matrix(0, ncol=1, nrow=length(close.pts.index[[i]]))
  amcmc[[i]]$var <- matrix(0, ncol=length(close.pts.index[[i]]), 
                           nrow=length(close.pts.index[[i]]))
}
amcmc.it <- 100

D <- function(lambda.star) {
  -2 * LogLike(lambda.star, Nk)
}


# Metropolis within Gibbs sampler to estimate lambda and beta coefficients
# Gibbs sampler to estimate sigma2 and lambda
MHGibbs <- function(ndraws, thin.factor, init.lambda, init.beta, temp.var) {
  # define variables to use throughout mh.gibbs
  vars <- c("HI_MAX","HI_MIN","T2MAX","T2MIN","SW_MIN","SW_MAX")
  var.index <- switch(temp.var,
                      HI_MAX = 1,
                      HI_MIN = 2,
                      T2MAX = 3,
                      T2MIN = 4,
                      SW_MIN = 5,
                      SW_MAX = 6)
  # for now, this only evaluates for HI_MAX
  H <- lapply(temp.data.nomiss, function(x) { as.matrix(x[,vars[var.index]]) })
  temp.means <- apply(sapply(H, function(x) { apply(x, 1, mean) }), 1, mean)
  temp.means <- cbind(1, temp.means)

  pb <- txtProgressBar(min=0, max=ndraws*thin.factor, style=3)
  n <- ifelse(thin.factor == 1, 1, 0)
  
  # initialize containers to hold drawsta
  lambda.star <- matrix(NA, nrow=ndraws, ncol=num.pred.locs)
  lvar <- numeric(ndraws)
  beta <- matrix(NA, nrow=ndraws, ncol=2)
  d.vals <- numeric(ndraws)
#   bvar <- numeric(ndraws)
  
  # create initial proposal and initialize variables 
  lambda.star[1, ] <- temp.lstar <- init.lambda
  beta[1, ] <- temp.beta <- init.beta
  d.vals[1] <- D(temp.lstar)
  temp.lvar <- 0.01
  lambda.var.a <- 0.01
  lambda.var.b <- 0.01
  
  
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
    # Get draws for lambda.var using the complete conditional
    lambda.a <- lambda.var.a + num.pred.locs/2
    lambda.b <- lambda.var.b + (1/2)*t(temp.lstar)%*%
      lambda.inverse.matern%*%(temp.lstar) 
    temp.lvar <- 1/rgamma(1, shape=lambda.a, rate=lambda.b)
    
    # Here we implement metropolis hastings to get draws for lambda
    new.lstar <- temp.lstar
    for (j in 1:num.blocks) {
      # If enough iterations have been completed, begin using adaptive MCMC techniques
      prop.var.const <- 5e-2
      if (i < amcmc.it) {
        prop.var <- prop.var.const*diag(num.pred.locs/num.blocks)
      } else {
        prop.var <- (2.4^2/(num.pred.locs/num.blocks))*
          (prop.var.const*diag(num.pred.locs/num.blocks)+amcmc[[j]]$var)
      }
      
      prop.lstar <- mvrnorm(1, new.lstar[close.pts.index[[j]]], prop.var) 
      prop.lstar.vec <- new.lstar
      prop.lstar.vec[close.pts.index[[j]]] <- prop.lstar
      
      log.MH <- LogLike(prop.lstar.vec, Nk) - LogLike(new.lstar, Nk) 
      log.MH <- log.MH + LogLambdaPrior(prop.lstar.vec, temp.lvar, lambda.inverse.matern, temp.means, temp.beta) - 
        LogLambdaPrior(new.lstar, temp.lvar, lambda.inverse.matern, temp.means, temp.beta)
      
      if (log(runif(1)) < log.MH) {
        new.lstar[close.pts.index[[j]]] <- prop.lstar
      }
      new.amcmc <- amcmcUpdate(new.lstar[close.pts.index[[j]]], 
                                amcmc[[j]]$mn, amcmc[[j]]$var, i-1)
      amcmc[[j]] <- new.amcmc
    }
    temp.lstar <- new.lstar
    
    # Metropolis hastings to get draws for beta
    s2 <- 100
    beta.var <- solve(1/s2 * diag(2) + t(temp.means)%*%lambda.inverse.matern%*%temp.means)
    beta.mean <- (1/temp.lvar)*beta.var%*%t(temp.means)%*%lambda.inverse.matern%*%temp.lstar
    temp.beta <- mvrnorm(1, beta.mean, beta.var)
    
    if (i %% thin.factor == 0) {
      n <- n + 1
      lvar[n] <- temp.lvar
      lambda.star[n, ] <- temp.lstar
      beta[n, ] <- temp.beta
      d.vals[n] <- D(temp.lstar)
    }
    
    setTxtProgressBar(pb, i)
  }
  
  return(list(delta=delta, lambda.star=lambda.star, lstar.var=lvar, beta=beta, race=race.draws, gender=gender.draws, age=age.draws, d.vals=d.vals))
}

init.beta <- c(-0.2,0)
init.lambda <- rep(log(1/num.pred.locs), num.pred.locs)

DIC <- function(draws) {
  lambda.star <- apply(draws$lambda.star, 2, mean)
  2 * mean(draws$d.vals) - D(lambda.star)
}

temp.vars <- c("HI_MAX", "HI_MIN", "T2MAX", "T2MIN", "SW_MAX", "SW_MIN")
draws <- lapply(temp.vars, function(x) { MHGibbs(100, 50, init.lambda, init.beta, x) })
DIC.vals <- lapply(draws, DIC)
names(draws) <- temp.vars
names(DIC.vals) <- temp.vars
DIC.vals
