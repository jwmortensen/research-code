library(LatticeKrig)
library(FNN)
library(MASS)
library(utils)
library(MCMCpack)
library(Rcpp)
library(RcppArmadillo)

source("AMCMCUpdate.R")
load("./RData/Spatial911PtPtrn.RData")
load("./RData/MortalityTempDataNoMiss.RData")

# Variables used throughout
num.pred.locs <- nrow(pp.grid)

# Break 1428 grid points into B blocks in order to speed up updating
num.blocks <- 21  # Chose this because it factors evenly into 1428
# Choose points to act as approximate centers for each block
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
death.locs <- data.frame(lat=death$D_LAT, lon=death$D_LONG)
nn <- get.knnx(pp.grid, death.locs, 1)
nn.index <- nn$nn.index

# Get dates where a death occurred
dates <- as.Date(paste(death$Year, death$Month, death$Day), "%Y %m %d")
unique.dates <- sort(unique(dates))
# Make sure the order is the same for both
temp.data.nomiss <- temp.data.nomiss[as.character(unique.dates)]
date.ind <- cbind(unique.dates, 1:length(unique.dates))

## Merge to get lambda indices. Note that there are 1384 unique out of 1389 total so 
## there is little to be gained from combining them.
colnames(date.ind) <- c("dates", "date.index")
loc.date <- cbind(nn.index, dates)
colnames(loc.date) <- c("location.index","dates")
merged <- merge(date.ind, loc.date)
date.loc.ind <- as.matrix(merged[,2:3])

# Calculate distance between elements for use in Gaussian process
space.dist <- rdist(pp.grid)
time.dist <- rdist(0:3)

# Fix nu
nu <- 3.5

# Generate matern matrices for varying values of phi
lambda.phi <- 500
lambda.matern <- Matern(space.dist, alpha=lambda.phi, nu=nu) 
lambda.inverse.matern <- solve(lambda.matern)

# Fix this because it wasn't able to estimate it from the data.
beta.phi <- 5
beta.matern <- Matern(time.dist, alpha=beta.phi, nu=nu)
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
amcmc.it <- 500
beta.amcmc.it <- 1000

# Functions used throughout
FacToNum <- function(x) {
  as.numeric(as.character(x))
}
sourceCpp("./estimateParameters.cpp")


# Metropolis within Gibbs sampler to estimate lambda and beta coefficients
# Gibbs sampler to estimate sigma2 and lambda
MHGibbs <- function(ndraws, thin.factor, init.lambda, init.beta) {
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
  beta <- matrix(NA, nrow=ndraws, ncol=num.lags)
  
  # create initial proposal and initialize variables 
  lambda.star[1, ] <- init.lambda
  beta[1, ] <- init.beta
  temp.lstar <- init.lambda
  temp.lvar <- 0.01
  lambda.var.a <- 0.01
  lambda.var.b <- 0.01
  temp.beta <- init.beta
  temp.bvar <- 0.1
  beta.var.a <- 0.01
  beta.var.b <- 0.01
  

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

    # Here we implement metropolis hastings to get draws for lambda
    new.lstar <- temp.lstar
    for (j in 1:num.blocks) {
      # If enough iterations have been completed, begin using adaptive MCMC techniques
      prop.var.const <- 1e-3
      if (i < amcmc.it) {
        prop.var <- prop.var.const*diag(num.pred.locs/num.blocks)
      } else {
       prop.var <- (2.4^2/(num.pred.locs/num.blocks))*
         (prop.var.const*diag(num.pred.locs/num.blocks)+amcmc[[j]]$var)
      }

      prop.lstar <- mvrnorm(1, new.lstar[close.pts.index[[j]]], prop.var) 
      prop.lstar.vec <- new.lstar
      prop.lstar.vec[close.pts.index[[j]]] <- prop.lstar
      
      log.MH <- LogLike(H, prop.lstar.vec, temp.beta, date.loc.ind) - LogLike(H, new.lstar, temp.beta, date.loc.ind) 
      log.MH <- log.MH + LogLambdaPrior(prop.lstar.vec, temp.lvar, lambda.inverse.matern) - 
        LogLambdaPrior(new.lstar, temp.lvar, lambda.inverse.matern)
      
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
    
    # Metropolis hastings to get draws for beta
    prop.var.const <- 1e-6
    if (i < beta.amcmc.it) {
      prop.var <- prop.var.const*diag(num.lags)
    } else {
      prop.var <-  (2.4^2/num.lags)*
        (prop.var.const*diag(num.lags)+beta.amcmc$var)
    }
    prop.beta <- mvrnorm(1, temp.beta, prop.var)
    log.MH <- LogLike(H, temp.lstar, prop.beta, date.loc.ind) - LogLike(H, temp.lstar, temp.beta, date.loc.ind)
    log.MH <- log.MH + LogBetaPrior(prop.beta, temp.bvar, beta.inverse.matern) -
      LogBetaPrior(temp.beta, temp.bvar, beta.inverse.matern)
    temp.beta <- if(log(runif(1)) < log.MH) prop.beta else temp.beta
    if (update) beta[n.draw, ] <- temp.beta

    new.amcmc <- AMCMC.update(temp.beta, 
                              beta.amcmc$mn, beta.amcmc$var, i-1)
    beta.amcmc <- new.amcmc

    setTxtProgressBar(pb, i)
  }
 
  return(list(delta=delta, lambda.star=lambda.star, beta=beta, race=race.draws, gender=gender.draws, age=age.draws))
}

init.beta <- c(0, 0, 0, 0)
init.lambda <- rep(log(1/num.pred.locs), num.pred.locs)

Rprof()
time <- system.time(draws2000 <- MHGibbs(2000, 1, init.lambda, init.beta))
Rprof(NULL)
# save(draws, file="./RData/MHDrawsHI_MAX.RData")

# draws: 1e-2 1e-5
pdf("draws.pdf")
PlotOutput(draws)
dev.off()
# draws1: 1e-3 1e-6
pdf("draws1.pdf")
PlotOutput(draws1)
dev.off()
# draws2: 1e-4 1e-7
pdf("draws2.pdf")
PlotOutput(draws2)
dev.off()
