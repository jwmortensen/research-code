library(parallel)
library(LatticeKrig)
library(FNN)

source("AMCMCUpdate.R")
load("./RData/TempDataNoMiss.RData")

# Create matrix to allow for plotting of values
x.grid <- matrix(hrldas.grid[,2], nrow=125)
y.grid <- matrix(hrldas.grid[,1], nrow=125)

# Calculate distance between elements for use in Gaussian process
D <- rdist(pp.grid)

# Fix nu
nu <- 3.5

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

# Generate matern matrices for varying values of phi
phi.seq <- seq(5, 500, length=20)
matern.matrices <- mclapply(phi.seq, function(phi) { Matern(D, alpha=phi, nu=nu)}, mc.cores=16) 
inverse.matern <- mclapply(matern.matrices, solve, mc.cores=16)

# Adaptive MCMC stuff
amcmc - vector("list", length=num.blocks)
for (i in 1:num.blocks) {
  amcmc[[i]]$mn <- matrix(0, ncol=1, nrow=length(close.pts.index[[i]]))
  amcmc[[i]]$var <- matrix(0, ncol=length(close.pts.index[[i]]), nrow=length(close.pts.index[[i]]))
}
amcmc.it <- 100

# Metropolis within Gibbs sampler to estimate lambda and beta coefficients
#mh.gibbs <- function(draws, s


lambda.star <- rep(1, 1428)
beta <- rep(1,4)
vars <- c("HI_MAX","HI_MIN","T2MAX","T2MIN","SW_MIN","SW_MAX")
postfix <- c(".0",".1",".2",".3")
lagged.vars <- paste(vars[1], postfix, sep="")

calc.log.lambda <- function(H, lambda.star, beta, vars) {
  lambda.star + as.matrix(H[,vars])%*%beta - log(sum(exp(lambda.star + as.matrix(H[,vars])%*%beta)))
}

log.lambda <- mclapply(temp.data.nomiss, calc.log.lambda, lambda.star=lambda.star, beta=beta, vars=lagged.vars, mc.cores=16)
