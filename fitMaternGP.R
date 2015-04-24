library(LatticeKrig)

fit.Matern.GP <- function(y,X,s,nu){
	
  ## Create a Sequence for Spatial Range
  D <- rdist(s)
  max.dist <- max(D)
  sr <- 88.58573 # We estimated this previously and fix it to improve computation time
  ## Create a Sequence for %Spatial
  pct.spatial <- seq(0,.95,length=25)
  
  ## log-likelihoods
  ll <- matrix(0,nrow=1,ncol=length(pct.spatial))
  s2.hats <- ll
  beta.hats <- matrix(0,nrow=ncol(X),ncol=length(pct.spatial))
  the.it <- 0
  for(pct in pct.spatial){
    the.it <- the.it +1
    R <- pct*Matern(D,nu=nu,alpha=sr)+(1-pct)*diag(nrow(s))
    R.chol <- t(chol(R))
    first.piece <- forwardsolve(R.chol,X)
    XpRinvX <- t(first.piece)%*%first.piece
    last.piece <- forwardsolve(R.chol,y)
    beta.hat <- solve(XpRinvX)%*%t(first.piece)%*%last.piece
    beta.hats[,the.it] <- beta.hat
    ss <- forwardsolve(R.chol,y-X%*%beta.hat)
    ss <- t(ss)%*%ss
    s2.hat <- ss/nrow(s)
    s2.hats[1,pct==pct.spatial] <- s2.hat
    log.like <- -0.5*nrow(s)*log(s2.hat)-sum(log(diag(R.chol)))-0.5*ss/(s2.hat)
    ll[1,pct==pct.spatial] <- log.like
  }
  
  ## Find maximum likelihood estimate
  the.mle <- which(ll==max(ll),arr.ind=TRUE)
  pct <- pct.spatial[the.mle[1,2]]
  the.mle <- which(ll==max(ll))
  s2.hat <- s2.hats[the.mle]
  beta.hat <- beta.hats[,the.mle]
  
  ## Return List
  return(list(phi=sr,sigma2=pct*s2.hat,tau2=(1-pct)*s2.hat,beta.hat=beta.hat))
	
}

	
