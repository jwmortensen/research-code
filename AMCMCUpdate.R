AMCMC.update <- function(draw,cur.mn,cur.var,cur.it){
	if(cur.it >0){
		mn <- ((cur.it-1)*cur.mn+draw)/cur.it
		if(cur.it==1){
			v <- matrix(0,nrow=length(draw),ncol=length(draw))
		} else {
			v <- (cur.it-2)*cur.var+(cur.it-1)*(cur.mn%*%t(cur.mn))+draw%*%t(draw)
			v <- (v-cur.it*(mn%*%t(mn)))/(cur.it-1)
		}
	} else {
		mn <- matrix(0,nrow=nrow(cur.mn),ncol=1)
		v <- matrix(0,nrow=nrow(draw),ncol=nrow(draw))
	}
	return(list(mn=mn,var=v))
}

## Test run of the function
# mu <- 0
# the.var <- 1
# x <- rnorm(100,0,1)
# mcmc.draws <- matrix(0,nrow=10000,ncol=2)
# amcmc <- list(mn=matrix(0,nrow=2,ncol=1),var=matrix(0,nrow=2,ncol=2))

# for(mcmc.it in 1:nrow(mcmc.draws)){
	# prop.var <- diag(rep(0.01^2,2))+(mcmc.it>250)*amcmc$var
	# prop <- c(mu,log(the.var))+t(chol(prop.var))%*%rnorm(2)
	# alph <- sum(dnorm(x,mean=prop[1],sd=sqrt(exp(prop[2])),log=TRUE)-dnorm(x,mean=mu,sd=sqrt(the.var),log=TRUE))
	# if(log(runif(1))<alph){
		# mu <- prop[1]
		# the.var <- exp(prop[2])
	# }
	# amcmc <- AMCMC.update(matrix(c(mu,log(the.var)),ncol=1),amcmc$mn,amcmc$var,mcmc.it)
	
	# mcmc.draws[mcmc.it,] <- c(mu,the.var)
# }
	