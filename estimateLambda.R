library(parallel)
load("./RData/TempDataNoMiss.RData")

lambda.star <- rep(1, 1428)
beta <- rep(1,4)
vars <- c("HI_MAX","HI_MIN","T2MAX","T2MIN","SW_MIN","SW_MAX")
postfix <- c(".0",".1",".2",".3")
lagged.vars <- paste(vars[1], postfix, sep="")

calc.log.lambda <- function(H, lambda.star, beta, vars) {
  lambda.star + as.matrix(H[,vars])%*%beta - log(sum(exp(lambda.star + as.matrix(H[,vars])%*%beta)))
}

response <- mclapply(temp.data.nomiss, calc.log.lambda, lambda.star=lambda.star, beta=beta, vars=lagged.vars, mc.cores=16)
