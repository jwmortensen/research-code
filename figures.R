load("./RData/Spatial911PtPtrn.RData")
load("./RData/MHDrawsHI_MAX_2015_05_19.RData")
load("./RData/TempDataNoMiss.RData")

PlotOutput <- function(draws) {
  # Plot lambda.star values
  plot(draws$lambda.star[, 1], type="l")
  plot(draws$lambda.star[, 948], type="l")
  plot(draws$lambda.star[, 707], type="l")
  plot(draws$lambda.star[, 1428], type="l")
  
  # plot beta values
  plot(draws$beta[, 1], type="l")
  plot(draws$beta[, 2], type="l")
  plot(draws$beta[, 3], type="l")
  plot(draws$beta[, 4], type="l")
}

PlotOutput(draws.new)

PostMeans <- function(draws, burnin) {
  post.delta <- list(mn=mean(draws$delta), 
                     lower.95=quantile(draws$delta, 0.025), 
                     upper.95=quantile(draws$delta, 0.975))
  post.lambda <- list(mn=apply(draws$lambda.star[-c(1:burnin), ], 2, mean),
                      lower.95=apply(draws$lambda.star[-c(1:burnin), ], 2, quantile, probs=0.025),
                      upper.95=apply(draws$lambda.star[-c(1:burnin), ], 2, quantile, probs=0.975))
  post.beta <- list(mn=apply(draws$beta[-c(1:burnin), ], 2, mean), 
                    lower.95=apply(draws$beta[-c(1:burnin), ], 2, quantile, probs=0.025), 
                    upper.95=apply(draws$beta[-c(1:burnin), ], 2, quantile, probs=0.975))
  post.race <- list(mn=apply(draws$race, 2, mean),
                    lower.95=apply(draws$race, 2, quantile, probs=0.025),
                    upper.95=apply(draws$race, 2, quantile, probs=0.975))
  post.gender <- list(mn=apply(draws$gender, 2, mean),
                      lower.95=apply(draws$gender, 2, quantile, probs=0.025),
                      upper.95=apply(draws$gender, 2, quantile, probs=0.975))
  post.age <- list(mn=apply(draws$age, 2, mean),
                   lower.95=apply(draws$age, 2, quantile, probs=0.025),
                   upper.95=apply(draws$age, 2, quantile, probs=0.975))
  list(delta=post.delta, lambda.star=post.lambda, beta=post.beta, race=post.race, gender=post.gender, age=post.age)
}
post.draws <- PostMeans(draws.mortality, 50)
CalcLambda <- function(lambda.star, beta, temps) {
  exp(lambda.star + temps%*%beta) / sum(exp(lambda.star + temps%*%beta))
}
hi_max.2006.08.20 <- as.matrix(temp.data.nomiss[["2006-08-20"]][which(names(temp.data.nomiss[["2006-08-20"]]) %in% c("HI_MAX.0", "HI_MAX.1", "HI_MAX.2", "HI_MAX.3"))])
lambda.w.beta <- CalcLambda(post.draws$lambda.star$mn, post.draws$beta$mn, hi_max.2006.08.20)
lambda.no.beta <- CalcLambda(post.draws$lambda.star$mn, c(0,0,0,0), hi_max.2006.08.20)

pdf("mortality_results1.pdf")
y.grid <- matrix(hrldas.grid[,1],nrow=125)
x.grid <- matrix(hrldas.grid[,2],nrow=125)
plot.grid <- matrix(NA, nrow=125, ncol=125)
plot.grid[kp.gp] <- lambda.no.beta
image.plot(x.grid, y.grid, plot.grid, axes=FALSE, frame.plot=TRUE, ann=FALSE,
           xlim=c(-95.79, -95.0), ylim=c(29.5, 30.15)) #, zlim=c(0, 0.012)) 
dev.off()
pdf("mortality_results2.pdf")
plot.grid2 <- matrix(NA, nrow=125, ncol=125)
plot.grid2[kp.gp] <- lambda.w.beta
image.plot(x.grid, y.grid, plot.grid2, axes=FALSE, frame.plot=TRUE, ann=FALSE, 
           xlim=c(-95.79, -95.0), ylim=c(29.5, 30.15)) #, zlim=c(0, 0.012))
dev.off()

head(death)
spat.death <- SpatialPointsDataFrame(cbind(death$D_LONG, death$D_LAT), death)
houston.death <- gIntersection(houston, spat.death, byid=TRUE)
