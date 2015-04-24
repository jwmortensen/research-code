library(LatticeKrig)

load("../data/Spatial911PtPtrn.RData")
load("../data/TempDataNoMiss.RData")
load("../data/TemperatureData_withMissing.RData")

# Create matrix to plot values
y.grid <- matrix(hrldas.grid[,1],nrow=125)
x.grid <- matrix(hrldas.grid[,2],nrow=125)
plot.grid <- matrix(NA, nrow=125, ncol=125)

numToCheck <- 2
for (i in seq(1, length(temp.data.nomiss), length=numToCheck)) {
  for (j in 3:ncol(temp.data.nomiss[[i]])) {
    plot.grid[kp.gp] <- temp.data[[i]][,j]
    pdf(paste(i,"compare",colnames(temp.data[[i]])[j],".pdf", sep=""))
    image.plot(x.grid, y.grid, plot.grid, axes=FALSE, frame.plot=TRUE, 
               main=paste("miss", colnames(temp.data[[i]])[j]),
               xlim=c(-95.79, -95.0), ylim=c(29.5, 30.15))
    plot.grid[kp.gp] <- temp.data.nomiss[[i]][,j]
    image.plot(x.grid, y.grid, plot.grid, axes=FALSE, frame.plot=TRUE, 
               main=paste("nomiss", colnames(temp.data.nomiss[[i]])[j]),
               xlim=c(-95.79, -95.0), ylim=c(29.5, 30.15))
    dev.off()
  }
}

