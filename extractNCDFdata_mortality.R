library(ncdf)
library(doParallel)
load("./RData/Spatial911PtPtrn.RData")

## Takes a Date object and returns vector with year, month, day
dateToYMD <- function(date) {
  full <- strsplit(as.character(date), "-")
  year <- full[[1]][1]
  month <- full[[1]][2]
  day <- full[[1]][3]
  c(year, month, day)
}

## Gets the values from an NCDF file for a given variable on a given day
## Converts it to a vector so that latitude and longitude coordinates are
## consistent with hrldas.grid. Note that this function also centers the 
## temperature data. If uncentered data is desired, this code needs to be
## changed.
getVarOnDay <- function(ncdf.file, varid, day) {
  var <- get.var.ncdf(ncdf.file, varid=varid)[,,as.numeric(day)]
  var.vec <- as.vector(t(var))
  var.vec
}

# Get unique dates where a 911 call occurred
unique.dates <- as.Date(unique(paste(death$Year, death$Month, death$Day)), "%Y %m %d")

# Get dates for three day lags
lag1 <- unique.dates - 1
lag2 <- unique.dates - 2
lag3 <- unique.dates - 3

filepath <- "../Data/Heat/RawNCDFFiles/houston_hrldas_"
endpiece <- "_alldays.nc"

list.names <- as.character(unique.dates)
vars <- c("HI_MAX", "HI_MIN", "T2MIN", "T2MAX", "SW_MIN", "SW_MAX")

# Register cores to use
numCores <- 16
registerDoParallel(cores=numCores)
GetTempData <- function() {
  foreach (i=1:length(unique.dates)) %dopar% {
    ######################################
    # Get variables for the observed day #
    ######################################
    ymd <- dateToYMD(unique.dates[i])
    ncdf.file <- open.ncdf(paste(filepath, ymd[1],  ymd[2],  endpiece, sep=""))
    var0 <- sapply(1:length(vars), function(x) { getVarOnDay(ncdf.file, vars[x], ymd[3]) })

    #######################
    # Get lag 1 variables #
    #######################
    ymd.l1 <- dateToYMD(lag1[[i]])
    # If the month is different, open up a new ncdf file
    if (ymd.l1[2] != ymd[2]) {
      close.ncdf(ncdf.file)
      ncdf.file <- open.ncdf(paste(filepath, ymd.l1[1], ymd.l1[2], endpiece, sep=""))
    }
    var1 <- sapply(1:length(vars), function(x) { getVarOnDay(ncdf.file, vars[x], ymd.l1[3]) })
    
    #######################
    # Get lag 2 variables #
    #######################
    ymd.l2 <- dateToYMD(lag2[[i]])
    # If the month for lag2 is different than the month for lag1, open up a new ncdf file
    if (ymd.l2[2] != ymd.l1[2]) {
      close.ncdf(ncdf.file) 
      ncdf.file <- open.ncdf(paste(filepath, ymd.l2[1], ymd.l2[2], endpiece, sep=""))
    }
    var2 <- sapply(1:length(vars), function(x) { getVarOnDay(ncdf.file, vars[x], ymd.l2[3]) })
  
    #######################
    # Get lag 3 variables #
    #######################
    ymd.l3 <- dateToYMD(lag3[[i]])
    # Get new ncdf if month has changed
    if (ymd.l3[2] != ymd.l2[2]) {
      close.ncdf(ncdf.file)
      ncdf.file <- open.ncdf(paste(filepath, ymd.l3[1], ymd.l3[2], endpiece, sep=""))
    }
    var3 <- sapply(1:length(vars), function(x) { getVarOnDay(ncdf.file, vars[x], ymd.l3[3]) })
  
    # Combine variables and insert as data into list
    close.ncdf(ncdf.file)
    df <- cbind(hrldas.grid, var0, var1, var2, var3)[kp.gp,]
    df[,-c(1:2)] <- apply(df[,-c(1:2)], 2, scale, center=TRUE, scale=FALSE)
    names(df) <- c("LAT", "LON", paste(vars, ".0", sep=""), paste(vars,".1",sep=""),
      paste(vars,".2",sep=""), paste(vars,".3",sep=""))
    df
  }
}
temp.data <- GetTempData()

names(temp.data) <- list.names
save(temp.data, file="./RData/MortalityTempData_withMissing.RData")
