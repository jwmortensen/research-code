library(ncdf)
library(doParallel)
load("../Spatial911PtPtrn.RData")

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
## consistent with hrldas.grid
getVarOnDay <- function(ncdf.file, varid, day) {
  var <- get.var.ncdf(ncdf.file, varid=varid)[,,as.numeric(day)]
  var.vec <- as.vector(t(var))
  var.vec
}

# Get unique dates where a 911 call occurred
split.dates <- strsplit(as.character(strptime(calls$DOY, format="%j")), "-")
months <- numeric(length(split.dates))
days <- numeric(length(split.dates))
for (i in 1:length(split.dates)) {
  months[i] <- split.dates[[i]][2]
  days[i] <- split.dates[[i]][3]
}
unique.dates <- as.Date(unique(paste(calls$Year, months, days)), "%Y %m %d")

# Get dates for three day lags
lag1 <- unique.dates - 1
lag2 <- unique.dates - 2
lag3 <- unique.dates - 3

filepath <- "../Data/Heat/RawNCDFFiles/houston_hrldas_"
endpiece <- "_alldays.nc"

list.names <- as.character(unique.dates)
vars <- c("HI_MAX", "HI_MIN")

# Register cores to use
numCores <- 4
registerDoParallel(cores=numCores)
temp.data <- foreach (i=1:length(unique.days)) %dopar% {
  ######################################
  # Get variables for the observed day #
  ######################################
  ymd <- dateToYMD(unique.dates[i])
  ncdf.file <- open.ncdf(paste(filepath, ymd[1],  ymd[2],  endpiece, sep=""))
  var0 <- matrix(NA, ncol=length(vars), nrow=nrow(hrldas.grid))
  for (j in 1:length(vars)) {
    var0[,j] <-  getVarOnDay(ncdf.file, vars[j], ymd[3])
  }
  #######################
  # Get lag 1 variables #
  #######################
  ymd.l1 <- dateToYMD(lag1[[i]])
  # If the month is different, open up a new ncdf file
  if (ymd.l1[2] != ymd[2]) {
    close.ncdf(ncdf.file)
    ncdf.file <- open.ncdf(paste(filepath, ymd.l1[1], ymd.l1[2], endpiece, sep=""))
  }
  var1 <- matrix(NA, ncol=length(vars), nrow=nrow(hrldas.grid))
  for (k in 1:length(vars)) {
    var1[,k] <-  getVarOnDay(ncdf.file, vars[k], ymd[3])
  }
  
  #######################
  # Get lag 2 variables #
  #######################
  var2 <- matrix(NA, ncol=length(vars), nrow=nrow(hrldas.grid))
  ymd.l2 <- dateToYMD(lag2[[i]])

  # If the month for lag2 is different than the month for lag1, open up a new ncdf file
  if (ymd.l2[2] != ymd.l1[2]) {
    close.ncdf(ncdf.file) 
    ncdf.file <- open.ncdf(paste(filepath, ymd.l2[1], ymd.l2[2], endpiece, sep=""))
  }
  for (l in 1:length(vars)) {
    var2[,l] <-  getVarOnDay(ncdf.file, vars[l], ymd[3])
  }

  #######################
  # Get lag 3 variables #
  #######################
  var3 <- matrix(NA, ncol=length(vars), nrow=nrow(hrldas.grid))
  ymd.l3 <- dateToYMD(lag3[[i]])

  # Get new ncdf if month has changed
  if (ymd.l3[2] != ymd.l2[2]) {
    close.ncdf(ncdf.file)
    ncdf.file <- open.ncdf(paste(filepath, ymd.l3[1], ymd.l3[2], endpiece, sep=""))
  }
  for (m in 1:length(vars)) {
    var3[,m] <-  getVarOnDay(ncdf.file, vars[m], ymd[3])
  }

  # Combine variables and insert as data into list
  close.ncdf(ncdf.file)
  df <- cbind(hrldas.grid, var0, var1, var2, var3)[kp.gp,]
  names(df) <- c("LAT", "LON", paste(vars, ".0", sep=""), paste(vars,".1",sep=""), paste(vars,".2",sep=""), paste(vars,".3",sep=""))
  df
}

names(time.temp.data) <- list.names
save(temp.data, file="TemperatureData_withMissing.RData")
