death <- read.csv("../data/deathData.csv", header=T)
### Eliminate unused months
death <- death[which(4 < death$Month & death$Month < 10),]
load("./RData/Spatial911PtPtrn.RData")

### Delete unused variables
death$p99_09_ID <- NULL
death$p99_05_ID <- NULL
death$CID <- NULL
death$D_DDATT <- NULL
death$DeathAcc <- NULL
death$D_MARRY <- NULL
death$Marital <- NULL
death$D_ALC_US <- NULL
death$Alcohol <- NULL
death$D_GEOCOD <- NULL
death$D_PLACTY <- NULL
death$DeathPlace <- NULL
death$D_PLOTH <- NULL
death$D_OCSTRE <- NULL
death$D_OCITNM <- NULL
death$D_OCZIP <- NULL
death$D_OCZIPEX <- NULL
death$IFINJURY <- NULL
death$Injury <- NULL


### Note that D_UNDCAU lists the ICD-10 codes for cause of death
### We only use those prefixed with:
###   I - Diseases of the circulatory system
###   J - Diseases of the respiratory system
###   N - Diseases of the genitourinary system
###   R - Symptoms and signs involving ciruclatory and respiratory systems
### Since these are the causes of death that have been linked to extreme heat.
death <- death[which(substr(death$D_UNDCAU, 1, 1) == "I" |
                     substr(death$D_UNDCAU, 1, 1) == "J" |
                     substr(death$D_UNDCAU, 1, 1) == "N" |
                     substr(death$D_UNDCAU, 1, 1) == "R"), ]


death.spat <- SpatialPoints(cbind(death$D_LONG, death$D_LAT))
microbenchmark(cont.ind <- ifelse(rowSums(gContains(houston, death.spat, byid=TRUE)) > 0, TRUE, FALSE), times=10)
death <- death[cont.ind, ]

write.csv(death, "./RData/cleanDeathData.csv", row.names=FALSE)
