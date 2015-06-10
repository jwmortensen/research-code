library(sp)
library(rgeos)
library(GISTools)
library(LatticeKrig)
load("./RData/Spatial911PtPtrn.RData")
load("./RData//AllTempData.RData")

FacToNum <- function(x) {
  as.numeric(as.character(x))
}

sp.pts <- SpatialPoints(cbind(pp.grid[,2], pp.grid[,1]))
sp.grd <- points2grid(sp.pts, tolerance=0.000763359)
grd.layer <- as.SpatialPolygons.GridTopology(sp.grd)
grd.cont <- ifelse(colSums(gContains(grd.layer, sp.pts, byid=T)) > 0, TRUE, FALSE)
grd.layer <- grd.layer[grd.cont]
grd.layer <- SpatialPolygonsDataFrame(grd.layer,
              data = data.frame(c(1:length(sp.pts))),
              match.ID = FALSE)
names(grd.layer) <- "ID"
row.names(grd.layer) <- paste("g", 1:length(sp.pts), sep="")

census <- readShapePoly("../data/Census2010/census2010.shp")
census$TotalPop <- FacToNum(census$TotalPop)
census$PCTover65 <- FacToNum(census$PCTover65)

cen.cont <- ifelse(colSums(gIntersects(census, grd.layer, byid=T)) > 0, TRUE, FALSE)
census.layer <- census[cen.cont, ]
plot(census.layer)

int.res <- gIntersection(grd.layer, census.layer, byid=T)
tmp <- strsplit(names(int.res), " ")

grid.id <- (sapply(tmp, "[[", 1))
census.id <- (sapply(tmp, "[[", 2))

int.areas <- gArea(int.res, byid=T)
census.areas <- gArea(census.layer, byid=T)
grid.areas <- gArea(grd.layer, byid=T)

# This method, where you use the proportion of the census tract area, only works for count data, not percentages
cens.index <- match(census.id, row.names(census.layer))
census.areas <- census.areas[cens.index]
census.prop <- zapsmall(int.areas/census.areas, 3)
pop <- zapsmall(census.layer$TotalPop[cens.index] * census.prop, 5)

# Calculate percentage over 65
grid.areas <- grid.areas[grid.id]
grid.prop <- zapsmall(int.areas/grid.areas, 3)
pct65 <- zapsmall(FacToNum(census.layer$PCTover65[cens.index]) * grid.prop, 1)
census.df <- data.frame(census.layer)
census.pop <- census.df$TotalPop[cens.index]
census.pct65 <- census.df$PCTover65[cens.index]

df <- data.frame(grid.id, census.id, pct65, census.pct65, pop, census.pop, census.prop, grid.prop)
int.layer.pct65 <- xtabs(df$pct65~df$grid.id)
int.layer.pop <- xtabs(df$pop~df$grid.id)

index <- as.numeric(gsub("g", "", names(int.layer.pct65)))
pct65 <- numeric(dim(data.frame(grd.layer))[1])
pct65[index] <- int.layer.pct65
pop <- numeric(dim(data.frame(grd.layer))[1])
pop[index] <- int.layer.pop
int.layer <- SpatialPolygonsDataFrame(grd.layer, data=data.frame(data.frame(grd.layer), pop, pct65, coordinates(grd.layer)), match.ID=FALSE)

# png("TotalPopulation.png", width=720, height=540)
# par(mfrow = c(1, 2))
# par(mar=c(0, 0, 0, 0))
# shades <- auto.shading(census.layer$TotalPop, n = 9, cols=brewer.pal(9, "Greens"))
# choropleth(census.layer, census.layer$TotalPop, shades)
# choro.legend(-95.7, 30.3, shades, title="Total Population", cex=0.8)
# shades <- auto.shading(int.layer$pop, n = 9, cols=brewer.pal(9, "Greens"))
# choropleth(int.layer, int.layer$pop, shades)
# choro.legend(-95.7, 30.275, shades, title="Total Population", cex=0.8)
# dev.off()

# png("PercentOver65.png", width=720, height=540)
# par(mfrow = c(1, 2))
# par(mar = c(0, 0, 0, 0))
# shades <- auto.shading(FacToNum(census.layer$PCTover65), n = 9, cols=brewer.pal(9, "Reds"))
# choropleth(census.layer, FacToNum(census.layer$PCTover65), shades)
# choro.legend(-95.7, 30.3, shades, title="Percent Over 65", cex=0.8)
# shades <- auto.shading(int.layer$pct65, n = 9, cols = brewer.pal(9, "Reds"))
# choropleth(int.layer, int.layer$pct65, shades)
# choro.legend(-95.7, 30.275, shades, title="Percent Over 65", cex=0.8)
# dev.off()

int.layer.df <- data.frame(int.layer)
names(int.layer.df) <- c("ID", "Population", "PercentOver65", "Longitude", "Latitude")
int.layer.df <- int.layer.df[order(int.layer.df$Longitude, int.layer.df$Latitude), ]

vars <- c("HI_MAX","HI_MIN","T2MAX","T2MIN","SW_MIN","SW_MAX")
getSpatialMeans <- function(temp.var) {
  H <- lapply(temp.data.nomiss, function(x) { as.matrix(x[, temp.var]) })
  spatial.means <- apply(sapply(H, function(x) { apply(x, 1, mean) }), 1, mean)
  spatial.means
}

spatial.means <- lapply(vars, getSpatialMeans)
names(spatial.means) <- vars

intercept.df <- cbind(pp.grid, int.layer.df$Population, int.layer.df$PercentOver65, spatial.means)
names(intercept.df) <- c("Latitude", "Longitude", "Population", "PercentOver65", vars)
save(intercept.df, file="./RData/InterceptData.RData")

