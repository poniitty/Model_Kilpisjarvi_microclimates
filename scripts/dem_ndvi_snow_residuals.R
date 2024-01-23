library(terra)
library(sf)
library(tidyverse)

dem <- rast("data/predictors/dem.tif")
names(dem) <- "dem"
ndvi <- rast("data/predictors/NDVI_max.tif")

s <- c(dem, ndvi)

waters <- st_read("/projappl/project_2007415/repos/Model_Kilpisjarvi_microclimates/data/vectors/waters.gpkg")
s <- mask(s, vect(waters), inverse = TRUE) 
names(s) <- c("dem","ndvi")

d <- spatSample(s, 10000, "random", na.rm = T)

m <- lm(ndvi ~ dem, d)
summary(m)

pr <- predict(dem, m)
ndvires <- round(ndvi-pr)
names(ndvires) <- "ndviresi"
plot(ndvires, range = c(-300, 300))
writeRaster(ndvires, "data/predictors/ndviresi.tif", overwrite = T, datatype = "INT2S")

# Snow

snow <- rast("data/predictors/scd.tif")

s <- c(dem, snow)

waters <- st_read("/projappl/project_2007415/repos/Model_Kilpisjarvi_microclimates/data/vectors/waters.gpkg")
s <- mask(s, vect(waters), inverse = TRUE) 
names(s) <- c("dem","snow")

d <- spatSample(s, 100000, "random", na.rm = T)

m <- lm(snow ~ dem, d)
summary(m)

pr <- predict(dem, m)
snowres <- round(snow-pr)
names(snowres) <- "snowresi"
plot(snowres, range = c(-50, 50))
writeRaster(snowres, "data/predictors/snowresi.tif", overwrite = T, datatype = "INT2S")

rm(s, dem, pr, ndvi, snow, snowres, ndvires)

