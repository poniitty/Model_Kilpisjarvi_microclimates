library(terra)
library(sf)
library(tidyverse)

########################################################################
# Create DEM

latest <- st_read("/appl/data/geo/mml/dem2m/2008_latest/dem2m.shp")

aoi <- st_read("data/coordinates/aoi.gpkg") %>% 
  st_transform(crs = st_crs(latest))

latest_t <- latest[aoi,]

# Merge files

f <- latest_t$path

rast.list <- list()
for(ii in 1:length(f)) { rast.list[ii] <- rast(f[ii]) }

rast.list <- sprc(rast.list)
rast.mosaic <- mosaic(rast.list)

rast.mosaic <- terra::crop(rast.mosaic, aoi)

# resample to the SCD resolution
scd <- rast("data/predictors/scd_whole.tif")
rast.mosaic <- project(rast.mosaic, scd)

plot(rast.mosaic)
plot(scd)

writeRaster(round(rast.mosaic*100), paste0("data/predictors/dem.tif"),
            datatype = "INT4U", overwrite = T)

unlink(list.files(tempdir(), full.names = T, recursive = T))

###############################################################
# Topo indices

library(Rsagacmd, lib.loc = "/projappl/project_2003061/Rpackages/")

saga_version("saga_cmd")
saga <- saga_gis(cores = future::availableCores(), raster_backend = "terra")

dem <- rast("data/predictors/demm.tif")

dem <- saga$ta_preprocessor$fill_sinks_xxl_wang_liu(elev = dem)
writeRaster(dem, "data/predictors/dem_wang.tif", overwrite = TRUE)
swi <- saga$ta_hydrology$saga_wetness_index(dem = dem, suction = 2, slope_type = 0)
plot(swi$twi)
writeRaster(round(swi$twi*100), "data/predictors/swi_suct2.tif", datatype = "INT2S", overwrite = T)
swi <- saga$ta_hydrology$saga_wetness_index(dem = dem, suction = 8, slope_type = 0)
plot(swi$twi)
writeRaster(round(swi$twi*100), "data/predictors/swi_suct8.tif", datatype = "INT2S", overwrite = T)
swi <- saga$ta_hydrology$saga_wetness_index(dem = dem, suction = 64, slope_type = 0)
plot(swi$twi)
writeRaster(round(swi$twi*100), "data/predictors/swi_suct64.tif", datatype = "INT2S", overwrite = T)
swi <- saga$ta_hydrology$saga_wetness_index(dem = dem, suction = 256, slope_type = 0)
plot(swi$twi)
writeRaster(round(swi$twi*100), "data/predictors/swi_suct256.tif", datatype = "INT2S", overwrite = T)

# Radiation

dem <- rast("data/predictors/dem.tif")/100

svi <- saga$ta_lighting$sky_view_factor(dem = dem, radius = 1000)
plot(svi$svf)

writeRaster(svi$svf, "data/predictors/svi.tif", overwrite = T)


svi <- rast("data/predictors/svi.tif")
for(month in 1:12){
  
  print(month)
  rad <- saga$ta_lighting$potential_incoming_solar_radiation(grd_dem = dem,
                                                             grd_svf = svi,
                                                             period = 1, 
                                                             location = 1,
                                                             day = paste0("2020-",month,"-15"),
                                                             hour_step = 1)
  
  writeRaster(round(rad$grd_total*1000), paste0("data/predictors/pisr_",month,".tif"), datatype = "INT2U", overwrite = T)
  saga_remove_tmpfiles()
}

# TPI

dem <- rast("data/predictors/dem.tif")/100

tpi <- saga$ta_morphometry$topographic_position_index_tpi(dem = dem, radius_max = 15)
plot(tpi)
writeRaster(round(tpi*100), "data/predictors/tpi15.tif", datatype = "INT2S", overwrite = T)

tpi <- saga$ta_morphometry$topographic_position_index_tpi(dem = dem, radius_max = 30)
plot(tpi)
writeRaster(round(tpi*100), "data/predictors/tpi30.tif", datatype = "INT2S", overwrite = T)

tpi <- saga$ta_morphometry$topographic_position_index_tpi(dem = dem, radius_max = 90)
plot(tpi)
writeRaster(round(tpi*100), "data/predictors/tpi90.tif", datatype = "INT2S", overwrite = T)

dem <- aggregate(dem, 5, mean, na.rm = T)

tpi <- saga$ta_morphometry$topographic_position_index_tpi(dem = dem, radius_max = 250)
plot(tpi)
writeRaster(round(tpi*100), "data/predictors/tpi250.tif", datatype = "INT2S", overwrite = T)

tpi <- saga$ta_morphometry$topographic_position_index_tpi(dem = dem, radius_max = 500)
plot(tpi)
writeRaster(round(tpi*100), "data/predictors/tpi500.tif", datatype = "INT2S", overwrite = T)

tpi <- saga$ta_morphometry$topographic_position_index_tpi(dem = dem, radius_max = 1000)
plot(tpi)
writeRaster(round(tpi*100), "data/predictors/tpi1000.tif", datatype = "INT2S", overwrite = T)


# SLOPE
slp <- saga$ta_morphometry$slope_aspect_curvature(elevation = dem)
slp$slope <- 180/pi*slp$slope
writeRaster(round(slp$slope*100), "data/predictors/slope.tif", 
            overwrite = T, datatype = "INT2U")

# Wind indices
wind <- saga$ta_morphometry$wind_exposition_index(dem = dem, maxdist = 0.5,step = 30)
writeRaster(round(wind*1000), "data/predictors/windexp.tif", 
            overwrite = T, datatype = "INT2U")
plot(rast("data/predictors/windexp.tif"))

wind <- saga$ta_morphometry$wind_exposition_index(dem = dem, maxdist = 0.1,step = 30)
writeRaster(round(wind*1000), "data/predictors/windexp100.tif", 
            overwrite = T, datatype = "INT2U")
plot(rast("data/predictors/windexp100.tif"))

# CostDistance to water
MASTER <- rast("data/predictors/dem.tif")/100

v1 <- st_read("data/vectors/waters.gpkg")
v2 <- st_read("data/vectors/rivers.gpkg")

v1$value <- 1
v2$value <- 1

r1 <- rasterize(as(v1, "SpatVector"), MASTER, field = "value")
r2 <- rasterize(as(v2, "SpatVector"), MASTER, field = "value")

d1 <- saga$grid_tools$proximity_grid(features = sum(c(r1,r2), na.rm = T))
wd <- rast(d1$distance)
plot(wd)
wd[is.na(MASTER)] <- NA
names(wd) <- "disttowater"
writeRaster(round(wd), "data/predictors/disttowater.tif", 
            overwrite = T, datatype = "INT2U")

# Water effect
slo <- saga$ta_morphometry$slope_aspect_curvature(elevation = MASTER, unit_slope = 1)

dc1 <- saga$grid_analysis$accumulated_cost(dest_grid = sum(c(r1,r2), na.rm = T), cost = slo$slope, dest_type = 1)
wd <- rast(dc1$accumulated)

wd[wd > 100] <- 100
wd <- abs(wd-100)
plot(wd)
writeRaster(round(wd), "data/predictors/watereffect.tif", 
            overwrite = T, datatype = "INT2U")

saga_remove_tmpfiles()

# Snow effect
# install.packages("whitebox", repos="http://R-Forge.R-project.org", lib = "/projappl/project_2003061/Rpackages/")
library(whitebox, lib.loc = "/projappl/project_2003061/Rpackages/")
whitebox::wbt_init()
whitebox::wbt_version()

dem <- rast("data/predictors/dem.tif")/100

invdem <- saga$grid_tools$invert_grid(dem)
writeRaster(invdem, "/scratch/project_2007415/temp/invdem.tif", 
            overwrite = T)

snow <- rast("data/predictors/scd_whole.tif")

weeks <- seq(from = 150, to = 258, by = 7)

snow3 <- rast()
for(i in weeks){
  print(i)
  
  snow2 <- ifel(snow >= i, 1, 0)
  snow2[snow2 < 1] <- NA
  writeRaster(snow2, "/scratch/project_2007415/temp/scd.tif", 
              overwrite = T, datatype = "INT1U")
  
  wbt_downslope_distance_to_stream(dem = "/scratch/project_2007415/temp/invdem.tif", 
                                   streams = "/scratch/project_2007415/temp/scd.tif", 
                                   output = "/scratch/project_2007415/temp/temp1.tif")
  r <- rast("/scratch/project_2007415/temp/temp1.tif")
  
  r[r > 100] <- 100
  r[is.na(r)] <- 100
  r[r == 0] <- NA
  
  snow3 <- c(snow3, r)
}
gc()

invsnow <- saga$grid_tools$invert_grid(sum(snow3, na.rm = T))
invsnow[is.na(dem)] <- NA
plot(invsnow)
writeRaster(invsnow, "data/predictors/snoweffect.tif", 
            overwrite = T)

# TWI

library(whitebox, lib.loc = "/projappl/project_2003061/Rpackages/")
whitebox::wbt_init()

r <- rast("data/predictors/slope.tif")
r[is.na(r)] <- 1040
plot(r)
writeRaster(r, "data/predictors/slope_degree.tif", overwrite = T)

wbt_d8_flow_accumulation(input = "data/predictors/dem_wang.tif",
                        output = "data/predictors/d8_sca.tif", 
                        out_type = "specific contributing area")
wbt_d_inf_flow_accumulation(input = "data/predictors/dem_wang.tif",
                         output = "data/predictors/Dinf_sca.tif", 
                         out_type = "Specific Contributing Area")
wbt_fd8_flow_accumulation(dem = "data/predictors/dem_wang.tif",
                            output = "data/predictors/FD8f_sca.tif", 
                            out_type = "specific contributing area")

wbt_wetness_index(sca = "data/predictors/d8_sca.tif",
                  slope = "data/predictors/slope_degree.tif",
                  output = "data/predictors/TWI_d8.tif")
wbt_wetness_index(sca = "data/predictors/Dinf_sca.tif",
                  slope = "data/predictors/slope_degree.tif",
                  output = "data/predictors/TWI_Dinf.tif")
wbt_wetness_index(sca = "data/predictors/FD8f_sca.tif",
                  slope = "data/predictors/slope_degree.tif",
                  output = "data/predictors/TWI_FD8f.tif")
# Post process

twi <- rast("data/predictors/TWI_d8.tif")
plot(twi)
twi <- round(twi*100)
writeRaster(twi, "data/predictors/TWI_d8.tif", overwrite = T, datatype = "INT2S")
twi <- rast("data/predictors/TWI_Dinf.tif")
plot(twi)
twi <- round(twi*100)
writeRaster(twi, "data/predictors/TWI_Dinf.tif", overwrite = T, datatype = "INT2S")
twi <- rast("data/predictors/TWI_FD8f.tif")
plot(twi)
twi <- round(twi*100)
writeRaster(twi, "data/predictors/TWI_FD8f.tif", overwrite = T, datatype = "INT2S")
