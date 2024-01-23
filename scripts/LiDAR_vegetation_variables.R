# install.packages("lidR", lib = "/projappl/project_2003061/Rpackages/")
library(lidR)
library(sf)
library(raster)
library(tidyverse)

get_lidr_threads()
set_lidr_threads(future::availableCores())
get_lidr_threads()

max_z <- 25 # in meters

latest <- st_read("/appl/data/geo/mml/laserkeilaus/2008_latest/2008_latest.shp") %>% 
  mutate(path = unlist(lapply(path, function(x)paste0("/",tail(str_split(x, "//")[[1]],1)))))

# Prepocess study area polygons
rois <- st_read("data/coordinates/aoi.gpkg") %>% 
                    st_transform(crs = st_crs(latest)) %>% 
  mutate(full_name = "KIL")

# rois <- rois %>% filter(area != "KIL") %>% filter(area != "RAS")

for(area_name in rois$full_name){
  # area_name <- "KIL"
  roi <- rois %>%
    filter(full_name == area_name) %>%
    st_transform(crs = st_crs(latest))
  
  pols <- st_make_grid(roi, cellsize = 3000) %>% 
    st_as_sf() %>% 
    mutate(id = 1:nrow(.))
  
  for(i in pols$id){
    # i <- 6
    print(area_name)
    print(i)
    
    if(!file.exists(paste0("/projappl/project_2007415/temp/chm_",i,".tif"))){
      roi_t <- pols %>% filter(id == i) %>% 
        st_buffer(10)
      
      latest_t <- latest[roi_t,]
      
      if(nrow(latest_t) > 0){
        lass <- readLAS(latest_t$path)
        
        # CLIP THE LAS
        lass <- clip_roi(lass, roi_t)
        if(!is.null(lass)){
          if(lass@header@PHB$`Number of point records` > 1000){
            
            dtm <- grid_terrain(lass, 1, tin())
            plot(dtm)
            
            lass <- normalize_height(lass, tin(), na.rm = TRUE)
            
            lass <- filter_poi(lass, Z <= max_z)
            
            chm <- grid_canopy(lass, res = 1, algorithm = dsmtin(max_edge = 5))
            
            chm[chm < 0] <- 0
            chm[is.na(chm)] <- 0
            
            roi_t <- pols %>% filter(id == i)
            
            chm <- crop(chm, roi_t, snap = "out")
            
            writeRaster(round(chm*100), paste0("/projappl/project_2007415/temp/chm_",i,".tif"),
                        format = "GTiff", datatype = "INT2U", overwrite = T)
          }
        }
      }
    }
  }
  
  #################################################################
  # Merge files
  
  vars <- unique(unlist(lapply(list.files("/projappl/project_2007415/temp/", pattern = "tif$"), function(x) str_split(x, "_")[[1]][1])))
  
  for(i in vars){
    print(i)
    
    f <- list.files("/projappl/project_2007415/temp", pattern = paste0(i,"_"), full.names = T)
    
    rast.list <- list()
    for(ii in 1:length(f)) { rast.list[ii] <- raster(f[ii]) }
    rast.list$fun <- mean
    rast.list$tolerance <- 0.5
    rast.mosaic <- do.call(mosaic,rast.list)
    
    plot(rast.mosaic, main = i)
    
    writeRaster(rast.mosaic, paste0("data/predictors/",i,"_",area_name,".tif"),
                format = "GTiff", datatype = dataType(raster(f[1])), overwrite = T)
    
  }
  
  unlink(list.files("/projappl/project_2007415/temp", full.names = T))
  
}


###########################################################
# Finetune the vegetation height with slope and NDVI

library(terra)

chm <- rast("data/predictors/chm.tif")
ndvi <- rast("data/predictors/NDVI_max.tif")
slp <- rast("data/predictors/slope.tif")

chm <- project(chm, ndvi)
slp <- project(slp, ndvi)

ndvi <- mask(ndvi, slp)
chm <- mask(chm, ndvi)

plot(ndvi)
plot(chm)
plot(slp)

chm[slp > 6000] <- 0
chm[ndvi < 400] <- 0
chm[slp > 3000 & ndvi < 600] <- 0


ccv5m <- ifel(chm > 500, 1, 0)
fm <- focalMat(ccv5m, 6, "circle")
fm[fm > 0] <- 1
ccv5m <- round(focal(ccv5m, fm, "mean", na.policy = "omit", na.rm = T)*100)
plot(ccv5m)

ccv2m <- ifel(chm > 200, 1, 0)
fm <- focalMat(ccv2m, 6, "circle")
fm[fm > 0] <- 1
ccv2m <- round(focal(ccv2m, fm, "mean", na.policy = "omit", na.rm = T)*100)
plot(ccv2m)

shrubs <- ifel(chm >= 25 & chm <= 100, 1, 0)
fm <- focalMat(shrubs, 6, "circle")
fm[fm > 0] <- 1
shrubs <- round(focal(shrubs, fm, "mean", na.policy = "omit", na.rm = T)*100)
shrubs[slp > 3000] <- 0
plot(shrubs)

writeRaster(round(chm), "data/predictors/chm.tif", overwrite = T, datatype = "INT2U")
writeRaster(ccv5m, "data/predictors/ccv5m.tif", overwrite = T, datatype = "INT1U")
writeRaster(ccv2m, "data/predictors/ccv2m.tif", overwrite = T, datatype = "INT1U")
writeRaster(shrubs, "data/predictors/shrubs.tif", overwrite = T, datatype = "INT1U")
