library(tidyverse)
library(sf)
library(terra)

tifs <- list.files("data/predictors/", pattern = ".tif$")

tifs <- tifs[!grepl("_sca", tifs)]
tifs <- tifs[!grepl("dem_", tifs)]
tifs <- tifs[!grepl("demm", tifs)]
tifs <- tifs[!grepl("predstack", tifs)]
tifs <- tifs[!grepl("slope_degree", tifs)]
tifs <- tifs[!grepl("svi.tif", tifs)]

p <- st_read("data/coordinates/combined_coords.gpkg")

ext_rinf <- function(x, p2, buff, dir){
  
  r <- rast(paste0(dir, x))
  names(r) <- gsub(".tif","",x)
  p3 <- p2 %>% 
    st_transform(crs = crs(r, proj = T)) %>% 
    st_buffer(buff)
  e1 <- extract(r, vect(p3), fun = mean, exact = TRUE, na.rm = T) %>% round
  
  df <- p2 %>% 
    st_drop_geometry() %>% 
    bind_cols(., e1) %>% 
    select(-ID)
  
  return(df)
}

all <- lapply(tifs, function(x){ext_rinf(x, p2 = p, buff = 1, dir = "data/predictors/")})

all <- all %>% 
  reduce(full_join, by = "site")
summary(all)

all %>% filter(!complete.cases(.))
all %>% select(-site) %>% cor(., use = "pairwise.complete.obs") %>% round(., 2)

write_csv(all, "data/predictors/env_data.csv")

# With 20-m bigger buffer

tifs <- c(paste0("pisr_", 1:12, ".tif"))

all <- lapply(tifs, function(x){ext_rinf(x, p2 = p, buff = 20, dir = "data/predictors/")})

all <- all %>% 
  reduce(full_join, by = "site")
summary(all)

names(all)[-1] <- unlist(lapply(names(all)[-1], function(x){
  paste0(str_split(x,"_")[[1]][1], "20",
         ifelse(grep("_",x), paste0("_", str_split(x,"_")[[1]][2]), ""))
}))

all %>% select(-site) %>% cor(., use = "pairwise.complete.obs") %>% round(., 2)

write_csv(all, "data/predictors/env_buf20m_data.csv")

# Planet PCA

pca <- rast("/scratch/project_2007415/Planet/KilpisjarviAll/extra/planet_PCA.tif")

names(pca) <- gsub("PC","PCA",names(pca))

pca <- extract(pca, p %>% 
                st_transform(crs = crs(pca, proj = T)) %>% 
                st_buffer(2) %>% vect,
              fun = mean, exact = TRUE, na.rm = T) %>% round

pca <- p %>% 
  st_drop_geometry() %>% 
  bind_cols(., pca) %>% 
  select(-ID)

pca %>% select(-site) %>% cor(., use = "pairwise.complete.obs") %>% round(., 2)

write_csv(pca, "data/predictors/planet_pca_data.csv")
