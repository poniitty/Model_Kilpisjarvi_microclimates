library(tidyverse)
library(sf)
library(readxl)

d <- read_csv("data/microclimate_variables.csv")

p1 <- st_read("data/coordinates/all_tomsts_2022.gpkg") %>% 
  rename(site = id) %>% 
  select(site) %>% 
  mutate(site = toupper(site)) %>% 
  mutate(site = gsub("_","",site)) %>% 
  mutate(site = ifelse(site == "SAA19", "SAA019", site))
p2 <- st_read("data/coordinates/roslin_plots.gpkg") %>% 
  st_transform(crs = 32634) %>% 
  select(site) %>% 
  mutate(site = toupper(site))
p3 <- read_excel("/scratch/project_2007415/microclim/Oulu/Haaps_TOMST_Kilp_Pekalle/Kilpisjarvi_TOMST-plots_Pekka.xlsx") %>% 
  filter(complete.cases(.)) %>% 
  st_as_sf(coords = c("long WGS-84", "lat WGS-84"), crs = 4326) %>% 
  select(ID) %>% 
  rename(site = ID,
         geom = geometry) %>% 
  st_transform(crs = 32634)

p <- bind_rows(p1, p2, p3)

d$site[!d$site %in% p$site]
p$site[!p$site %in% d$site]

p <- p %>% 
  filter(site %in% d$site)

st_write(p, "data/coordinates/combined_coords.gpkg", append = FALSE)
