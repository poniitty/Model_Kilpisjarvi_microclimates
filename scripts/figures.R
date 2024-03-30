library(tidyverse)
library(lubridate)
library(zoo)
library(patchwork)
library(terra)
library(tidyterra, lib.loc = "/projappl/project_2003061/Rpackages/")
library(viridis)
library(colorspace)
library(parallel)
library(sf)

setwd("/projappl/project_2007415/repos/Model_Kilpisjarvi_microclimates/")

# Map of the study sites

library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")

r <- rast("data/predictors/dem.tif")/100

w <- st_read("data/vectors/waters.gpkg") %>% 
  mutate(area = as.numeric(st_area(.))) %>% 
  filter(area > 100000) %>% 
  st_crop(r)

p <- st_read("data/coordinates/combined_coords.gpkg") %>% 
  st_transform(crs = st_crs(r)) %>% 
  mutate(Study_design = ifelse(startsWith(site, "MAL"), 1, NA),
         Study_design = ifelse(startsWith(site, "AIL"), 1, Study_design),
         Study_design = ifelse(startsWith(site, "MI"), 2, Study_design),
         Study_design = ifelse(startsWith(site, "RA"), 3, Study_design),
         Study_design = ifelse(startsWith(site, "Haap"), 4, Study_design),
         Study_design = ifelse(startsWith(site, "K"), 5, Study_design)) %>% 
  mutate(Study_design = factor(Study_design)) %>% 
  st_crop(st_bbox(r))

ggplot() +
  geom_spatraster(data = r, maxcell = 500000) +
  scale_fill_gradientn(colours = grey.colors(100),
                       na.value = "white") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_light() +
  geom_sf(data = p, aes(color = Study_design, group = Study_design), size = 1) +
  scale_color_viridis_d(option = "B") +
  geom_sf(data = w, color = "lightblue", fill = "lightblue")


ggplot(data = world) +
  geom_sf() +
  geom_sf(data = p %>% st_bbox %>% st_as_sfc %>% st_centroid(), color = "red", size = 4) +
  theme_bw() +
  coord_sf(ylim = c(6000000, 8000000), xlim = c(-500000,1200000), 
           crs = st_crs(p))


#######################################################

d <- bind_rows(read_csv("../MICROCLIMATES/output/KIL/tomst_data_cleaned.csv") %>% 
                 rename(probl = error_tomst) %>% 
                 select(-tomst_id),
               read_csv("data/tomst_data/Oulu_tomst_data.csv"),
               read_csv("data/tomst_data/Roslin_tomst_data.csv")) %>% 
  mutate(datetime = with_tz(datetime, tzone = "Etc/GMT-2")) %>% 
  filter(minute(datetime) %in% c(0,30,60)) %>% 
  mutate(date = as_date(datetime)) %>% 
  relocate(date, .after = datetime) %>% 
  filter(date >= "2019-09-01") %>% 
  filter(date <= "2023-08-31") %>% 
  filter(!startsWith(site, "L"),
         !startsWith(site, "AE"),
         !startsWith(site, "SE"),
         !startsWith(site, "SAL"),
         !startsWith(site, "X"))

# Calibration function for the moisture count values for unknown soils from Kopecky et al. 2020
cal_funNA <- function(x) {((-1.34e-8) * (x^2) + (2.50e-4) * x + (-1.58e-1))*100 }

# Calibrate the moisture values
d %>% mutate(moist = round(cal_funNA(moist),1)) %>% 
  mutate(moist = ifelse(T1 < 1 | probl == 1, NA, moist)) -> d

d <- d %>% 
  relocate(moist, .after = T3)

d <- d %>% select(-T2)

# Fill missing observations (max length of three consecutive NAs) with running average
# Linear interpolation
d %>% group_by(site) %>% 
  mutate(across(T1:moist, ~na.approx(.x, datetime, maxgap = 3, na.rm = F))) -> d

########################################################################
# AGGREGATE TO DAILY VALUES

notNA_prop <- function(x){ round(sum(is.finite(x)/length(x))*100,1) }
gc()

system.time({
  daily <- mclapply(unique(d$site), function(x){
    print(x)
    suppressWarnings({
      d %>% 
        filter(site == x) %>% 
        group_by(site, date) %>%
        summarise(across(T1:moist, ~notNA_prop(.x), 
                         .names = "{.col}_prop"),
                  across(T1:T3, list(mean = mean), 
                         na.rm = T, .names = "{.col}_{.fn}"),
                  moist_mean = mean(moist, na.rm = T),
                  probl = max(probl, na.rm = T),
                  .groups = "drop") -> daily
    })
    
    return(daily)
  }, mc.cores = future::availableCores()) %>% bind_rows()
})
gc()


# T1

T1 <- daily %>% 
  select(site, date, probl, T1_prop, T1_mean) %>% 
  filter(T1_prop >= 90) %>% 
  filter(probl %in% c(0,3,5:8,10)) %>% 
  pivot_wider(id_cols = date, names_from = site, values_from = T1_mean) %>% 
  arrange(date)

T1 <- T1 %>% 
  mutate(across(where(is.numeric), ~!is.na(.x))) %>% 
  rowwise() %>% 
  mutate(T1 = mean(c_across(where(is.logical)))) %>% 
  select(date, T1)

# T3

T3 <- daily %>% 
  select(site, date, probl, T3_prop, T3_mean) %>% 
  filter(T3_prop >= 90) %>% 
  filter(probl %in% c(0,6,9,11)) %>% 
  pivot_wider(id_cols = date, names_from = site, values_from = T3_mean) %>% 
  arrange(date)

T3 <- T3 %>% 
  mutate(across(where(is.numeric), ~!is.na(.x))) %>% 
  rowwise() %>% 
  mutate(T3 = mean(c_across(where(is.logical)))) %>% 
  select(date, T3)

# Moist

moist <- daily %>% 
  select(site, date, probl, moist_prop, moist_mean) %>% 
  filter(moist_prop >= 90) %>% 
  filter(probl %in% c(0,3:5,7,11)) %>% 
  filter(month(date) %in% 6:8) %>% 
  pivot_wider(id_cols = date, names_from = site, values_from = moist_mean) %>% 
  arrange(date)

moist <- moist %>% 
  mutate(across(where(is.numeric), ~!is.na(.x))) %>% 
  rowwise() %>% 
  mutate(moist = mean(c_across(where(is.logical)))) %>% 
  select(date, moist)

#

all <- full_join(T1, T3) %>% 
  full_join(., moist) %>% 
  pivot_longer(cols = T1:moist, names_to = "var", values_to = "Proportion")

all %>% 
  mutate(var = ifelse(var == "moist", "Moisture", var)) %>% 
  ggplot(aes(y = Proportion, x = date)) +
  geom_bar(stat = "identity", width = 1, color = "gray30") +
  facet_grid(rows = vars(var)) +
  theme_bw() +
  ylab("Proportion of sites with data") + xlab("Date")


###############################################################
#

d <- read_csv("data/imputed_tomst_data.csv")

gg1 <- d %>% 
  group_by(date) %>% 
  summarise(Tmin = min(T1_min, na.rm = T),
            Tmean = mean(T1_mean, na.rm = T),
            Tmax = max(T1_max, na.rm = T)) %>% 
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymax=Tmax, ymin=Tmin), fill="#291208",alpha=.6) +
  geom_line(aes(y = Tmean)) +
  theme_bw() + ylab("Soil temperature") + xlab("Date") +
  ylim(-36, 42)

gg2 <- d %>% 
  group_by(date) %>% 
  summarise(Tmin = min(T3_min, na.rm = T),
            Tmean = mean(T3_mean, na.rm = T),
            Tmax = max(T3_max, na.rm = T)) %>% 
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymax=Tmax, ymin=Tmin), fill="#227988",alpha=.6) +
  geom_line(aes(y = Tmean)) +
  theme_bw() + ylab("Near-surface air temperature") + xlab("Date") +
  ylim(-36, 42)

gg3 <- d %>% 
  group_by(date) %>% 
  summarise(Tmin = min(moist, na.rm = T),
            Tmean = mean(moist, na.rm = T),
            Tmax = max(moist, na.rm = T)) %>% 
  mutate(across(Tmin:Tmax, ~ifelse(is.infinite(.x), NA, .x))) %>% 
  mutate(across(Tmin:Tmax, ~ifelse(.x < 0, 0, .x))) %>% 
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymax=Tmax, ymin=Tmin), fill="#002090",alpha=.7) +
  geom_line(aes(y = Tmean)) +
  theme_bw() + ylab("Soil moisture") + xlab("Date") +
  ylim(0, 62)

gg1 / gg2 / gg3
ggsave("visuals/time_series.pdf", width = 13, height = 7, units = "in")


###########################################################################
# Maps

list.files("output/predictions/")

r <- rast("output/predictions/T3_GDD3.tif")
gg1 <- ggplot() +
  geom_spatraster(data = r, maxcell = 1000000) +
  scale_fill_gradientn(colours = diverge_hsv(n = 100, v = 0.7, s = 1.3), na.value = "white",
                       limits = c(681, 1310)) +
  ggtitle("T3_GDD3") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_light()
r2 <- crop(r, ext(492000, 495000, 7660000, 7663000))
gg12 <- ggplot() +
  geom_spatraster(data = r2, maxcell = 100000) +
  scale_fill_gradientn(colours = diverge_hsv(n = 100, v = 0.7, s = 1.3), na.value = "white",
                       limits = c(681, 1310), guide = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void()

r <- rast("output/predictions/T3_max.tif")
gg2 <- ggplot() +
  geom_spatraster(data = r, maxcell = 1000000) +
  scale_fill_gradientn(colours = heat_hcl(n = 100,  h = c(0, 90), c = c(80, 30),
                                          l = c(30, 90), power = c(1/5, 1.5), rev = T), 
                       na.value = "white",
                       limits = c(20.28, 34.90)) +
  ggtitle("T3_max") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_light()
r2 <- crop(r, ext(492000, 495000, 7660000, 7663000))
gg22 <- ggplot() +
  geom_spatraster(data = r2, maxcell = 100000) +
  scale_fill_gradientn(colours = heat_hcl(n = 100,  h = c(0, 90), c = c(80, 30),
                                          l = c(30, 90), power = c(1/5, 1.5), rev = T), 
                       na.value = "white",
                       limits = c(20.28, 34.90), guide = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void()


r <- rast("output/predictions/T1_FDD.tif")
gg3 <- ggplot() +
  geom_spatraster(data = r, maxcell = 1000000) +
  scale_fill_gradientn(colours = heat_hcl(n = 100,  h = c(340, 243), c = c(24, 63),
                                          l = c(77, 43), power = 4, rev = T), 
                       na.value = "white",
                       limits = c(-1156.82, -38.79)) +
  ggtitle("T1_FDD") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_light()
r2 <- crop(r, ext(492000, 495000, 7660000, 7663000))
gg32 <- ggplot() +
  geom_spatraster(data = r2, maxcell = 100000) +
  scale_fill_gradientn(colours = heat_hcl(n = 100,  h = c(340, 243), c = c(24, 63),
                                          l = c(77, 43), power = 4, rev = T), 
                       na.value = "white",
                       limits = c(-1156.82, -38.79), guide = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void()


r <- rast("output/predictions/moist_mean_summer.tif")
gg4 <- ggplot() +
  geom_spatraster(data = r, maxcell = 1000000) +
  scale_fill_gradientn(colours = diverge_hcl(n = 100,  h = c(230, 50), c = 150,
                                             l = c(70, 40), power = 1.5, rev = T),
                       na.value = "white",
                       limits = c(12, 52.77)) +
  ggtitle("moist_mean_summer") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_light()
r2 <- crop(r, ext(492000, 495000, 7660000, 7663000))
gg42 <- ggplot() +
  geom_spatraster(data = r2, maxcell = 100000) +
  scale_fill_gradientn(colours = diverge_hcl(n = 100,  h = c(230, 50), c = 150,
                                             l = c(70, 40), power = 1.5, rev = T),
                       na.value = "white",
                       limits = c(12, 52.77), guide = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void()


gg <- (gg1  + gg2) / (gg3 + gg4)
gg
ggsave("visuals/maps.pdf", device = "pdf", height = 12, width = 12)
gg2 <- (gg12  + gg22) / (gg32 + gg42)
gg2
ggsave("visuals/zoom_maps.pdf", device = "pdf", height = 8, width = 8)
