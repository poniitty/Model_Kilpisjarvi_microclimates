library(sf)
library(tidyverse)

aoi <- st_read("data/coordinates/aoi.gpkg")

list.files("/appl/data/geo/mml/maastotietokanta/2023/gpkg", full.names = T)
st_layers("/appl/data/geo/mml/maastotietokanta/2023/gpkg/MTK-vakavesi_23-03-02.gpkg")
st_layers("/appl/data/geo/mml/maastotietokanta/2023/gpkg/MTK-virtavesi_23-03-02.gpkg")

p <- bind_rows(read_sf("/appl/data/geo/mml/maastotietokanta/2023/gpkg/MTK-vakavesi_23-03-02.gpkg", "jarvi"),
               read_sf("/appl/data/geo/mml/maastotietokanta/2023/gpkg/MTK-virtavesi_23-03-02.gpkg", "virtavesialue"))

waters <- st_crop(p, aoi %>% st_transform(crs = st_crs(p)))
waters <- waters %>% st_transform(crs = st_crs(aoi))
plot(st_geometry(waters))

st_write(waters, "data/vectors/waters.gpkg")

# Smaller rivers

p <- read_sf("/appl/data/geo/mml/maastotietokanta/2023/gpkg/MTK-virtavesi_23-03-02.gpkg", "virtavesikapea")

rivers <- st_crop(p, aoi %>% st_transform(crs = st_crs(p)))
rivers <- rivers %>% st_transform(crs = st_crs(aoi))
plot(st_geometry(rivers))

st_write(rivers, "data/vectors/rivers.gpkg")


# Roads

st_layers("/appl/data/geo/mml/maastotietokanta/2023/gpkg/MTK-tie_23-03-02.gpkg")

p <- read_sf("/appl/data/geo/mml/maastotietokanta/2023/gpkg/MTK-tie_23-03-02.gpkg", "tieviiva")

roads <- st_crop(p, aoi %>% st_transform(crs = st_crs(p)))
roads <- roads %>% st_transform(crs = st_crs(aoi))
plot(st_geometry(roads))

st_write(roads, "data/vectors/roads.gpkg")

# buildings

st_layers("/appl/data/geo/mml/maastotietokanta/2023/gpkg/MTK-rakennus_23-03-02.gpkg")

p <- read_sf("/appl/data/geo/mml/maastotietokanta/2023/gpkg/MTK-rakennus_23-03-02.gpkg", "rakennus")

buildings <- st_crop(p, aoi %>% st_transform(crs = st_crs(p)))
buildings <- buildings %>% st_transform(crs = st_crs(aoi))
plot(st_geometry(buildings))

st_write(buildings, "data/vectors/buildings.gpkg")
