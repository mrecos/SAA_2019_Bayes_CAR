library("sf")
library("tidyverse")
library("mapview")
library("FedData")
library("raster")
library("velox")
library("fitdistrplus")

pass <- st_read("./SITE_DATA/PASS_w_info_EPSG4326.shp")
phys <- st_read("./SITE_DATA/Physiographic_Sections_of_Pennsylvania.shp")

# phys_R3 <- phys %>% 
#   filter(SECTION_ %in% c("High Plateau Section", "Northwestern Glaciated Plateau Section")) %>% 
#   st_transform(crs = 26918)

phys_R3 <- phys %>% 
  mutate(dissolve = 1) %>% 
  group_by(dissolve) %>% 
  summarise() %>% 
  st_transform(crs = 26918)

# pass_R3 <- pass %>%
#   filter(Region == 2) %>%
#   st_transform(crs = 26918) %>% # UTM 18N meters
#   mutate(area_m2 = as.numeric(st_area(.)))

pass_R3 <- pass %>%
  st_transform(crs = 26918) %>% # UTM 18N meters
  mutate(area_m2 = as.numeric(st_area(.)))

mapview(phys_R3) + mapview(pass_R3)

pass_fishnet_X <- st_make_grid(st_as_sfc(st_bbox(st_buffer(phys_R3, 50))), 
                             cellsize = 10000, square = FALSE) %>% 
  st_intersection(., phys_R3) %>%
  st_sf() %>% 
  mutate(fishnet_id = row_number())

pass_R3_int <-  st_intersection(pass_fishnet_X, pass_R3)

pass_R3_2 <- pass_R3_int %>% 
    group_by(fishnet_id) %>% 
    summarise(count = n()) %>% 
    dplyr::select(fishnet_id, count)

pass_fishnet_X2 <- pass_fishnet_X %>% 
  left_join(., st_drop_geometry(pass_R3_2), by = "fishnet_id") %>% 
  mutate(count = if_else(is.na(count), 0, as.numeric(count)),
         presence = if_else(count > 0, 1, 0))

# mapview(pass_fishnet_X2, zcol = "count") + mapview(pass_R3)

descdist(pass_fishnet_X2$count)










