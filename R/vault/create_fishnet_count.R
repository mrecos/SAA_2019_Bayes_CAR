library("sf")
library("tidyverse")
library("mapview")
library("FedData")
library("raster")
library("velox")

pass <- st_read("./SITE_DATA/PASS_w_info_EPSG4326.shp")

pass1 <- pass %>% 
  filter(test == 1) %>% 
  st_transform(crs = 26918) %>% # UTM 18N meters
  mutate(area_m2 = as.numeric(st_area(.)))


pass_fishnet <- st_make_grid(st_as_sfc(st_bbox(st_buffer(pass1, 500))), 
                             cellsize = 2000, square = FALSE) %>%
  st_sf() %>% 
  mutate(fishnet_id = row_number(),
         f_area_m2 = as.numeric(st_area(.)))

elevation <- raster("C:/R_local/SAA_2019_Bayes_CAR/RAW/NED/13/USGS_NED_13_n41w076_ArcGrid/grdn41w076_13")

# mapview(pass_fishnet, alpha.regions = 0) + mapview(pass1, color = "red", col.regions = "red")

pass_intersect <- st_intersection(pass_fishnet, pass1)

pass2 <- pass_intersect %>% 
  mutate(intersect_area_m2 = as.numeric(st_area(.)),
         eligible = str_detect(NATIONAL_R, "SHPO: Eligible|Keeper: Eligible|Listed|NHL")) %>% 
  group_by(fishnet_id) %>% 
  summarise(count = n(),
            eligible = if_else(sum(eligible) > 0, 1, 0),
            pcnt_area = sum(intersect_area_m2)/mean(f_area_m2)) %>% 
  dplyr::select(fishnet_id, count, eligible,pcnt_area)

# mapview(pass_fishnet, alpha.regions = 0) + mapview(pass2, color = "red", col.regions = "red")

pass_fishnet2 <- pass_fishnet %>% 
  left_join(., st_drop_geometry(pass2), by = "fishnet_id") %>% 
  mutate(count = if_else(is.na(count), 0, as.numeric(count)),
         presence = if_else(count > 0, 1, 0),
         eligible = if_else(is.na(eligible), 0, eligible),
         pcnt_area = if_else(is.na(pcnt_area), 0, pcnt_area))

mapview(pass_fishnet2, zcol = "presence")

########### EXTRACT RASTER
vx <- velox(elevation) # cast raster to velox
pass_fishnet2$mean_elev <- vx$extract(ll(pass_fishnet2), fun = mean, small = TRUE)
pass_fishnet2$med_elev  <- vx$extract(ll(pass_fishnet2), fun = median, small = TRUE)
pass_fishnet2$max_elev  <- vx$extract(ll(pass_fishnet2), fun = max, small = TRUE)
pass_fishnet2$min_elev  <- vx$extract(ll(pass_fishnet2), fun = min, small = TRUE)

mapview(pass_fishnet2, zcol = "mean_elev")
