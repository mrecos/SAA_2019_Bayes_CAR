library("raster")
library("velox")
library("sf")
library(mapview)
library(tidyverse)

ll <- function(dat, crs = 4326){
  st_transform(dat, crs)
}
UTM18 <- function(dat, crs = 32618){
  st_transform(dat, crs)
}
scale_this <- function(x){
  as.vector(scale(x))
}

fishnet_cell_size <- 500

shp_loc <- "D:/PROJECTS_RESEARCH/PASS_DistReg_data/SHP"
raster_loc <- "D:/PROJECTS_RESEARCH/PASS_DistReg_data/variables"

sites <- read_sf(file.path(shp_loc,"pass_pnts","PASS_poly_AEAC.shp")) %>% 
  st_transform(., crs = 4326)

all_ed_h2 <- raster(file.path(raster_loc,"tiff", "all_ed_h2.tif"))
all_cd_h7 <- raster(file.path(raster_loc,"tiff", "all_cd_h7.tif"))
all_std_32c <- raster(file.path(raster_loc,"tiff", "all_std_32c.tif"))
all_elev_2_drainh <- raster(file.path(raster_loc,"tiff", "all_elev_2_drainh.tif"))

env_vars <- raster::stack(all_ed_h2, all_cd_h7, all_elev_2_drainh, all_std_32c)

# cast to Velox
all_ed_h2_vx   <- velox(all_ed_h2) # cast raster to velox
all_cd_h7_vx   <- velox(all_cd_h7) # cast raster to velox
all_std_32c_vx <- velox(all_std_32c) # cast raster to velox
all_elev_2_drainh_vx <- velox(all_elev_2_drainh) # cast raster to velox

# make fishnet and join site attributes
sites_fishnet <- st_make_grid(st_as_sfc(st_bbox(st_buffer(UTM18(sites),1000))),
                                  cellsize = fishnet_cell_size, square = FALSE) %>%
  st_sf() %>%
  dplyr::mutate(fishnet_id = dplyr::row_number(),
                f_area_m2  = as.numeric(st_area(.))) %>%
  ll(.)

sites_intersect <- st_intersection(sites_fishnet, sites)

sites2 <- sites_intersect %>% 
  mutate(intersect_area_m2 = as.numeric(st_area(.)),
         eligible = str_detect(NATIONAL_R, "SHPO: Eligible|Keeper: Eligible|Listed|NHL")) %>% 
  group_by(fishnet_id) %>% 
  summarise(count = n(),
            eligible = if_else(sum(eligible) > 0, 1, 0),
            pcnt_area = sum(intersect_area_m2)/mean(f_area_m2)) %>% 
  dplyr::select(fishnet_id, count, eligible,pcnt_area)

# mapview(sites_fishnet, alpha.regions = 0) + mapview(sites2, color = "red", col.regions = "red")

sites_fishnet2 <- sites_fishnet %>% 
  left_join(., st_drop_geometry(sites2), by = "fishnet_id") %>% 
  mutate(count = if_else(is.na(count), 0, as.numeric(count)),
         presence = if_else(count > 0, 1, 0),
         eligible = if_else(is.na(eligible), 0, eligible),
         pcnt_area = if_else(is.na(pcnt_area), 0, pcnt_area))

# mapview(sites_fishnet2, zcol = "presence")

## Join raster Data to fishnet w/ sites info
### Make features
sites_fishnet3 <- sites_fishnet2 %>% 
  mutate(mean_ed_h2   = as.numeric(all_ed_h2_vx$extract(., fun = mean, small = TRUE)),
         mean_cd_h7   = as.numeric(all_cd_h7_vx$extract(., fun = mean, small = TRUE)),
         mean_std_32c = as.numeric(all_std_32c_vx$extract(., fun = mean, small = TRUE)),
         mean_elev_2_drainh = as.numeric(all_elev_2_drainh_vx$extract(., fun = mean, small = TRUE)))
pairs(as.matrix(st_drop_geometry(sites_fishnet3[,c(3,6,7:10)])))

## Prepare vars/colnames for Stan
sites_fishnet_stan <- sites_fishnet3 %>% 
  mutate(mean_ed_h2 = scale_this(mean_ed_h2),
         mean_cd_h7 = scale_this(mean_cd_h7),
         mean_std_32c = scale_this(mean_std_32c),
         mean_elev_2_drainh = scale_this(mean_elev_2_drainh)) %>% 
  dplyr::select(count, mean_val1 = mean_ed_h2, mean_val2 = mean_cd_h7, 
                       mean_val3 = mean_std_32c, mean_val4 = mean_elev_2_drainh)
