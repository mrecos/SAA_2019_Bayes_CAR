library("raster")
library("velox")
library("sf")
library(mapview)

ll <- function(dat, crs = 4326){
  st_transform(dat, crs)
}


shp_loc <- "D:/PROJECTS_RESEARCH/PASS_DistReg_data/SHP"
raster_loc <- "D:/PROJECTS_RESEARCH/PASS_DistReg_data/variables/"

sites <- read_sf(file.path(shp_loc,"pass_pnts","pass_pnts_test_area.shp")) %>% 
  st_transform(., crs = 4326)

var <- "cd_h7"
up <- raster(file.path(raster_loc,"subarea", "upland_section_6", paste0(var, ".grd")))
up <- raster::projectRaster(up, crs = st_crs(sites)$proj4string)
rv <- raster(file.path(raster_loc,"subarea", "riverine_section_6", paste0(var, ".grd")))
rv <- raster::projectRaster(rv, crs = st_crs(sites)$proj4string)

# sites_bbox <- st_bbox(sites)
sites_buffer <- ll(st_buffer(sites %>% st_transform(.,crs=26918),1000))

up <- crop(up, extent(sites_buffer))
rv <- crop(rv, extent(sites_buffer))

raster::writeRaster(up, filename = file.path(raster_loc,"tiff",paste0("up_",var,".tiff")))
raster::writeRaster(rv, filename = file.path(raster_loc,"tiff",paste0("rv_",var,".tiff")))

# 
# mapview(rv6_ed_h2) + mapview(sites)
# 
# 
# # xx <- raster::extract(elevation, as(ll(pass_fishnet), "Spatial"), small = TRUE) # TOO SLOW!
# 
# vx <- velox(elevation) # cast raster to velox
# pass_fishnet2$mean_elev <- vx$extract(ll(pass_fishnet2), fun = mean, small = TRUE)
# pass_fishnet2$med_elev  <- vx$extract(ll(pass_fishnet2), fun = median, small = TRUE)
# pass_fishnet2$max_elev  <- vx$extract(ll(pass_fishnet2), fun = max, small = TRUE)
# pass_fishnet2$min_elev  <- vx$extract(ll(pass_fishnet2), fun = min, small = TRUE)