library("FedData")
library("raster")
library("velox")
library("sf")

ll <- function(dat, crs = 4326){
  st_transform(dat, crs)
}

# elevation <- get_ned(as(ll(pass_fishnet), "Spatial"), label = "test", res = "13")

elevation <- raster("C:/R_local/SAA_2019_Bayes_CAR/RAW/NED/13/USGS_NED_13_n41w076_ArcGrid/grdn41w076_13")

# xx <- raster::extract(elevation, as(ll(pass_fishnet), "Spatial"), small = TRUE) # TOO SLOW!

vx <- velox(elevation) # cast raster to velox
pass_fishnet2$mean_elev <- vx$extract(ll(pass_fishnet2), fun = mean, small = TRUE)
pass_fishnet2$med_elev  <- vx$extract(ll(pass_fishnet2), fun = median, small = TRUE)
pass_fishnet2$max_elev  <- vx$extract(ll(pass_fishnet2), fun = max, small = TRUE)
pass_fishnet2$min_elev  <- vx$extract(ll(pass_fishnet2), fun = min, small = TRUE)