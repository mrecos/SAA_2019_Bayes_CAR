library("NLMR")
library("raster")
library("landscapetools")
library("velox")
library("sf")
library("tidyverse")


### Number of rows and columns in prediction rasters
## needed for making simulated rasters, as well as for predicting real-world rasters
cols = 100
rows = 100
fisnet_cell_size = 10

### RANDOM ENV and RANDOM CULTURAL ------------------------------------
# raster
RenvRcult1 <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = 1)
vx_RenvRcult1 <- velox(RenvRcult1) # cast raster to velox
RenvRcult2 <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = 1)
vx_RenvRcult2 <- velox(RenvRcult2) # cast raster to velox
RenvRcult3 <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = 1)
vx_RenvRcult3 <- velox(RenvRcult3) # cast raster to velox
RenvRcult4 <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = 1)
vx_RenvRcult4 <- velox(RenvRcult4) # cast raster to velox
# table
RenvRcult_fishnet <- st_make_grid(st_as_sfc(st_bbox(RenvRcult1)), 
                             cellsize = fisnet_cell_size, square = FALSE) %>%
  st_sf() %>% 
  mutate(fishnet_id = row_number(),
         count = rpois(nrow(.),8),
         mean_val1 = as.numeric(vx_RenvRcult1$extract(., fun = mean, small = TRUE)),
         mean_val2 = as.numeric(vx_RenvRcult2$extract(., fun = mean, small = TRUE)),
         mean_val3 = as.numeric(vx_RenvRcult3$extract(., fun = mean, small = TRUE)),
         mean_val4 = as.numeric(vx_RenvRcult4$extract(., fun = mean, small = TRUE)))
hex_plot <- gather(st_drop_geometry(RenvRcult_fishnet), id, val, -fishnet_id, -count) %>% 
  left_join(., RenvRcult_fishnet, by = "fishnet_id") %>% 
  st_sf() %>% 
  dplyr::select(id, val)
show_landscape(list(
  "mean_val1" = RenvRcult1,
  "mean_val2" = RenvRcult2,
  "mean_val3" = RenvRcult3,
  "mean_val4" = RenvRcult4), unique_scales = TRUE) +
  geom_sf(data = hex_plot, inherit.aes = FALSE, color = "gray30", alpha = 0)

pairs(as.matrix(st_drop_geometry(RenvRcult_fishnet[,3:6])))


# or is this semi-structured? structured may be totally homogeneous
### STRUCTURED ENV and RANDOM CULTURAL ------------------------------------
# raster
SenvRcult_grad <- NLMR::nlm_distancegradient(cols,rows, origin = c(20, 30, 10, 15))
SenvRcult_grad <- abs(1-(SenvRcult_grad+0.01))*2 # invert and weight
SenvRcult1 <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = cols*2) * SenvRcult_grad
SenvRcult1 <- raster::scale(SenvRcult1)
vx_SenvRcult1 <- velox(SenvRcult1) # cast raster to velox
SenvRcult2 <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = cols*2) * -SenvRcult_grad
SenvRcult2 <- raster::scale(SenvRcult2)
vx_SenvRcult2 <- velox(SenvRcult2) # cast raster to velox
SenvRcult3 <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = cols*2) * SenvRcult_grad
SenvRcult3 <- raster::scale(SenvRcult3)
vx_SenvRcult3 <- velox(SenvRcult3) # cast raster to velox
SenvRcult4 <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = cols*2) * -SenvRcult_grad
SenvRcult4 <- raster::scale(SenvRcult4)
vx_SenvRcult4 <- velox(SenvRcult4) # cast raster to velox
# table
SenvRcult_fishnet <- st_make_grid(st_as_sfc(st_bbox(SenvRcult1)), 
                                  cellsize = fisnet_cell_size, square = FALSE) %>%
  st_sf() %>% 
  mutate(fishnet_id = row_number(),
         count = rpois(nrow(.),8),
         mean_val1 = as.numeric(vx_SenvRcult1$extract(., fun = mean, small = TRUE)),
         mean_val2 = as.numeric(vx_SenvRcult2$extract(., fun = mean, small = TRUE)),
         mean_val3 = as.numeric(vx_SenvRcult3$extract(., fun = mean, small = TRUE)),
         mean_val4 = as.numeric(vx_SenvRcult4$extract(., fun = mean, small = TRUE)))

hex_plot <- gather(st_drop_geometry(SenvRcult_fishnet), id, val, -fishnet_id, -count) %>% 
  left_join(., SenvRcult_fishnet, by = "fishnet_id") %>% 
  st_sf() %>% 
  dplyr::select(id, val)

show_landscape(list(
  "mean_val1" = SenvRcult1,
  "mean_val2" = SenvRcult2,
  "mean_val3" = SenvRcult3,
  "mean_val4" = SenvRcult4), unique_scales = TRUE) +
  geom_sf(data = hex_plot, inherit.aes = FALSE, color = "gray30", alpha = 0)

pairs(as.matrix(st_drop_geometry(SenvRcult_fishnet[,3:6])))






s_var1 <- rescale_sim_raster(s_var1r, 50, 10) 
s_var2 <- rescale_sim_raster(s_var1r, 3, 2) 
b_var1r <- NLMR::nlm_gaussianfield(cols,rows,autocorr_range = 20)
b_var1 <- rescale_sim_raster(b_var1r, 100, 20) 
b_var2 <- rescale_sim_raster(b_var1r, 6, 3) 
### Create a site-present trend surface  (sim data only)
trend_coords <- sim_trend(cols, rows, n = 3)
coords <- trend_coords$coords
trend <- trend_coords$trend
inv_trend <- abs(1-trend)
var1 <- (s_var1 * trend) + (b_var1 * inv_trend)
var2 <- (s_var2 * trend) + (b_var2 * inv_trend)
#### end simulated data creation ####

### Create raster stack of predictor variables
pred_var_stack <- raster::stack(var1, var2)
names(pred_var_stack) <- c("var1","var2")
### scale rasters to training data
pred_var_stack_scaled <- scale_prediction_rasters(pred_var_stack, params, verbose = 0)
### Predict raster (single chunk, not in parallel) 
pred_rast <- KLR_raster_predict(pred_var_stack_scaled, ngb = ngb, params, split = FALSE, ppside = NULL,
                                progress = FALSE, parallel = FALSE)
### plot with simulated sites
rasterVis::levelplot(pred_rast, margin = FALSE, par.settings=viridisTheme()) +
  layer(sp.points(sp.points(SpatialPoints(coords), pch=15, cex = 2.25, col = "red")), columns=1)
