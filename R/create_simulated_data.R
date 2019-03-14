library("NLMR")
library("raster")
library("landscapetools")
library("velox")
library("sf")
library("tidyverse")


### Number of rows and columns in prediction rasters
## needed for making simulated rasters, as well as for predicting real-world rasters
cols = 200
rows = 200
fisnet_cell_size = 10

### RANDOM ENV and RANDOM CULTURAL ------------------------------------
### Both env and counts are uniform
# raster
RenvRcult1 <- NLMR::nlm_random(cols,rows)
RenvRcult1 <- raster::scale(RenvRcult1)
vx_RenvRcult1 <- velox(RenvRcult1) # cast raster to velox
RenvRcult2 <- NLMR::nlm_random(cols,rows)
RenvRcult2 <- raster::scale(RenvRcult2)
vx_RenvRcult2 <- velox(RenvRcult2) # cast raster to velox
RenvRcult3 <- NLMR::nlm_random(cols,rows)
RenvRcult3 <- raster::scale(RenvRcult3)
vx_RenvRcult3 <- velox(RenvRcult3) # cast raster to velox
RenvRcult4 <- NLMR::nlm_random(cols,rows)
RenvRcult4 <- raster::scale(RenvRcult4)
vx_RenvRcult4 <- velox(RenvRcult4) # cast raster to velox
# table
RenvRcult_fishnet <- st_make_grid(st_as_sfc(st_bbox(RenvRcult1)), 
                             cellsize = fisnet_cell_size, square = FALSE) %>%
  st_sf() %>% 
  mutate(fishnet_id = row_number(),
         # logit rho ~ -1.5
         count = rpois(nrow(.),8),                     # random by G.D.P
         # logit rho very negative
         # count = sample(0:20,nrow(.), replace = TRUE), # or uniform across expected counts?
         # logit rho ~ 0
         # count = sample(0:1,nrow(.), replace = TRUE),    # or uniform across (0,1) ints?
         # count = runif(nrow(.)),                       # or unform (0,1)? # NO!, b/c Poisson
         mean_val1 = as.numeric(vx_RenvRcult1$extract(., fun = mean, small = TRUE)),
         mean_val2 = as.numeric(vx_RenvRcult2$extract(., fun = mean, small = TRUE)),
         mean_val3 = as.numeric(vx_RenvRcult3$extract(., fun = mean, small = TRUE)),
         mean_val4 = as.numeric(vx_RenvRcult4$extract(., fun = mean, small = TRUE)))
# hex_plot <- gather(st_drop_geometry(RenvRcult_fishnet), id, val, -fishnet_id, -count) %>%
#   left_join(., RenvRcult_fishnet, by = "fishnet_id") %>%
#   st_sf() %>%
#   dplyr::select(id, val)
# show_landscape(list(
#   "mean_val1" = RenvRcult1,
#   "mean_val2" = RenvRcult2,
#   "mean_val3" = RenvRcult3,
#   "mean_val4" = RenvRcult4), unique_scales = TRUE) +
#   geom_sf(data = hex_plot, inherit.aes = FALSE, color = "gray30", alpha = 0)
# 
pairs(as.matrix(st_drop_geometry(RenvRcult_fishnet[,2:6])))


### STRUCTURED ENV and STRCTURED CULTURAL ------------------------------------
# raster
# started with this as "structure" == uniformity, may need to be about tight spatial correlation
# SenvScult1 <- NLMR::nlm_random(cols,rows) # random b/ it will be set to a value of 1
# SenvScult1[] <- 1

SenvScult1 <- NLMR::nlm_distancegradient(cols,rows, origin = c(90, 110, 90, 110))
SenvScult1 <- SenvScult1 + NLMR::nlm_random(cols,rows)
vx_SenvScult1 <- velox(SenvScult1) # cast raster to velox
SenvScult2 <- NLMR::nlm_distancegradient(cols,rows, origin = c(90, 110, 90, 110))
SenvScult2 <- SenvScult2 + NLMR::nlm_random(cols,rows)
vx_SenvScult2 <- velox(SenvScult2) # cast raster to velox
SenvScult3 <- NLMR::nlm_distancegradient(cols,rows, origin = c(90, 110, 90, 110))
SenvScult3 <- SenvScult3 + NLMR::nlm_random(cols,rows)
vx_SenvScult3 <- velox(SenvScult3) # cast raster to velox
SenvScult4 <- NLMR::nlm_distancegradient(cols,rows, origin = c(90, 110, 90, 110))
SenvScult4 <- SenvScult4 + NLMR::nlm_random(cols,rows)
vx_SenvScult4 <- velox(SenvScult4) # cast raster to velox


SenvScult_sites <- NLMR::nlm_distancegradient(cols,rows, origin = c(90, 110, 90, 110))
SenvScult_sites <- SenvScult_sites + NLMR::nlm_random(cols,rows)
vx_SenvScult_sites <- velox(SenvScult_sites) # cast raster to velox

# table
## Strucutred sites will follow the same gradient, but need int counts... I guess
SenvScult_fishnet <- st_make_grid(st_as_sfc(st_bbox(SenvScult1)), 
                                  cellsize = fisnet_cell_size, square = FALSE) %>%
  st_sf() %>% 
  mutate(fishnet_id = row_number(),
         count = ,
         mean_val1 = as.numeric(vx_SenvScult1$extract(., fun = mean, small = TRUE)),
         mean_val2 = as.numeric(vx_SenvScult2$extract(., fun = mean, small = TRUE)),
         mean_val3 = as.numeric(vx_SenvScult3$extract(., fun = mean, small = TRUE)),
         mean_val4 = as.numeric(vx_SenvScult4$extract(., fun = mean, small = TRUE)))

# previous method of uniformity
  # st_sf() %>% 
  # mutate(fishnet_id = row_number(),
  #        count = 1,
  #        mean_val1 = 1,
  #        mean_val2 = 1,
  #        mean_val3 = 1,
  #        mean_val4 = 1)

# hex_plot <- gather(st_drop_geometry(SenvScult_fishnet), id, val, -fishnet_id, -count) %>% 
#   left_join(., SenvScult_fishnet, by = "fishnet_id") %>% 
#   st_sf() %>% 
#   dplyr::select(id, val)
# 
# show_landscape(list(
#   "mean_val1" = SenvScult1,
#   "mean_val2" = SenvScult1,
#   "mean_val3" = SenvScult1,
#   "mean_val4" = SenvScult1), unique_scales = TRUE) +
#   geom_sf(data = hex_plot, inherit.aes = FALSE, color = "gray30", alpha = 0)
# 
# pairs(as.matrix(st_drop_geometry(SenvScult_fishnet[,2:6])))


### RADNOM ENV and STRCTURED CULTURAL ------------------------------------
# raster
RenvScult1 <- NLMR::nlm_random(cols,rows)
RenvScult1 <- raster::scale(RenvScult1)
vx_RenvScult1 <- velox(RenvScult1) # cast raster to velox
RenvScult2 <- NLMR::nlm_random(cols,rows)
RenvScult2 <- raster::scale(RenvScult2)
vx_RenvScult2 <- velox(RenvScult2) # cast raster to velox
RenvScult3 <- NLMR::nlm_random(cols,rows)
RenvScult3 <- raster::scale(RenvScult3)
vx_RenvScult3 <- velox(RenvScult3) # cast raster to velox
RenvScult4 <- NLMR::nlm_random(cols,rows)
RenvScult4 <- raster::scale(RenvScult4)
vx_RenvScult4 <- velox(RenvScult4) # cast raster to velox
# table
RenvScult_fishnet <- st_make_grid(st_as_sfc(st_bbox(RenvScult1)), 
                                  cellsize = fisnet_cell_size, square = FALSE) %>%
  st_sf() %>% 
  mutate(fishnet_id = row_number(),
         count = 1,
         mean_val1 = as.numeric(vx_RenvScult1$extract(., fun = mean, small = TRUE)),
         mean_val2 = as.numeric(vx_RenvScult2$extract(., fun = mean, small = TRUE)),
         mean_val3 = as.numeric(vx_RenvScult3$extract(., fun = mean, small = TRUE)),
         mean_val4 = as.numeric(vx_RenvScult4$extract(., fun = mean, small = TRUE)))
# hex_plot <- gather(st_drop_geometry(RenvScult_fishnet), id, val, -fishnet_id, -count) %>% 
#   left_join(., RenvScult_fishnet, by = "fishnet_id") %>% 
#   st_sf() %>% 
#   dplyr::select(id, val)
# show_landscape(list(
#   "mean_val1" = RenvScult1,
#   "mean_val2" = RenvScult2,
#   "mean_val3" = RenvScult3,
#   "mean_val4" = RenvScult4), unique_scales = TRUE) +
#   geom_sf(data = hex_plot, inherit.aes = FALSE, color = "gray30", alpha = 0)
# 
# pairs(as.matrix(st_drop_geometry(RenvScult_fishnet[,2:6])))



### STRUCTURED ENV and RANDOM CULTURAL ------------------------------------
# raster
SenvRcult1 <- NLMR::nlm_random(cols,rows)
SenvRcult1[] <- 1
# table
SenvRcult_fishnet <- st_make_grid(st_as_sfc(st_bbox(SenvRcult1)), 
                                  cellsize = fisnet_cell_size, square = FALSE) %>%
  st_sf() %>% 
  mutate(fishnet_id = row_number(),
         # count = rpois(nrow(.),8),
         count = sample(0:1,nrow(.), replace = TRUE),    # or uniform across (0,1) ints?
         mean_val1 = 1,
         mean_val2 = 1,
         mean_val3 = 1,
         mean_val4 = 1)
# hex_plot <- gather(st_drop_geometry(SenvRcult_fishnet), id, val, -fishnet_id, -count) %>% 
#   left_join(., SenvRcult_fishnet, by = "fishnet_id") %>% 
#   st_sf() %>% 
#   dplyr::select(id, val)
# show_landscape(list(
#   "mean_val1" = SenvRcult1,
#   "mean_val2" = SenvRcult1,
#   "mean_val3" = SenvRcult1,
#   "mean_val4" = SenvRcult1), unique_scales = TRUE) +
#   geom_sf(data = hex_plot, inherit.aes = FALSE, color = "gray30", alpha = 0)
# 
pairs(as.matrix(st_drop_geometry(SenvRcult_fishnet[,2:6])))



# or is this semi-structured? structured may be totally homogeneous            WORKING!!!
### SEMI-STRUCTURED ENV and RANDOM CULTURAL ------------------------------------ ???????
# raster
SenvRcult_grad <- NLMR::nlm_distancegradient(cols,rows, origin = c(20, 30, 10, 15))
SenvRcult_grad <- abs(1-(SenvRcult_grad+0.01))*2 # invert and weight
SenvRcult1 <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = cols*2) * SenvRcult_grad
SenvRcult1 <- raster::scale(SenvRcult1)
vx_SenvRcult1 <- velox(SenvRcult1) # cast raster to velox
SenvRcult2 <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = cols*2) * SenvRcult_grad
SenvRcult2 <- raster::scale(SenvRcult2)
vx_SenvRcult2 <- velox(SenvRcult2) # cast raster to velox
SenvRcult3 <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = cols*2) * SenvRcult_grad
SenvRcult3 <- raster::scale(SenvRcult3)
vx_SenvRcult3 <- velox(SenvRcult3) # cast raster to velox
SenvRcult4 <- NLMR::nlm_gaussianfield(cols,rows, autocorr_range = cols*2) * SenvRcult_grad
SenvRcult4 <- raster::scale(SenvRcult4)
vx_SenvRcult4 <- velox(SenvRcult4) # cast raster to velox
# table
SenvRcult_fishnet <- st_make_grid(st_as_sfc(st_bbox(SenvRcult1)), 
                                  cellsize = fisnet_cell_size, square = FALSE) %>%
  st_sf() %>% 
  mutate(fishnet_id = row_number(),
         count = 1, # sample(0:1,nrow(.), replace = TRUE),    # or uniform across (0,1) ints?
         mean_val1 = as.numeric(vx_SenvRcult1$extract(., fun = mean, small = TRUE)),
         mean_val2 = as.numeric(vx_SenvRcult2$extract(., fun = mean, small = TRUE)),
         mean_val3 = as.numeric(vx_SenvRcult3$extract(., fun = mean, small = TRUE)),
         mean_val4 = as.numeric(vx_SenvRcult4$extract(., fun = mean, small = TRUE)))

# hex_plot <- gather(st_drop_geometry(SenvRcult_fishnet), id, val, -fishnet_id, -count) %>% 
#   left_join(., SenvRcult_fishnet, by = "fishnet_id") %>% 
#   st_sf() %>% 
#   dplyr::select(id, val)
# 
show_landscape(list(
  "mean_val1" = SenvRcult1,
  "mean_val2" = SenvRcult2,
  "mean_val3" = SenvRcult3,
  "mean_val4" = SenvRcult4), unique_scales = TRUE) #+
#   geom_sf(data = hex_plot, inherit.aes = FALSE, color = "gray30", alpha = 0)
# 
# pairs(as.matrix(st_drop_geometry(SenvRcult_fishnet[,2:6])))



########## Testing ability to detect correlation...
# intersect fishnet to all raster, for each cell compute corr, put more sites in cells with higher corr??

# get all raster values within each cell
x1 <- vx_SenvRcult1$extract(SenvRcult_fishnet, small = TRUE)
x2 <- vx_SenvRcult2$extract(SenvRcult_fishnet, small = TRUE)
x3 <- vx_SenvRcult3$extract(SenvRcult_fishnet, small = TRUE)
x4 <- vx_SenvRcult4$extract(SenvRcult_fishnet, small = TRUE)

# cast raster values to list column of vectors, map vectors to mean correlation, join to fishnet by id
xt <- tibble::enframe(x1) %>% 
  rename("x1" = value) %>% 
  mutate(x2 = tibble::enframe(x2) %>% pull(value),
         x3 = tibble::enframe(x3) %>% pull(value),
         x4 = tibble::enframe(x4) %>% pull(value),
         fishnet_id = as.integer(name)) %>% 
  mutate(cor = pmap_dbl(list(x1,x2,x3,x4), function(a,b,c,d) mean(cor(cbind(a,b,c,d))))) %>% 
  dplyr::select(fishnet_id, cor) %>% 
  left_join(SenvRcult_fishnet, ., by = "fishnet_id")

new <- cbind(xt$fishnet_id, sapply(-st_drop_geometry(xt[,c(3,7)]), rank)) %>% 
  data.frame() %>% 
  mutate(x = 0 - mean_val1 - cor) %>% 
  arrange(desc(x))

mapview::mapview(xt, zcol = "cor")







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
