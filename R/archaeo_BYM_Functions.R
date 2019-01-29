st_drop_geometry <- function(x) {
  if(inherits(x,"sf")) {
    x <- st_set_geometry(x, NULL)
    class(x) <- 'data.frame'
  }
  return(x)
}
ll <- function(dat, crs = 4326){
  st_transform(dat, crs)
}