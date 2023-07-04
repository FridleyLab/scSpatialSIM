## code to prepare `existing_window_function` dataset goes here

usethis::use_data(existing_window_function, overwrite = TRUE)

example_window_function = function(){
  w <- spatstat.geom::owin(c(-1,1), c(-1,1), mask=matrix(TRUE, 100,100))
  X <- spatstat.geom::raster.x(w)
  Y <- spatstat.geom::raster.y(w)
  wm <- spatstat.geom::owin(w$xrange, w$yrange, mask=(X^2 + Y^2 <= 1))
  return(wm)
}
