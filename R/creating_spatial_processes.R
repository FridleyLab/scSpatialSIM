#' Generate Spatial Point Pattern
#'
#' Generate a spatial point pattern within the simulation object's window using a Poisson point process.
#'
#' @param sim_object A 'SpatSimObj' containing a window.
#' @param lambda The intensity of the point pattern Default is 25.
#' @param ... Additional arguments passed to 'rpoispp'.
#'
#' @return The updated 'sim_object' with a simulated point process added to the 'Processes' slot.
#'
#' @details This function generates a spatial point process within the window of the 'sim_object'
#'  using a Poisson point pattern with intensity 'lambda'. The simulated point pattern is added
#'  to the 'Patterns' slot of the 'sim_object'. Additional arguments can be passed to the 'rpoispp' function.
#'
#' @export
#' @examples
#' sim_object <- CreateSimulationObject()
#' sim_object <- GenerateSpatialPattern(sim_object, lambda = 30)
GenerateSpatialPattern = function(sim_object, lambda = 25, ...){
  if(!methods::is(sim_object, "SpatSimObj")) stop("`sim_object` must be of class 'SpatSimObj'")
  if(is.null(lambda)) stop("Need an intensity in order to simulate points")

  window = sim_object@Window
  sims = sim_object@Sims
  spatial_dfs = lapply(seq(sims), function(hld){
    spatstat.random::rpoispp(lambda, win = window, nsim = 1, ...)
  })

  sim_object@Patterns = spatial_dfs
  return(sim_object)
}
