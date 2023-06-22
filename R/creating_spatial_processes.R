#' Generate Spatial Point Process
#'
#' Generate a spatial point process within the simulation object's window using a Poisson point process.
#'
#' @param sim_object A 'Spatial Simulation Object' containing a window.
#' @param lambda The intensity of the point process. Default is 25.
#' @param ... Additional arguments passed to 'rpoispp'.
#'
#' @return The updated 'sim_object' with a simulated point process added to the 'Processes' slot.
#'
#' @details This function generates a spatial point process within the window of the 'sim_object'
#'  using a Poisson point process with intensity 'lambda'. The simulated point process is added
#'  to the 'Processes' slot of the 'sim_object'. Additional arguments can be passed to the 'rpoispp' function.
#'
#' @export
#' @examples
#' sim_object <- create_simObject()
#' sim_object <- GenerateSpatialProcess(sim_object, lambda = 30)
GenerateSpatialProcess = function(sim_object, lambda = 25, ...){
  if(!is(sim_object, "Spatial Simulation Object")) stop("`sim_object` must be of class 'Spatial Simulation Object'")
  if(is.null(lambda)) stop("Need an intensity in order to simulate points")

  window = sim_object@Window
  sims = sim_object@Sims
  spatial_dfs = lapply(seq(sims), function(hld){
    spatstat.random::rpoispp(lambda, win = window, nsim = 1, ...)
  })

  sim_object@Processes = spatial_dfs
  return(sim_object)
}
