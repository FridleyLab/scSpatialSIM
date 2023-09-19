#' Generate Spatial Point Pattern
#'
#' Generate a spatial point pattern within the simulation object's window using a Poisson point process.
#'
#' @param sim_object A 'SpatSimObj' containing a window.
#' @param lambda The intensity of the point pattern Default is 25.
#' @param ... Additional arguments passed to 'rpoispp'.
#' @param overwrite boolean indicating whether or not to replace point patterns if they exist in object
#' @param gridded boolean value to whether or not simulate the point pattern in a grid. See details for more.
#' @param grid_shift the amount to move alternative columns down when gridded; between -0.5 and 0.5
#'
#' @return The updated 'sim_object' with a simulated point process added to the 'Processes' slot.
#'
#' @details This function generates a spatial point process within the window of the 'sim_object'
#'  using a Poisson point pattern with intensity 'lambda'. The simulated point pattern is added
#'  to the 'Patterns' slot of the 'sim_object'. Additional arguments can be passed to the 'rpoispp' function.
#'
#'  The `gridded` parameter is used for simulating point patterns that would represent som spatial transcriptomic
#'  technologies such as visium, where rather than like cell being randomly distributed in a sample, the spots where
#'  data is measured is evenly spaced.
#'
#' @export
#' @examples
#' sim_object <- CreateSimulationObject()
#' sim_object <- GenerateSpatialPattern(sim_object, lambda = 30)
GenerateSpatialPattern = function(sim_object, lambda = 25, ..., overwrite = FALSE, gridded = FALSE, grid_shift = 0.5){
  # stop conditions
  if(!methods::is(sim_object, "SpatSimObj")) stop("`sim_object` must be of class 'SpatSimObj'")
  if(is.null(lambda)) stop("Need an intensity in order to simulate points")
  if(!is.empty(sim_object, "Patterns") & overwrite == FALSE) stop("Already have point patterns and `overwrite == FALSE`")
  if(!is.empty(sim_object, "Patterns") & overwrite == TRUE) message("Overwriting existing point patterns")


  window = sim_object@Window
  sims = sim_object@Sims

  if(gridded){
    if(grid_shift < -0.5 | grid_shift > 0.5){
      stop("Max shift in either positive or negative direction is 0.5")
    }
    #find the area to calculate the number of cells
    area = spatstat.geom::area(window)
    #get the window size to identify the region for simulation
    x_size = diff(window$xrange)
    y_size = diff(window$yrange)
    #calculate the number of cells with lambda and the area
    n_cells = lambda * area #should work!
    #calculate cells per side
    n_cells_side = floor(sqrt(n_cells))
    #identify the spacing in both directions
    x_spacing = x_size / (n_cells_side + 1)
    y_spacing = y_size / (n_cells_side + 1)
    #find x_positions
    x_pos = seq(x_spacing/2, x_size, by = x_spacing)
    #find y_positions
    y_pos = seq(y_spacing/2, y_size, by = y_spacing)
    #get x indices where y needs to shift
    x_vals = x_pos[seq(1, length(x_pos), by = 2)]
    #expand grid for positions
    new_raw_df = data.frame(expand.grid(x = x_pos, y = y_pos)) %>%
      dplyr::mutate(y = ifelse(x %in% x_vals, y + (y_spacing * grid_shift), y))

    #make the list of spatial dataframes like spatstat
    spatial_dfs = lapply(seq(sims), function(hld){
      sp_df = list(window = window,
           n = nrow(new_raw_df),
           x = new_raw_df$x + min(window$xrange),
           y = new_raw_df$y + min(window$yrange),
           markformat = 'none')
      class(sp_df) = 'ppp'
      return(sp_df)
    })
  } else {
    spatial_dfs = lapply(seq(sims), function(hld){
      spatstat.random::rpoispp(lambda, win = window, nsim = 1, ...) %>%
        data.frame()
    })
  }



  sim_object@Patterns = spatial_dfs
  return(sim_object)
}
