#' Update the simulation window in a SpatSimObj
#'
#' This function updates the simulation window of a \code{SpatSimObj} by replacing
#' the existing window with a new one.
#'
#' @param sim_object A \code{SpatSimObj} object
#' @param window A new \code{owin} object representing the updated simulation window
#' @return The updated \code{SpatSimObj} object
#' @details The \code{UpdateSimulationWindow()} function checks that the input \code{sim_object} is of class
#' 'SpatSimObj', that the input \code{window} is not null and is of class 'owin'.
#' If these checks pass, the function updates the simulation window in the input \code{sim_object} and
#' returns the updated \code{SpatSimObj} object.
#'
#' @examples
#' # Create a simulation object
#' sim_obj <- CreateSimulationObject()
#'
#' # Update the simulation window
#' new_window <- spatstat.geom::owin(c(0, 5), c(0, 5))
#' updated_sim_obj <- UpdateSimulationWindow(sim_obj, window = new_window)
#'
#' @seealso \code{\link{CreateSimulationObject}}
#'
#' @export

UpdateSimulationWindow = function(sim_object, window = NULL){
  if(!methods::is(sim_object, "SpatSimObj")) stop("`sim_list` must be of class 'SpatSimObj'")
  if(is.null(window)) stop("In order to update the simulation window, window must not be null")
  if(!methods::is(window, "owin")) stop("`window` must be of class 'owin'")
  sim_object@Window = window
  return(sim_object)
}
