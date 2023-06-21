#' Update the simulation window in a Spatial Simulation Object
#'
#' This function updates the simulation window of a \code{Spatial Simulation Object} by replacing
#' the existing window with a new one.
#'
#' @param sim_object A \code{Spatial Simulation Object} object
#' @param window A new \code{owin} object representing the updated simulation window
#' @return The updated \code{Spatial Simulation Object} object
#' @details The \code{UpdateSimulationWindow()} function checks that the input \code{sim_object} is of class
#' 'Spatial Simulation Object', that the input \code{window} is not null and is of class 'owin'.
#' If these checks pass, the function updates the simulation window in the input \code{sim_object} and
#' returns the updated \code{Spatial Simulation Object} object.
#'
#' @examples
#' # Create a simulation object
#' sim_obj <- create_simObject()
#'
#' # Update the simulation window
#' new_window <- owin(c(0, 5), c(0, 5))
#' updated_sim_obj <- UpdateSimulationWindow(sim_obj, window = new_window)
#'
#' @seealso \code{\link{create_simObject}}
#'
#' @export

UpdateSimulationWindow = function(sim_object, window = NULL){
  if(!is(sim_object, "Spatial Simulation Object")) stop("`sim_list` must be of class 'Spatial Simulation Object'")
  if(is.null(window)) stop("In order to update the simulation window, window must not be null")
  if(!is(window, "owin")) stop("`window` must be of class 'owin'")
  sim_object@Window = window
  return(sim_object)
}
