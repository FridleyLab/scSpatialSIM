#' Create a spatial simulation object.
#'
#' This function creates a \code{SpatSimObj} for
#' spatial simulations. The object contains information about the
#' simulation window, the number of simulations to perform, and
#' lists of cells, tumor/stroma, holes, and spatial files.
#'
#' @details
#'
#' The simulation window is represented by an object of class
#' \code{owin}, which specifies the extent and shape of the spatial
#' domain in which the simulations will be performed. If no window
#' is provided, the function creates a rectangular window of size
#' (0,10) in both x and y directions.
#'
#' The \code{sims} argument specifies the number of simulations to
#' perform. If it is set to \code{NULL} or less than 1, the function
#'  defaults to 3.
#'
#' The \code{cell_types} argument specifies the number of cell
#' types to include in the simulation. By default, the function
#'  creates a single cell type, represented by an object of class
#'   \code{Cell}.
#'
#' The \code{SpatSimObj} is composed of the following classes:
#'
#' \itemize{
#'  \item A \code{Window} object of class \code{owin}.
#'  \item An integer \code{Sims} specifying the number of simulations to perform.
#'  \item A list of \code{Cells} of class \code{Cell}.
#'  \item A \code{Tissue} object of class \code{Tumor/Stroma}, representing the tumor and stroma components of the simulation.
#'  \item A \code{Holes} object of class \code{Holes}, representing holes in the simulation.
#'  \item A list of \code{Spatial Files} containing any spatial data associated with the simulation.
#' }
#'
#' @param window An object of class \code{owin} representing the simulation
#'  window. If \code{NULL}, defaults to a rectangular window of size (0,10)
#'  in both x and y directions.
#' @param sims The number of simulations to perform. If \code{NULL} or l
#' ess than 1, defaults to 3.
#' @param cell_types The number of cell types. Defaults to 1.
#'
#' @return A \code{SpatSimObj} containing the simulation
#' window, the number of simulations to perform, and lists of cells,
#' tumor/stroma, holes, and spatial files.
#'
#' @export
#'
#' @examples
#' CreateSimulationObject()
CreateSimulationObject = function(window = NULL,
                            sims = NULL,
                            cell_types = 1){
  #require(spatstat)
  if(is.null(window)){
    message("No `window` specified - defaulting to x (0, 10); y (0, 10)")
    window = spatstat.geom::owin(c(0,10), c(0,10))
  }
  if(is.null(cell_types)){
    stop("Please create provide at least 1 cell types.")
  }

  if(is.null(sims) || sims<1){
    sims = 3
  }

  if(!methods::is(window, "owin")){
    sims = length(window)
  }

  datClass = methods::new("SpatSimObj", Window = window, Sims = sims)
  c_list = lapply(seq(cell_types), function(x) methods::new("Cell"))
  datClass@Cells = c_list
  return(datClass)
}



setClass("Points", slot = list(
  Patterns = "list"
),
prototype = list(
  Patterns = list()
))

setClass("Tumor/Stroma", slots = list(
  Parameters = "list",
  `Simulationed Kernels` = "list",
  `Density Grids` = "list"
),
prototype = list(
  Parameters = list(k = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, sdmin = 1/2, sdmax = 2),
  `Simulationed Kernels` = list(),
  `Density Grids` = list()
))

setClass("Holes", slots = list(
  Parameters = "list",
  `Simulationed Kernels` = "list",
  `Density Grids` = "list"
),
prototype = list(
  Parameters = list(xmin = 0, xmax = 10, ymin = 0, ymax = 10, sdmin = 1/2, sdmax = 2, hole_prob = c(0.2, 0.35)),
  `Simulationed Kernels` = list(),
  `Density Grids` = list()
))

setClass("Cell", slots = list(
  Parameters = "list",
  `Simulationed Kernels` = "list",
  `Density Grids` = "list"
),
prototype = list(
  Parameters = list(k = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, sdmin = 1/2, sdmax = 2, probs = c(0, 1)),
  `Simulationed Kernels` = list(),
  `Density Grids` = list()
))

setClass("owin", slots = list(
  type = "character",
  xrange = "numeric",
  yrange = "numeric",
  units = "list"
))

setClass("SpatSimObj", slots = list(
  Window = "owin",
  Sims = "numeric",
  Patterns = "list",
  `Tissue` = "Tumor/Stroma",
  Holes = "Holes",
  Cells = "list",
  `Spatial Files` = "list"
))


