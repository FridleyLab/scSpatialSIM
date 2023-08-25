#' Generate Tissue
#'
#' This function generates a simulated tissue using a specified number of
#' clusters and spatial parameters for each pattern in the simulation object.
#' The tissue is represented by a grid of points with probabilities of
#' belonging to tissue 1 or tissue 2, based on a Gaussian kernel density
#' estimate calculated for each pattern
#'
#' @param sim_object A `SpatSimObj` created with
#'   \code{\link{CreateSimulationObject}}.
#' @param k Number of clusters to generate for each pattern
#' @param xmin Minimum x-coordinate for cluster centers.
#' @param xmax Maximum x-coordinate for cluster centers.
#' @param ymin Minimum y-coordinate for cluster centers.
#' @param ymax Maximum y-coordinate for cluster centers.
#' @param sdmin Minimum standard deviation for cluster kernels.
#' @param sdmax Maximum standard deviation for cluster kernels.
#' @param force Logical, whether to force generation of tissue even if the
#'   generated cluster centers would fall outside the simulation window. If
#'   \code{FALSE}, an error will be thrown if cluster centers are outside the
#'   window.
#' @param density_heatmap Logical, whether to calculate a density heatmap for
#'   the simulated tissue. If \code{TRUE}, a grid of points will be generated
#'   covering the entire simulation window, and the probability of each grid
#'   point belonging to tissue 1 will be calculated based on the generated tissue
#'   probability.
#' @param step_size Grid step size for the density heatmap.
#' @param cores Number of cores to use for parallel processing of density
#'   calculations.
#' @param overwrite boolean whether to overwrite if tissue kernels already exist
#'
#' @return A modified 'Spatial Simulation Object' with updated tissue grids and
#'   assigned tissue types for each simulated pattern.
#'
#' @details This function generates a simulated tissue for each pattern in the
#'   simulation object by first generating k clusters within the specified x
#'   and y ranges and with a standard deviation within the specified range.
#'   Then, a Gaussian kernel density estimate is calculated for each pattern
#'   using the generated clusters as center points and the specified standard
#'   deviation as kernel size. The density estimates represent the probability
#'   of each point in the simulation window belonging to tissue 1 or tissue 2.
#'   If \code{density_heatmap = TRUE}, a density heatmap will be
#'   calculated using a grid of points covering the entire simulation window.
#'   Finally, for each simulated point, the probability of belonging to
#'   tissue 1 is calculated based on the kernel density estimate, and the tissue
#'   type is assigned with
#'   probability proportional to the probability of belonging to tissue 1.
#'
#' @examples
#' # Create a simulation object with a window and point pattern
#' sim_object <- CreateSimulationObject()
#'
#' #simulate points
#' sim_object <- GenerateSpatialPattern(sim_object, lambda = 20)
#'
#' # Generate tissue with default parameters
#' sim_object <- GenerateTissue(sim_object)
#'
#' @export
GenerateTissue = function(sim_object, k = NA,
                          xmin = NA, xmax = NA, ymin = NA, ymax = NA,
                          sdmin = 1/2, sdmax = 2,
                          force = FALSE,
                          density_heatmap = FALSE, step_size = 1, cores = 1, overwrite = FALSE){
  #stop conditions
  if(!methods::is(sim_object, "SpatSimObj")) stop("`sim_object` must be of class 'SpatSimObj'")
  if(any(is.null(c(k, xmin, xmax, ymin, ymax, sdmin, sdmax)))) stop("Cannot have `NULL` parameters")

  if(!is.empty(sim_object@Tissue, "Simulated Kernels") & overwrite == FALSE) stop("Already have tissue kernels and `overwrite == FALSE`")
  if(!is.empty(sim_object@Tissue, "Simulated Kernels") & overwrite == TRUE){
    message("Overwriting existing tissue kernels")
    message("Resetting Tissue slots")
    #tissue
    sim_object@Tissue@`Simulated Kernels` = list()
    sim_object@Tissue@`Density Grids` = list()
    # #holes
    # sim_object@Holes@`Simulated Kernels` = list()
    # sim_object@Holes@`Density Grids` = list()
    # #cells
    # for(i in seq(sim_object@Cells)){
    #   sim_object@Cells[[i]]@`Simulated Kernels` = list()
    #   sim_object@Cells[[i]]@`Density Grids` = list()
    # }
    #spatial files
    sim_object@`Spatial Files` = lapply(sim_object@`Spatial Files`, function(spat){
      spat %>%
        dplyr::select(-dplyr::contains("Tissue"))
    })
    #letting know finished
    message("Reset...Continuing.")
  }

  #create parameter vector
  params = list(k = k,
             xmin = xmin, xmax = xmax,
             ymin = ymin, ymax = ymax,
             sdmin = sdmin, sdmax = sdmax)
  #if no parameters are input then use the initialized
  params = mapply(replace_na, sim_object@Tissue@Parameters, params, SIMPLIFY = FALSE)
  #update initialized paramters with custom input from user
  sim_object@Tissue@Parameters <- params
  #get the window size
  win_limits = c(sim_object@Window$xrange, sim_object@Window$yrange)
  #check whether the parameters would simulate outside window
  if(any((unlist(params[c(2, 4)]) < win_limits[c(1,3)]) |
         (unlist(params[c(3, 5)]) > win_limits[c(2,4)])) & force == FALSE){
    stop("x and y range outside simulation window limits")
  }
  #inform user parameter window inside simulation window
  if(any(c(unlist(params[c(2, 4)]) > win_limits[c(1,3)],
           unlist(params[c(3, 5)]) < win_limits[c(2,4)]))){
    message("x and y range inside window boundary")
  }
  #produce kernel parameter list for k clusters in each simulated pattern
  sim_object@Tissue@`Simulated Kernels` = lapply(seq(sim_object@Sims), function(hld){
    do.call(gaussian_kernel, params)
  })
  #make gric based on user step size for their window
  grid = expand.grid(x = seq(win_limits[1], win_limits[2], step_size),
                     y = seq(win_limits[3], win_limits[4], step_size))

  if(density_heatmap){
    message("Computing density heatmap")
    sim_object@Tissue@`Density Grids` = pbmcapply::pbmclapply(sim_object@Tissue@`Simulated Kernels`, function(gauss_tab){
      cbind(grid,
            prob = CalculateGrid(grid, gauss_tab, cores = cores))
    })
  }

  if(is.empty(sim_object, "Spatial Files")){
    sim_object@`Spatial Files` = lapply(sim_object@Patterns, data.frame)
  }

  message("Computing tissue probability")
  sim_object@`Spatial Files` = pbmcapply::pbmclapply(seq(sim_object@`Spatial Files`), function(spat_num){
    df = cbind(sim_object@`Spatial Files`[[spat_num]],
          `Tissue Probability` = CalculateGrid(sim_object@`Spatial Files`[[spat_num]],
                               sim_object@Tissue@`Simulated Kernels`[[spat_num]], cores = cores) * 0.9)
    df$`Tissue Assignment` = ifelse(stats::rbinom(nrow(df), size = 1, prob = df$`Tissue Probability`) == 1, "Tissue 1", "Tissue 2")
    return(df)
  })

  return(sim_object)
}
