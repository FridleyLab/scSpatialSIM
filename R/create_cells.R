#' Generate Cell Positivity
#'
#' Generate the probability of a cell being positive given a set of simulation parameters for each file in a SpatSimObj.
#'
#' @param sim_object A \code{SpatSimObj} object containing the simulated data.
#' @param k An integer specifying the number of clusters for each simulated patterns
#' @param xmin A numeric value specifying the minimum x value for the kernel.
#' @param xmax A numeric value specifying the maximum x value for the kernel.
#' @param ymin A numeric value specifying the minimum y value for the kernel.
#' @param ymax A numeric value specifying the maximum y value for the kernel.
#' @param sdmin A numeric value specifying the minimum standard deviation for the kernel.
#' @param sdmax A numeric value specifying the maximum standard deviation for the kernel.
#' @param probs A numeric vector of length 2 specifying the minimum and maximum probability values for scaling the kernel values.
#' @param Force A logical value indicating whether to force simulation parameters to be within the simulation window limits.
#' @param density_heatmap A logical value indicating whether to compute a density heatmap for each cell.
#' @param step_size A numeric value specifying the step size for the grid of points within the window.
#' @param cores An integer value specifying the number of cores to use for parallel computation.
#' @param correlation A value between -1 and 1 for how related a second or more cell type is to the first
#'
#' @return Returns the original \code{sc.SpatialSIM} object with additional generated data added to each cell object.
#'
#' @details The function generates the probability of a cell being positive given a set of simulation parameters f
#' or each file in a \code{sc.SpatialSIM} object. It creates a kernel parameter list for \code{k} clusters
#' in each simulated pattern and computes the probability of each point in the grid of points within the
#' window for each cell. The function also computes a density heatmap for each cell if \code{density_heatmap} is set to \code{TRUE}.
#'
#'
#' @export
GenerateCellPositivity = function(sim_object, k = NA,
                          xmin = NA, xmax = NA, ymin = NA, ymax = NA,
                          sdmin = 1/2, sdmax = 2, probs = c(0, 1),
                          Force = FALSE,
                          density_heatmap = FALSE, step_size = 1, cores = 1,
                          correlation = 0){
  if(!methods::is(sim_object, "SpatSimObj")) stop("`sim_object` must be of class 'SpatSimObj'")
  if(any(is.null(c(k, xmin, xmax, ymin, ymax, sdmin, sdmax)))) stop("Cannot have `NULL` parameters")

  #create parameter vector
  params = list(k = k,
             xmin = xmin, xmax = xmax,
             ymin = ymin, ymax = ymax,
             sdmin = sdmin, sdmax = sdmax,
             probs = probs)
  ncells = length(sim_object@Cells)
  #dummy variable to prevent console printing
  dmb = lapply(seq(ncells), function(cell){
    #if no parameters are input then use the initialized
    params = mapply(replace_na, sim_object@Cells[[cell]]@Parameters, params, SIMPLIFY = FALSE)
    #add updated parameters to the object cell
    sim_object@Cells[[cell]]@Parameters <<- params
    #get the window size
    win_limits = c(sim_object@Window$xrange, sim_object@Window$yrange)
    #check whether the parameters would simulate outside window
    if(any((unlist(params[c(2, 4)]) < win_limits[c(1,3)]) |
           (unlist(params[c(3, 5)]) < win_limits[c(2,4)])) & Force == FALSE){
      stop("x and y range outside simulation window limits")
    }
    #inform user parameter window inside simulation window
    if(any(c(unlist(params[c(2, 4)]) > win_limits[c(1,3)],
             unlist(params[c(3, 5)]) > win_limits[c(2,4)]))){
      message("x and y range inside window boundary")
    }
    #produce kernel parameter list for k clusters in each simulated pattern
    if(cell == 1){
      sim_object@Cells[[cell]]@`Simulationed Kernels` <<- lapply(seq(sim_object@Sims), function(hld){
        do.call(gaussian_kernel, utils::head(params, -1))
      })
    } else {
      sim_object@Cells[[cell]]@`Simulationed Kernels` <<- sim_object@Cells[[1]]@`Simulationed Kernels`
    }

    #make gric based on user step size for their window
    grid = expand.grid(x = seq(win_limits[1], win_limits[2], step_size),
                       y = seq(win_limits[3], win_limits[4], step_size))

    if(density_heatmap){
      if(cell == 1){
        message(paste0("Computing density heatmap for Cell ", cell))
        sim_object@Cells[[cell]]@`Density Grids` <<- pbmcapply::pbmclapply(sim_object@Cells[[cell]]@`Simulationed Kernels`, function(gauss_tab){
          cbind(grid,
                prob = CalculateGrid(grid, gauss_tab, cores = cores))
        })
      } else {
        message(paste0("Adjusting density heatmap for Cell ", cell))
        sim_object@Cells[[cell]]@`Density Grids` <<- pbmcapply::pbmclapply(sim_object@Cells[[1]]@`Density Grids`, function(grid_tab){
          if(correlation == 0){
            grid_tab$prob = stats::runif(nrow(grid_tab), min = 0, max = 1)
          } else if(correlation < 0){
            grid_tab$prob = grid_tab$prob * correlation + 1
          } else {
            grid_tab$prob = grid_tab$prob * correlation
          }
          return(grid_tab)
        })
      }


    }

    if(is.empty(sim_object, "Spatial Files")){
      sim_object@`Spatial Files` <<- lapply(sim_object@Patterns, data.frame)
    }

    message(paste0("Computing probability for Cell ", cell))
    sim_object@`Spatial Files` <<- pbmcapply::pbmclapply(seq(sim_object@`Spatial Files`), function(spat_num){
      vec = CalculateGrid(sim_object@`Spatial Files`[[spat_num]],
                                                      sim_object@Cells[[cell]]@`Simulationed Kernels`[[spat_num]], cores = cores)
      #if the cell is other than the first, adjust it based on first cell and correlation
      if(cell != 1){
        if(correlation == 0){
          vec = stats::runif(length(vec), min = 0, max = 1)
        } else if(correlation < 0){
          vec = vec * correlation + 1
        } else {
          vec = vec * correlation
        }
      }
      #make table with probabilities and positive/negative
      df = data.frame(col1 = scale_probs(vec * 0.9, params$probs))
      df$col2 = ifelse(stats::rbinom(nrow(df), size = 1, prob = df$col1) == 1, "Positive", "Negative")

      names(df) = c(paste("Cell", cell, "Probability"), paste("Cell", cell, "Assignment"))

      return(cbind(sim_object@`Spatial Files`[[spat_num]],
                   df))
    })
  })

  return(sim_object)
}
