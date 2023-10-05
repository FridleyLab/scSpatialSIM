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
#' @param probs Either a vector of c(low probability, high probability) for all cell types or data frame where each row
#' is the low and high probabilities for the cell type. If data frame, number of rows must equal number of cells
#' @param Force A logical value indicating whether to force simulation parameters to be within the simulation window limits.
#' @param density_heatmap A logical value indicating whether to compute a density heatmap for each cell.
#' @param step_size A numeric value specifying the step size for the grid of points within the window.
#' @param cores An integer value specifying the number of cores to use for parallel computation.
#' @param shift A value between 0 and 1 for how related a second or more cell type is to the first
#' @param random whether or not to randomly generate kernels for cells 2 or more, uf TRUE, shift is not used
#' @param overwrite boolean whether to overwrite existing cell kernels and assignments if present
#' @param use_window boolean whether to use the simulation window to set x and y limits
#'
#' @return Returns the original \code{scSpatialSIM} object with additional generated data added to each cell object.
#'
#' @details The function generates the probability of a cell being positive given a set of simulation parameters f
#' or each file in a \code{scSpatialSIM} object. It creates a kernel parameter list for \code{k} clusters
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
                          shift = 0, random = FALSE, overwrite = FALSE,
                          use_window = FALSE){
  if(!methods::is(sim_object, "SpatSimObj")) stop("`sim_object` must be of class 'SpatSimObj'")
  if(any(is.null(c(k, xmin, xmax, ymin, ymax, sdmin, sdmax)))) stop("Cannot have `NULL` parameters")

  if(!is.empty(sim_object@Cells[[1]], "Simulated Kernels") & overwrite == FALSE) stop("Already have cell kernels and `overwrite == FALSE`")

  ncells = length(sim_object@Cells)
  if(methods::is(probs, "data.frame"))
    if(nrow(probs) != ncells) stop("`probs` should either be data.frame with nrow length of Cells or a vector of length 2")
  if(!is.empty(sim_object@Cells[[1]], "Simulated Kernels") & overwrite == TRUE){
    message("Overwriting existing cell kernels")
    # #tissue
    # sim_object@Tissue@`Simulated Kernels` = list()
    # sim_object@Tissue@`Density Grids` = list()
    # #holes
    # sim_object@Holes@`Simulated Kernels` = list()
    # sim_object@Holes@`Density Grids` = list()
    #cells
    message("Resetting Cell slots")
    for(i in seq(sim_object@Cells)){
      sim_object@Cells[[i]]@`Simulated Kernels` = list()
      sim_object@Cells[[i]]@`Density Grids` = list()
    }
    #spatial files
    sim_object@`Spatial Files` = lapply(sim_object@`Spatial Files`, function(spat){
      spat %>%
        dplyr::select(-dplyr::contains("Cell"))
    })
    #letting know finished
    message("Reset...Continuing.")
  }

  #check if using window is TRUE
  if(use_window){
    xmin = sim_object@Window$xrange[1]
    xmax = sim_object@Window$xrange[2]
    ymin = sim_object@Window$yrange[1]
    ymax = sim_object@Window$yrange[2]
  }

  #create parameter vector
  if(methods::is(probs, "vector")){
    probs2 = data.frame(matrix(rep(probs, ncells), nrow = ncells, byrow = TRUE)) %>% stats::setNames(c("Low", "High"))
  } else {
    probs2 = probs
  }

  params_overall = list(k = k,
                        xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax,
                        sdmin = sdmin, sdmax = sdmax,
                        probs = probs2)

  #make sure that the shift is in bounds
  if((shift < 0 | shift > 1) & ncells > 1) stop("supply an appropriate shift")
  #dummy variable to prevent console printing
  dmb = lapply(seq(ncells), function(cell){
    #subset the params to the specific cells parameters
    params = params_overall
    params$probs = params$probs[cell,] %>% as.numeric()
    #if no parameters are input then use the initialized
    params = mapply(replace_na, sim_object@Cells[[cell]]@Parameters, params, SIMPLIFY = FALSE)
    #add updated parameters to the object cell
    sim_object@Cells[[cell]]@Parameters <<- params
    #get the window size
    win_limits = c(sim_object@Window$xrange, sim_object@Window$yrange)
    #check whether the parameters would simulate outside window
    if(any((unlist(params[c(2, 4)]) < win_limits[c(1,3)]) |
           (unlist(params[c(3, 5)]) > win_limits[c(2,4)])) & Force == FALSE){
      stop("x and y range outside simulation window limits")
    }
    #inform user parameter window inside simulation window
    if(any(c(unlist(params[c(2, 4)]) > win_limits[c(1,3)],
             unlist(params[c(3, 5)]) < win_limits[c(2,4)]))){
      message("x and y range inside window boundary")
    }
    #produce kernel parameter list for k clusters in each simulated pattern
    if(cell == 1 | random == TRUE){
      sim_object@Cells[[cell]]@`Simulated Kernels` <<- lapply(seq(sim_object@Sims), function(hld){
        do.call(gaussian_kernel, utils::head(params, -1))
      })
    } else {
      #shift kernel from initial if wanted otherwise random make new ones?
      if(shift == 0){
        sim_object@Cells[[cell]]@`Simulated Kernels` <<- sim_object@Cells[[1]]@`Simulated Kernels`
      } else {
        sim_object@Cells[[cell]]@`Simulated Kernels` <<- lapply(seq(sim_object@Sims), function(hld){
          kern = sim_object@Cells[[1]]@`Simulated Kernels`[[hld]]
          #if there are less than 3 centers, just return the kernel
          #will invert below
          #if(nrow(kern)<3) return(kern)
          #if more than 3, we can move the centers
          gaussian_kernel_shift(kern, shift, win_limits)
          #do.call(gaussian_kernel_shift, utils::head(params, -1))
        })
      }

    }

    #make gric based on user step size for their window
    grid = expand.grid(x = seq(win_limits[1], win_limits[2], step_size),
                       y = seq(win_limits[3], win_limits[4], step_size))

    if(density_heatmap){
      if(cell == 1 | shift != 0 | random == TRUE){
        message(paste0("Computing density heatmap for Cell ", cell))
        sim_object@Cells[[cell]]@`Density Grids` <<- pbmcapply::pbmclapply(sim_object@Cells[[cell]]@`Simulated Kernels`, function(gauss_tab){
          cbind(grid,
                prob = CalculateGrid(grid, gauss_tab, cores = cores))
        })
      } else {
        message(paste0("Copying density heatmap for Cell ", cell))
        sim_object@Cells[[cell]]@`Density Grids` <<- pbmcapply::pbmclapply(sim_object@Cells[[1]]@`Density Grids`, function(grid_tab){
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
                                                      sim_object@Cells[[cell]]@`Simulated Kernels`[[spat_num]], cores = cores)
      #if the cell is other than the first, adjust it based on first cell and correlation
      # if(cell != 1){
      #   if(correlation == 0){
      #     vec = stats::runif(length(vec), min = 0, max = 1)
      #   } else if(correlation < 0){
      #     vec = vec * correlation + 1
      #   } else {
      #     vec = vec * correlation
      #   }
      # }
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
