#' summary function for SpatialSimulationObject
#'
#' @param object of class `SpatialSimulationObject`
#' @param ... nothing else to pass to summary if object is a SpatialSimulationObject
#'
#' @method summary SpatialSimulationObject
#'
#' @export
summary.SpatialSimulationObject <- function(object, ...){
  nsim = object@Sims
  win = matrix(c(object@Window$xrange, object@Window$yrange), nrow = 2, byrow = F)
  s_processes = length(object@Processes)
  t_kernels = length(object@Tissue@`Simulationed Kernels`)
  h_kernels = length(object@Holes@`Simulationed Kernels`)
  c_types = length(object@Cells)
  c_kernels=length(object@Cells[[1]]@`Simulationed Kernels`)

  cat("Spatial Simulation object for", nsim,"simulations. Currently, there are:\n")
  cat(paste0("\tWindow: x (", win[1], ",", win[2], "); y (", win[3], ",", win[4],")\n"), sep="")
  cat("\t", s_processes, " spatial processes\n", sep="")
  cat("\t", t_kernels, " tissue kernels\n", sep="")
  cat("\t", h_kernels, " hole kernels\n", sep="")
  cat("\t", c_kernels, " cell kernels for ", c_types, " cell types\n", sep="")
}

#' plot function for SpatialSimulationObject
#'
#' @param x of class `SpatialSimulationObject`
#' @param ... other things to pass to the plot method for SpatialSimulationObjects including `nrow`, `ncol` for the number of rows and columns of plots,
#' `which` processes to plot, and `what` which currently only works with "Processes" but may be updated in the future
#'
#' @method plot SpatialSimulationObject
#'
#'
#' @export
plot.SpatialSimulationObject <- function(x, ...){ #
  if(!exists("what")) what = "Processes"
  if(methods::is(nrow, "function")) nrow = 1
  if(methods::is(ncol, "function")) ncol = 1
  if(methods::is(which, "function")) which = 1

  if(what == "Processes"){
    if(length(which) == 1){
      plot(x@Processes[[which]], main = which)

    } else {
      graphics::par(mfrow = c(nrow, ncol))
      hld = lapply(which, function(i){
        plot(x@Processes[[i]], main = i)
      })
    }
  } else {
    cat("please be patient while method is being updated")
  }
}
