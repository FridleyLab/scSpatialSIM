#' Method `summary` for objects in the `mIFsim` package
#'
#' This provides an easy summary of what has already been executed in the Spatial Simulation Object.
#' Things like the number of simulations, number of cell types desired, as well as information about what has been simulated
#' are included in the sumamry output.
#'
#' @param object Spatial Simulation Object
#'
#' @usage \\method{summary}{`Spatial Simulation Object`}
#'
#' @export
setMethod("summary",
          signature = "Spatial Simulation Object",
          function(object) {
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
          })

#' Method `plot` for objects of class `Spatial Simulation Object`
#'
#' When proceeding though the simulation of multiplex data, it is nice to view some of the processes and results of the process.
#' Using the plot function, the `plot` from spatstat is used (in the case of plotting just the point process). When plotting
#' more detail from the `Spatial Simulation Object`, please use `PlotSimulation`
#'
#' @param x Spatial Simulation Object
#' @param nrow number of rows in the output plot
#' @param ncol number of columns in the output plot
#' @param which which spatial pattern(s) to plot
#' @param what currently only accepting "Processes" but potentially more things using the spatstat plot function in the future
#'
#' @usage \\method{plot}{`Spatial Simulation Object`}
#'
#' @export
setMethod("plot",
          signature = "Spatial Simulation Object",
          function(x, nrow = 1, ncol = 1, which = 1, what = "Processes"){ #
            if(what == "Processes"){
              if(length(which) == 1){
                plot(x@Processes[[which]], main = which)

              } else {
                par(mfrow = c(nrow, ncol))
                hld = lapply(which, function(i){
                  plot(x@Processes[[i]], main = i)
                })
              }
            } else {
              cat("please be patient while method is being updated")
            }
          })
