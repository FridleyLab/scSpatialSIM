#' summary function for SpatSimObj
#'
#' @param object of class `SpatSimObj`
#' @param ... nothing else to pass to summary if object is a SpatSimObj
#'
#' @method summary SpatSimObj
#'
#' @export
summary.SpatSimObj <- function(object, ...){
  nsim = object@Sims
  win = matrix(c(object@Window$xrange, object@Window$yrange), nrow = 2, byrow = F)
  s_patterns = length(object@Patterns)
  t_kernels = length(object@Tissue@`Simulated Kernels`)
  h_kernels = length(object@Holes@`Simulated Kernels`)
  c_types = length(object@Cells)
  c_kernels=length(object@Cells[[1]]@`Simulated Kernels`)

  cat("Spatial Simulation object for", nsim,"simulated images. Currently, there are:\n")
  cat(paste0("\tWindow: x (", win[1], ",", win[2], "); y (", win[3], ",", win[4],")\n"), sep="")
  cat("\t", s_patterns, " spatial point patterns\n", sep="")
  cat("\t", t_kernels, " tissue kernels\n", sep="")
  cat("\t", h_kernels, " hole kernels\n", sep="")
  cat("\t", c_kernels, " cell kernels for ", c_types, " cell types\n", sep="")
}

#' plot function for SpatSimObj
#'
#' @param x of class `SpatSimObj`
#' @param ... other things to pass to the plot method for SpatSimObj including `nrow`, `ncol` for the number of rows and columns of plots,
#' `which` processes to plot, and `what` which currently only works with "Processes" but may be updated in the future
#'
#' @method plot SpatSimObj
#'
#'
#' @export
plot.SpatSimObj <- function(x, ...){ #
  if(!exists("what")) what = "Patterns"
  if(methods::is(nrow, "function")) nrow = 1
  if(methods::is(ncol, "function")) ncol = 1
  if(methods::is(which, "function")) which = 1

  if(what == "Patterns"){
    if(length(which) == 1){
      plot(x@Patterns[[which]], main = which)

    } else {
      graphics::par(mfrow = c(nrow, ncol))
      hld = lapply(which, function(i){
        plot(x@Patterns[[i]], main = i)
      })
    }
  } else {
    cat("please be patient while method is being updated")
  }
}

setMethod(f = "show",
          signature = "SpatSimObj",
          definition = function(object){
            summary(object)
          })
