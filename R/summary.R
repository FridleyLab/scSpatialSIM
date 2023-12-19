#' summary function for SpatSimObj
#'
#' @param object of class `SpatSimObj`
#' @param ... nothing else to pass to summary if object is a SpatSimObj
#'
#' @method summary SpatSimObj
#' @returns summary of the SpatSimObj to the terminal
#'
#' @export
summary.SpatSimObj <- function(object, ...){
  nsim = object@Sims
  win = matrix(c(object@Window$xrange, object@Window$yrange), nrow = 2, byrow = FALSE)
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
#' `which` pattern to plot, and `what` which currently only works with "Processes" but may be updated in the future
#'
#' @method plot SpatSimObj
#' @returns basic x-y ggplot object
#'
#' @export
plot.SpatSimObj <- function(x, ...){ #
  params = as.list(substitute(list(...)))
  if(!("ncol" %in% names(params))) params$ncol = 1
  if(!("nrow" %in% names(params))) params$ncol = 1
  if(!("which" %in% names(params))) params$ncol = 1
  params$what = "Patterns"

  basic_plot = function(x, p){
    x@Patterns[[p]] %>%
      ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x = x, y = y)) +
      ggplot2::labs(title = p) +
      ggplot2::coord_equal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }

  if(params$what == "Patterns"){#
    if(length(params$which) == 1){
      basic_plot(x, 1)
    } else {
      hld = lapply(eval(params$which), function(i){
        #print(i)
        basic_plot(x, i)
      })
      ggpubr::ggarrange(plotlist = hld, ncol = params$ncol, nrow = params$nrow)
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
