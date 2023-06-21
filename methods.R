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
