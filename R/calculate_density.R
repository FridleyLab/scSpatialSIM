#' Compute Simulation Heatmaps
#'
#' @param sim_object object created with `CreateSimulationObject`
#' @param steps which simulation steps to compute heatmaps for (Tissue, Holes, Cells, or All)
#' @param which which simulation to compute it for
#' @param step_size resolution of heatmap
#' @param cores number of cpu cores
#'
#' @return a new `SpatSimObj` with probability densities calculated
#' @export
#'
CalculateDensity = function(sim_object, steps = NULL, which = "all", step_size = 1, cores = 1){
  if(!methods::is(sim_object, "SpatSimObj")) stop("`sim_object` must be of class 'SpatSimObj'")
  #slots with parameters
  slots = c('Tissue', 'Holes', 'Cells')
  #make the input step names case-correct
  step_adj = sapply(steps, function(x){
    paste0(toupper(substr(x, 1, 1)),
           tolower(substr(x, 2, nchar(x))))
  })

  #check for whether user input appropriate slot names
  #for all, just use the slot names to pull from
  if("All" %in% step_adj){
    step_adj = c('Tissue', 'Holes', 'Cells')
  } else {
    #otherwise get the names submitted
    step_adj = intersect(step_adj, slots)
  }
  #stop if names submitted don't work
  if(length(step_adj) == 0) stop("Provide a step from which to extract parameters: `Tissue`, `Holes`, `Cells`, or `All`")

  #environment variables
  win_limits = c(sim_object@Window$xrange, sim_object@Window$yrange)
  #make gric based on user step size for their window
  grid = expand.grid(x = seq(win_limits[1], win_limits[2], step_size),
                     y = seq(win_limits[3], win_limits[4], step_size))
  #set the cell IDs for which to generate grids
  if(which == "all") which = seq(sim_object@Sims)

  #looping over steps
  names(step_adj) = step_adj

  #run grids
  out = parallel::mclapply(step_adj, function(step){
    #if not working with the cell kernels
    if(step != "Cells"){
      message(step)
      #make empty class
      if(step == "Tissue"){
        cl = new("Tissue1/Tissue2")
      } else {
        cl = new("Holes")
      }
      #get the slot data
      s = methods::slot(sim_object, step)
      if(length(s@`Simulated Kernels`) == 0 ){
        #let user know that these haven't been used for any simulation step yet
        message(paste0("\t", step, " has not yet been simulated"))
        next
      }
      #make new
      cl@Parameters = s@Parameters
      cl@`Simulated Kernels` = s@`Simulated Kernels`
      cl@`Density Grids` = lapply(which, function(w){
        cat(w, "\n")
        cbind(grid,
              prob = CalculateGrid(grid, s@`Simulated Kernels`[[w]], cores = cores))
      })

      #assign slot back to bid data
      if(step == "Tissue"){
        sim_object@Tissue <<- cl
      } else {
        sim_object@Holes <<- cl
      }
    }

    if(step == "Cells"){
      message(step)
      for(cell in seq(sim_object@Cells)){
        #make new cell class object
        cl = new("Cell")
        s = methods::slot(sim_object, step)[[cell]]
        #fill in with existing parameters
        cl@Parameters = s@Parameters
        cl@`Simulated Kernels` = s@`Simulated Kernels`
        cl@`Density Grids` = lapply(which, function(w){
          cat(w, "\n")
          cbind(grid,
                prob = CalculateGrid(grid, s@`Simulated Kernels`[[w]], cores = cores))
        })

        sim_object@Cells[[cell]] <<- cl
      }
    }
  })

  return(sim_object)
}
