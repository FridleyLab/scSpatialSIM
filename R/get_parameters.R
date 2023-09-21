#' Get Simulation Parameters
#'
#' @param sim_object simulation object created with \code{\link{CreateSimulationObject}}
#' @param steps which parameters to extract from the `SpatSimObj`
#'
#' @return a list of S3 clas `SimParams` containing parameters stored in the `sim_object`
#' @export
#'
#' @details
#' This function will return any paramters that are stored in the simulation object. If no simulation steps
#' have been run, then this will return the default parameters. The defaults are over written if new parameters
#' are provided at each step.
#'
#'
#' @examples
#'
#' #create simulation object
#' sim_obj = CreateSimulationObject()
#' #extract default paramters for the Tissue simulation step
#' defs = ExtractParameters(sim_obj, "Tissue")
ExtractParameters = function(sim_object, steps = NULL){
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

  #looping over steps
  names(step_adj) = step_adj
  out = lapply(step_adj, function(step){
    #if lookings at cells, need to loop one deeper
    if(step != "Cells"){
      s = methods::slot(sim_object, step)
      if(length(s@`Simulated Kernels`) == 0 )
        #let user know that these haven't been used for any simulation step yet
        message(paste0("\t", step, " has not yet been simulated"))
      p = s@Parameters
      #return the parameters
      return(p)
    }

    #if cells then we need to dig further
    if(step == "Cells"){
      cells = lapply(sim_object@Cells, function(s){
        if(length(s@`Simulated Kernels`) == 0 )
          #let user know that these haven't been used for any simulation step yet
          message("\tCells have not yet been simulated")
        p = s@Parameters
        #return the parameters
        return(p)
      })

      #give cells different names
      names(cells) = paste0("Cell Type ", seq(cells))

      #return cell paramter list
      return(cells)
    }
  })

  class(out) = "SimParams"
  return(out)
}
