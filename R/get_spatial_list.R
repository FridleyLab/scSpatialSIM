#' Get Spatial Files from a Spatial Simulation Object
#'
#' This function extracts the 'Spatial Files' slot from a Spatial Simulation Object
#' and removes probability columns while converting 'Positive' and 'Negative' in
#' cell assignment columns to 1 and 0, respectively.
#'
#' The output of this function creates a list of the spatial files formatted in a way
#' that would allow direct import into a mIF object from the package 'spatialTIME'
#'
#' @param sim_object A Spatial Simulation Object
#' @return A list of data frames, one for each simulated cell type, with cleaned columns
#'
#' @export
CreateSpatialList =  function(sim_object){
  #require(dplyr)
  if(!methods::is(sim_object, "SpatialSimulationObject")) stop("`sim_object` must be of class 'SpatialSimulationObject'")

  if(is.null(sim_object@`Spatial Files`)){
    stop("There are no 'Spatial Files' associated with this 'SpatialSimulationObject'")
  }

  # Make a copy of the `Spatial Files` slot to avoid changing the original object
  spatial_files <- sim_object@`Spatial Files`

  # Loop over each data frame in the `spatial_files` list
  for(i in seq_along(spatial_files)){
    # Remove probability columns
    spatial_files[[i]] <- dplyr::select(spatial_files[[i]], -dplyr::contains("Probability"))

    # Replace "Positive" with 1 and "Negative" with 0 in assignment columns
    assignment_cols <- grep("Cell", colnames(spatial_files[[i]]))
    spatial_files[[i]][, assignment_cols] <-
      lapply(spatial_files[[i]][, assignment_cols],
             function(x) ifelse(x == "Positive", 1, 0))
  }

  return(spatial_files)
}
