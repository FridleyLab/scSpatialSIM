#' Get Spatial Files from a `SpatSimObj`
#'
#' This function extracts the 'Spatial Files' slot from a Spatial Simulation Object
#' and removes probability columns while converting 'Positive' and 'Negative' in
#' cell assignment columns to 1 and 0, respectively.
#'
#' The output of this function creates a list of the spatial files formatted in a way
#' that would allow direct import into a mIF object from the package 'spatialTIME'
#'
#' @param sim_object A `SpatSimObj`
#' @param single_df boolean as to whether to collapse the output list of data frames to a single data frame or not. default is `FALSE`
#' @param multihit_action string of value 'random', 'drop', or 'keep' for cells that are positive for multiple cell assignments
#' @return A list of data frames, one for each simulated cell type, with cleaned columns
#'
#' @export
CreateSpatialList =  function(sim_object, single_df = FALSE, multihit_action = "random"){
  #require(dplyr)
  if(!methods::is(sim_object, "SpatSimObj")) stop("`sim_object` must be of class 'SpatSimObj'")

  if(is.null(sim_object@`Spatial Files`)){
    stop("There are no 'Spatial Files' associated with this 'SpatSimObj'")
  }

  if(length(multihit_action) != 1 &
     multihit_action %notin% c('random', 'drop', 'keep'))
    stop("`multihit_action` should be appropriate single element string")

  # Make a copy of the `Spatial Files` slot to avoid changing the original object
  spatial_files <- sim_object@`Spatial Files`
  #get cell numbers
  assignment_cols <- grep("Cell", colnames(sim_object@`Spatial Files`[[1]]), value = T) %>%
    grep("Prob", ., value = T, invert = T)

  spatial_files = pbmcapply::pbmclapply(seq_along(spatial_files), function(i){
    # Remove probability columns
    tmp_dat <- dplyr::select(spatial_files[[i]], -dplyr::contains("Probability"))

    # Replace "Positive" with 1 and "Negative" with 0 in assignment columns
    tmp_dat[, assignment_cols] <-
      lapply(tmp_dat[, assignment_cols],function(x) ifelse(x == "Positive", 1, 0)) %>%
      unlist()
    #return data
    return(tmp_dat)
  }, mc.cores = 1)
  #give each spatial data frame a psuedo sample name
  names(spatial_files) = paste("Spatial Data", seq_along(spatial_files))

  #for multi-hit cells
  if(length(assignment_cols) > 1){
    if(multihit_action == "random"){ #randomly select positive
      spatial_files = lapply(spatial_files, multihit_random, assignment_cols)
    } else if(multihit_action == "drop"){ #drop mutlihit cells
      spatial_files = lapply(spatial_files, multihit_drop, assignment_cols)
    }
    #if multihit_action == 'keep' then do nothing because already that way

  }

  #for functional data analysis collapse to single data frame with the names in an ID column
  if(single_df){
    #collapse all down to single dataframe
    all_df = lapply(names(spatial_files), function(spat_name){
      spatial_files[[spat_name]] %>%
        dplyr::mutate(`Image Name` = spat_name, .before = 1)
    }) %>%
      do.call(dplyr::bind_rows, .)
    #return the single dataframe to the user
    return(all_df)
  }

  return(spatial_files)
}
