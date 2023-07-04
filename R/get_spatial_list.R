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
#' @return A list of data frames, one for each simulated cell type, with cleaned columns
#'
#' @export
CreateSpatialList =  function(sim_object, single_df = FALSE){
  #require(dplyr)
  if(!methods::is(sim_object, "SpatSimObj")) stop("`sim_object` must be of class 'SpatSimObj'")

  if(is.null(sim_object@`Spatial Files`)){
    stop("There are no 'Spatial Files' associated with this 'SpatSimObj'")
  }

  # Make a copy of the `Spatial Files` slot to avoid changing the original object
  spatial_files <- sim_object@`Spatial Files`

  spatial_files = pbmcapply::pbmclapply(seq_along(spatial_files), function(i){
    # Remove probability columns
    tmp_dat <- dplyr::select(spatial_files[[i]], -dplyr::contains("Probability"))

    # Replace "Positive" with 1 and "Negative" with 0 in assignment columns
    assignment_cols <- grep("Cell", colnames(tmp_dat))
    tmp_dat[, assignment_cols] <-
      lapply(tmp_dat[, assignment_cols],function(x) ifelse(x == "Positive", 1, 0)) %>%
      unlist()
    #return data
    return(tmp_dat)
  })
  #give each spatial data frame a psuedo sample name
  names(spatial_files) = paste("Spatial Data", seq_along(spatial_files))

  if(single_df){
    all_df = lapply(names(spatial_files), function(spat_name){
      spatial_files[[spat_name]] %>%
        dplyr::mutate(`Image Name` = spat_name, .before = 1)
    }) %>%
      do.call(dplyr::bind_rows, .)
    return(all_df)
  }

  return(spatial_files)
}
