#' Summarise Spatial
#'
#' @param spatial_list list of spatial data frames with `markers` column names
#' @param markers names of columns, probably cell types, that contain 1s and 0s representing positive/negative assignments
#'
#' @return data frome with summary counts and proportions for the markers in each spatial data frame
#' @export
#'
#'
SummariseSpatial = function(spatial_list, markers){
  #find cell frequency for each data frame in the spatial list
  out = lapply(seq(spatial_list),  function(spat_num){
    spat = spatial_list[[spat_num]] %>%
      #find total number of rows in dataframe
      dplyr::mutate(`Total Cells` = dplyr::n()) %>%
      #group by to maintain total cell count
      dplyr::group_by(`Total Cells`) %>%
      #find number of cells for each marker
      dplyr::summarise(dplyr::across(!!markers, ~sum(.x))) %>%
      #calculate the proportion of the different markers in the spatial data
      dplyr::mutate(dplyr::across(!!markers, ~ .x/`Total Cells` * 100, .names = "% {.col}")) %>%
      #add a column for which spatial data frame it came from
      dplyr::mutate(`Sample Tag` = paste("Spatial Data", spat_num), .before = 1)
    return(spat)
  }) %>%
    #bind the new data together
    do.call(dplyr::bind_rows, .)
  #return data
  return(out)
}


