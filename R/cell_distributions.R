#check if cells have been simulated
#' Generate Characteristic Distributions of Cells
#'
#' @param spatial_data object of either class `list` or `data.frame`. Can be created from `SpatSimObj` with \code{\link{CreateSpatialList}}
#' @param positive_mean,negative_mean number for mean of which to center the distribution of the positive and negative cell types.
#' Can be single number of vector with length matching number of Cells.
#' @param positive_sd,negative_sd number for the standard deviation around the positive cell type mean. Can be single value or same length as number of Cells.
#'
#' @return object with a class that is the same as input `spatial_data` with new columns containing distributions for positive/negative assigned cells
#' @export
#'
#' @examples
#' #create simulation object
#' spatial_data = CreateSimulationObject(sims = 1, cell_types = 1) %>%
#'   #produce the point pattern
#'   GenerateSpatialPattern() %>%
#'   #make tissues
#'   GenerateTissue(density_heatmap = FALSE, step_size = 0.1, cores = 1) %>%
#'   #create positive and negative cells
#'   GenerateCellPositivity(k = 4, sdmin = 3, sdmax = 5,
#'   density_heatmap = FALSE, step_size = 1, cores = 1, probs = c(0.0, 0.1), shift = 0) %>%
#'   #convert to a list of spatial data frames
#'   CreateSpatialList(single_df = FALSE)
#' spat_data_distribution = GenerateDistributions(spatial_data)
GenerateDistributions = function(spatial_data,
                                 positive_mean = 10,
                                 negative_mean = 2,
                                 positive_sd = 2,
                                 negative_sd = 1){
  #make sure object is of right class
  if(!(methods::is(spatial_data, "list") | methods::is(spatial_data, "data.frame"))) stop("Object must be data frame or list of data frames.")
  #get class of the spatial data
  dat_class = class(spatial_data)
  #number of cells
  if(dat_class == "list"){
    cells = length(grep("Cell", colnames(spatial_data[[1]])))
  } else {
    cells = length(grep("Cell", colnames(spatial_data)))
  }

  #check variable lengths
  if(FALSE %in% ((sapply(list(positive_mean, positive_sd, negative_mean, negative_sd), length) == cells) |
    (sapply(list(positive_mean, positive_sd, negative_mean, negative_sd), length) == 1)))
    stop("Values must either be of length 1 or same as number of cells")

  #if list iterate over all data frames
  if(dat_class == "list"){
    spat_list = pbmcapply::pbmclapply(spatial_data, function(spat){
      make_dis(spat, positive_mean, negative_mean, positive_sd, negative_sd)
    }, mc.cores = 1)
    return(spat_list)
  } else {
    #if data frame, perform once
    spat_df = make_dis(spatial_data, positive_mean, negative_mean, positive_sd, negative_sd)
    return(spat_df)
  }
}
