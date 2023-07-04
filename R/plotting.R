#' Plot Simulation
#'
#' Plot different aspects of a SpatSimObj
#'
#' @param sim_object A `SpatSimObj`
#' @param nrow Number of rows of plots (only applicable when more than one plot is made)
#' @param ncol Number of columns of plots (only applicable when more than one plot is made)
#' @param which Index of the elements of the SpatSimObj to be plotted
#' @param what What to plot ("tissue heatmap", "hole heatmap", or "whole core")
#'
#' @details
#' The `PlotSimulation` function is used to plot different aspects of a SpatSimObj
#' The function takes a `sim_object` as its first argument, which should be an object of class
#' "Spatial Simulation Object". The function can then be used to plot different aspects of the
#' simulation, such as heatmaps of the tissue or holes, or a plot of the whole core with assigned cells colored by type.
#'
#' When `what` is set to "tissue heatmap" or "hole heatmap", the function will plot heatmaps of
#' the specified tissue or hole. When `what` is set to "whole core", the function will plot the
#' entire core with assigned cells colored by type. Only a single element of the `sim_object` can
#' be plotted when `what` is set to "whole core".
#'
#' When more than one plot is made, `nrow` and `ncol` can be used to specify the number of rows
#' and columns of the plot grid, respectively.
#'
#' @return A plot or a grid of plots, depending on the input arguments
#'
#' @export
#'
#' @examples
#' # create a SpatSimObj
#' sim_object <- CreateSimulationObject()
#' sim_object = GenerateSpatialProcess(sim_object)
#' sim_object = GenerateTissue(sim_object, density_heatmap = TRUE, step_size = 1, cores = 1)
#' # plot a heatmap of tissue 1
#' PlotSimulation(sim_object, which = 1, what = "tissue heatmap")
PlotSimulation = function(sim_object, nrow = 1, ncol = 1, which = 1, what = "tissue heatmap"){
  #require(ggplot2)
  if(!methods::is(sim_object, "SpatSimObj")) stop("`sim_object` must be of class 'SpatSimObj'")
  density_plot = function(dat){
    dat %>%
      ggplot2::ggplot() +
      ggplot2::geom_contour_filled(ggplot2::aes(x = x, y = y, z = 100*prob)) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::guides(fill=ggplot2::guide_legend(title="Probability\nof Stroma"))
  }
  if(what == "tissue heatmap"){
    if(length(which) == 1){
      density_plot(sim_object@Tissue@`Density Grids`[[which]]) +
        ggplot2::labs(title = which)
    } else {
      p = lapply(seq(which), function(d){
        density_plot(sim_object@Tissue@`Density Grids`[[d]]) +
          ggplot2::labs(title = d)
      })
      ggpubr::ggarrange(plotlist = p, ncol = ncol, nrow = nrow, common.legend = TRUE)
    }
  } else if(what == "hole heatmap"){
    if(length(sim_object@Holes@`Simulationed Kernels`) == 0) stop("Holes have not been simulated")
    if(length(which) == 1){
      density_plot(sim_object@Holes@`Density Grids`[[which]]) +
        ggplot2::labs(title = which)
    } else {
      p = lapply(seq(which), function(d){
        density_plot(sim_object@Holes@`Density Grids`[[d]]) +
          ggplot2::labs(title = d)
      })
      ggpubr::ggarrange(plotlist = p, ncol = ncol, nrow = nrow, common.legend = TRUE)
    }
  } else if(what == "whole core"){
    if(length(which) == 1){
      if(length(sim_object@Cells[[1]]@`Simulationed Kernels`) == 0){
        stop("Need to simulate cells")
      }
      df = sim_object@`Spatial Files`[[which]]
      cell_cols = grep("Cell", colnames(df), value = T) %>%
        grep("Assignment", ., value = T)
      df = df %>%
        dplyr::mutate_at(cell_cols, ~ ifelse(.x == "Positive", 1, 0))
      if(length(cell_cols) > 1){
        df = df %>%
          dplyr::mutate(Multi_Hit = ifelse(rowSums(.[cell_cols]) > 1, 1, 0)) %>%
          dplyr::mutate_at(cell_cols, ~ ifelse(Multi_Hit == 1, 0, .x)) %>%
          dplyr::mutate(Background = ifelse(rowSums(.[c(cell_cols, "Multi_Hit")]) > 0, 0, 1)) %>%
          dplyr::filter(dplyr::if_any(dplyr::matches("Hole Assignment"), ~.x == "Keep")) %>% #{if("Hole Assignment" %in% names(.)) `Hole Assignment` else NULL} == "Keep"
          dplyr::select(x, y, `Tissue Assignment`, !!cell_cols, Multi_Hit, Background)
      } else {
        df = df %>%
          dplyr::mutate(Background = ifelse(rowSums(.[c(cell_cols)]) > 0, 0, 1)) %>%
          dplyr::filter(dplyr::if_any(dplyr::matches("Hole Assignment"), ~.x == "Keep")) %>%
          dplyr::select(x, y, `Tissue Assignment`, !!cell_cols, Background)
      }
      df = df %>%
        tidyr::gather("Cell Type", "Positive", -x, -y, -`Tissue Assignment`) %>%
        dplyr::filter(Positive == 1)
      ggplot2::ggplot() +
        ggplot2::geom_point(data = df %>% dplyr::filter(`Cell Type` == 'Background'),
                            ggplot2::aes(x = x, y = y, shape = `Tissue Assignment`), color = 'gray') +
        ggplot2::geom_point(data = df %>% dplyr::filter(`Cell Type` != 'Background'),
                            ggplot2::aes(x = x, y = y, color = `Cell Type`, shape = `Tissue Assignment`)) +
        ggplot2::theme_bw() +
        tune::coord_obs_pred() +
        ggplot2::scale_shape_manual(values=c(3, 16)) +
        ggplot2::facet_wrap(~`Tissue Assignment`) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title = which)
    } else {
      stop("For whole core, only a single spatial file can be plotted at once")
    }
  } else {
    cat("please be patient while method is being updated")
  }
}
