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
#' be plotted when `what` is set to "whole core". `what` equal to "tissue points", "hole points", or
#' "tissue hole points" will result in a point plot of the respective assignments on points.
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
#' sim_object = GenerateSpatialPattern(sim_object)
#' sim_object = GenerateTissue(sim_object, density_heatmap = TRUE, step_size = 1, cores = 1)
#' # plot a heatmap of tissue 1
#' PlotSimulation(sim_object, which = 1, what = "tissue heatmap")
PlotSimulation = function(sim_object, nrow = 1, ncol = 1, which = 1, what = "tissue heatmap"){
  #require(ggplot2)
  if(!methods::is(sim_object, "SpatSimObj")) stop("`sim_object` must be of class 'SpatSimObj'")

  if(what == "tissue heatmap" &
     sum(!(which %in% seq(sim_object@Tissue@`Density Grids`)))>0) stop("`which` not in Tissue Density Grids")
  if(what == "tissue points" &
     sum(!(which %in% seq(sim_object@`Spatial Files`))>0) &
     (TRUE %in% grepl("Tissue", colnames(sim_object@`Spatial Files`[[1]])))) stop("`which` not in Spatial Files or Tissue not simulated")

  if(what == "hole heatmap" &
     sum(!(which %in% seq(sim_object@Tissue@`Density Grids`)))>0) stop("`which` not in Tissue Density Grids")
  if(what == "hole points" &
     sum(!(which %in% seq(sim_object@`Spatial Files`))>0) &
     (TRUE %in% grepl("Hole", colnames(sim_object@`Spatial Files`[[1]])))) stop("`which` not in Spatial Files or Holes not simulated")

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
  } else if(what == "tissue points"){
    if(length(which) == 1){
      tissue_point_plot(sim_object@`Spatial Files`[[which]]) +
        ggplot2::labs(title = which)
    } else {
      p = lapply(seq(which), function(d){
        tissue_point_plot(sim_object@`Spatial Files`[[d]]) +
          ggplot2::labs(title = d)
      })
      ggpubr::ggarrange(plotlist = p, ncol = ncol, nrow = nrow, common.legend = TRUE)
    }
  } else if(what == "hole points"){
    if(length(which) == 1){
      hole_point_plot(sim_object@`Spatial Files`[[which]]) +
        ggplot2::labs(title = which)
    } else {
      p = lapply(seq(which), function(d){
        hole_point_plot(sim_object@`Spatial Files`[[d]]) +
          ggplot2::labs(title = d)
      })
      ggpubr::ggarrange(plotlist = p, ncol = ncol, nrow = nrow, common.legend = TRUE)
    }
  } else if(what == "tissue hole points"){
    if(length(which) == 1){
      tissue_hole_point_plot(sim_object@`Spatial Files`[[which]]) +
        ggplot2::labs(title = which)
    } else {
      p = lapply(seq(which), function(d){
        tissue_hole_point_plot(sim_object@`Spatial Files`[[d]]) +
          ggplot2::labs(title = d)
      })
      ggpubr::ggarrange(plotlist = p, ncol = ncol, nrow = nrow, common.legend = TRUE)
    }
  } else if(what == "whole core"){
    if(length(which) == 1){
      if(length(sim_object@`Spatial Files`[[1]]) == 0){
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
        ggplot2::coord_equal() +
        ggplot2::scale_shape_manual(values=c(3, 16)) +
        ggplot2::facet_wrap(~`Tissue Assignment`) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(title = which)
    } else {
      stop("For whole core, only a single spatial file can be plotted at once")
    }
  } else {
    stop("Please be patient while method is being updated")
  }
}

#helpers
density_plot = function(dat){
  dat %>%
    ggplot2::ggplot() +
    ggplot2::geom_contour_filled(ggplot2::aes(x = x, y = y, z = 100*prob)) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::guides(fill=ggplot2::guide_legend(title="Probability\nSurface"))
}
tissue_point_plot = function(dat){
  dat %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y, shape = `Tissue Assignment`)) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}
hole_point_plot = function(dat){
  dat %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y, alpha = `Hole Assignment`)) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}
tissue_hole_point_plot = function(dat){
  dat %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y, shape = `Tissue Assignment`, alpha = `Hole Assignment`)) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}
