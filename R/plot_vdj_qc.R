#' Scatterplot of two features of interest from colData with annotated thresholds and counts based on filters
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param x_filters list of filters for each SCE object in sce_list from scater::isOutlier, or a vector with attribute named thresholds that is a vector with min and max allowed values
#' @param x numeric column from colData that is in all objects in sce_list
#' @param y discrete column from colData that is in all objects in sce_list to split histograms by
#' @param color column from colData that is in all objects in sce_list
#' @param shape column from colData that is in all objects in sce_list
#' @param x_log whether to use log x-axis
#' @param text_size font size for annotations
#' @param facet_rows columns to facet on
#' @param facet_columns columns to facet on
#' @param facet_type either "wrap" or "grid", same as ggplot
#' @param ... other params passed into either facet_wrap or facet_grid, depending on facet_type parameter
#'
#' @import ggplot2
#' @importFrom purrr map map2
#' @importFrom dplyr summarise select group_by arrange desc
#'
#' @return
#' @export
#'
#' @examples
#' NULL
plot_vdj_gex_univariate_qc = function(sce_list,
                                      x_filters = NULL,
                                      vdj_filters = NULL,
                                      x,
                                      ...) {
  plot = plot_gex_univariate_qc(sce_list,
                                x_filters = x_filters,
                                x = x,
                                y = "count_TRA_TRB_IGL_IGK_IGH",
                                ...)

  order = attributes(vdj_filters[[1]])$allowed_values %>%
    arrange(desc(allowed), count_TRA_TRB_IGL_IGK_IGH)

  plot$data$count_TRA_TRB_IGL_IGK_IGH = factor(plot$data$count_TRA_TRB_IGL_IGK_IGH,
                                               levels = order$count_TRA_TRB_IGL_IGK_IGH)

  plot = plot + geom_hline(yintercept = sum(order$allowed) + 0.90,
                           linetype = "dashed")

  return(plot)
}