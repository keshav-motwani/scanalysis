#' Plot pairwise scatterplot of cell-level data
#'
#' @param sce_list List of SingleCellExperiment objects to plot
#' @param assay Assay to obtain data from (ex: counts, logcounts)
#' @param alt_exp Alternate experiment to obtain data from
#' @param x Numeric features to plot on x axis - can be from reducedDims, colData, or assay data
#' @param y Numeric features to plot on y axis - can be from reducedDims, colData, or assay data
#' @param color Column from reducedDims, colData, or assay data to color by
#' @param shape Column from reducedDims, colData, or assay data to shape by
#' @param size Column from reducedDims, colData, or assay data to size by
#' @param ... Other parameters to be passed to ggexp::plot_pairwise_scatterplot
#'
#' @importFrom ggexp plot_pairwise_scatterplot
#' @importFrom purrr imap_dfr
#' @importFrom dplyr mutate
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' NULL
plot_pairwise_scatterplot = function(sce_list,
                                     assay = "logcounts",
                                     alt_exp = NULL,
                                     x,
                                     y = x,
                                     color = NULL,
                                     shape = NULL,
                                     size = NULL,
                                     facet_rows = c(),
                                     facet_columns = c(),
                                     ...) {
  if (is.null(names(sce_list)))
    names(sce_list) = paste0("sample_", 1:length(sce_list))

  data = imap_dfr(sce_list,
                  ~ get_cell_features(
                    .x,
                    c(x, y, color, size, facet_rows, facet_columns),
                    assay,
                    alt_exp
                  ) %>%
                    mutate(., .sample = .y))

  plot = plot_pairwise_scatterplot(
    data = data,
    x = intersect(x, colnames(data)),
    y = intersect(y, colnames(data)),
    color = color,
    shape = shape,
    size = size,
    facet_rows = facet_rows,
    facet_columns = facet_columns,
    ...
  )

  return(plot)
}
