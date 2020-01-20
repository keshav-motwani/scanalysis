#' Plot feature distributions from SingleCellExperiment objects
#'
#' @param sce_list List of SingleCellExperiment objects
#' @param assay Assay to obtain data from (counts, logcounts, etc)
#' @param alt_exp Alternate experiment to obtain data from
#' @param features Numeric features to plot - can be from reducedDims, colData, or assay data
#' @param x colData variable to plot on x-axis
#' @param color colData variable to color points by
#' @param facet_rows Variables from colData to facet on, can also include ".sample" or ".feature" as described below
#' @param facet_columns Variables from colData to facet on, can also include ".sample" or ".feature" as described below
#' @param ... Parameters to ggexp::plot_distributions
#'
#' @importFrom dplyr mutate
#' @importFrom purrr imap_dfr
#' @importFrom tidyr pivot_longer
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' NULL
plot_distributions = function(sce_list,
                              assay,
                              alt_exp = NULL,
                              features,
                              x,
                              color = x,
                              facet_rows = c(".feature"),
                              facet_columns = c(".sample"),
                              ...) {
  if (is.null(names(sce_list))) {
    names(sce_list) = paste0("sample_", 1:length(sce_list))
  }

  original_features = features
  features = c(x, color, features, facet_rows, facet_columns)

  id_cols = unique(c("barcode", x, color, facet_rows, facet_columns, ".sample"))

  data = purrr::imap_dfr(sce_list,
                         ~ get_cell_features(.x, features, assay, alt_exp) %>%
                           mutate(., .sample = .y))

  id_cols = intersect(id_cols, colnames(data))

  data = tidyr::pivot_longer(data,
                             cols = -id_cols,
                             names_to = ".feature",
                             values_to = "value")

  levels = intersect(original_features, data$.feature)

  data$.feature = factor(data$.feature, levels = levels)

  plot = ggexp::plot_distributions(
    data = data,
    x = x,
    y = "value",
    color = color,
    facet_rows = facet_rows,
    facet_columns = facet_columns,
    ...
  ) + labs(y = NULL)

  return(plot)
}
