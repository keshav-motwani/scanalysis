#' Plot feature distributions from SingleCellExperiment objects
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param features numeric features to plot - can be from reducedDims, colData, or assay data
#' @param x colData variable to plot on x-axis
#' @param color colData variable to color points by
#' @param fill colData variable to fill by
#' @param facet_rows variables from colData to facet on, can also include ".sample" or ".feature" as described below
#' @param facet_columns variables from colData to facet on, can also include ".sample" or ".feature" as described below
#' @param assay assay to obtain data from (counts, logcounts, etc)
#' @param alt_exp alternate experiment to obtain data from
#' @param ... other parameters passed to ggexp::plot_distributions
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
plot_feature_distributions = function(sce_list,
                         features,
                         x,
                         color = x,
                         fill = x,
                         facet_rows = c(".feature"),
                         facet_columns = c(".sample"),
                         assay = "logcounts",
                         alt_exp = NULL,
                         ...) {
  if (is.null(names(sce_list))) {
    names(sce_list) = paste0("sample_", 1:length(sce_list))
  }

  original_features = features
  features = c(x, color, fill, features, facet_rows, facet_columns)

  id_cols = unique(c("barcode", x, color, facet_rows, facet_columns, ".sample"))

  data = purrr::imap_dfr(
    sce_list,
    ~ get_cell_features(.x, features, assay, alt_exp) %>%
      mutate(., .sample = .y)
  )

  id_cols = intersect(id_cols, colnames(data))

  data = tidyr::pivot_longer(
    data,
    cols = -id_cols,
    names_to = ".feature",
    values_to = "value"
  )

  levels = intersect(original_features, data$.feature)

  data$.feature = factor(data$.feature, levels = levels)

  plot = ggexp::plot_distributions(
    data = data,
    x = x,
    y = "value",
    color = color,
    fill = fill,
    facet_rows = facet_rows,
    facet_columns = facet_columns,
    ...
  ) + labs(y = NULL)

  return(plot)
}

#' @export
#' @rdname plot_feature_distributions
plot_features = plot_feature_distributions