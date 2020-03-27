#' Plot pairwise scatterplot of cell-level data
#'
#' @param sce_list list of SingleCellExperiment objects to plot
#' @param x numeric features to plot on x axis - can be from reducedDims, colData, or assay data
#' @param y numeric features to plot on y axis - can be from reducedDims, colData, or assay data
#' @param color column from reducedDims, colData, or assay data to color by
#' @param shape column from reducedDims, colData, or assay data to shape by
#' @param size column from reducedDims, colData, or assay data to size by
#' @param facet_rows columns from colData to facet on
#' @param facet_columns columns from colData to facet on
#' @param assay assay to obtain data from (ex: counts, logcounts)
#' @param alt_exp alternate experiment to obtain data from
#' @param ... other parameters to be passed to ggexp::plot_pairwise_scatterplot
#'
#' @importFrom ggexp plot_pairwise_scatterplot
#' @importFrom ggrepel geom_label_repel
#' @importFrom purrr imap_dfr
#' @importFrom dplyr mutate
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' NULL
plot_pairwise_features = function(sce_list,
                                  x,
                                  y = x,
                                  assay = "logcounts",
                                  alt_exp = NULL,
                                  color = NULL,
                                  shape = NULL,
                                  label = NULL,
                                  facet_rows = c(),
                                  facet_columns = c(),
                                  text_size = 3,
                                  ...) {
  if (is.null(names(sce_list)))
    names(sce_list) = paste0("sample_", 1:length(sce_list))

  data = imap_dfr(sce_list,
                  ~ get_cell_features(
                    .x,
                    c(x, y, color, shape, label, facet_rows, facet_columns),
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
    facet_rows = facet_rows,
    facet_columns = facet_columns,
    ...
  )

  if (!is.null(label)) {
    annotations = plot$data %>%
      group_by(.dots = c(facet_rows, facet_columns, label)) %>%
      summarize(x = median(.xvalue),
                y = median(.yvalue))

    plot = plot + geom_label_repel(
      data = annotations,
      aes_string(
        x = "x",
        y = "y",
        label = label,
        color = color,
        inherit.aes = FALSE
      ),
      label.padding = unit(0.1, "lines"),
      alpha = 1,
      fill = "white",
      size = text_size
    )
  }

  return(plot)
}
