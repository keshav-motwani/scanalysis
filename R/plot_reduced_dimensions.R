#' Plot reduced dimensional plot with multiple features
#'
#' Features are min-max normalized per feature, and the range of each feature is annotated per facet to consolidate multiple features into one color scale.
#'
#' If multiple SingleCellExperiments are provided in the sce_list, and you want to facet by this, you can add ".sample" to one of the faceting variables, as this is implicitly added in.
#'
#' In almost all cases, you would want to facet by feature. This is by default included in facet_columns, but if you add additional faceting variables to facet_columns, be sure to also include ".feature". To be clear, ".feature" can also be included in facet_rows. This is not included within the code to allow for flexibility if you want features to be faceted in the rows or columns.
#'
#' @param sce_list List of SingleCellExperiment objects to plot
#' @param type Name of reducedDim attribute to plot
#' @param assay Assay to obtain data from (ex: counts, logcounts)
#' @param alt_exp Alternate experiment to obtain data from
#' @param features Features to plot - can be from reducedDims, colData, or assay data, but note that all must be either numeric or categorical for one plot
#' @param alpha Alpha for points
#' @param pt_size Size of points
#' @param facet_rows Variables from colData to facet on, can also include ".sample" or ".feature" as described below
#' @param facet_columns Variables from colData to facet on, can also include ".sample" or ".feature" as described below
#' @param facet_type Either "wrap" or "grid", same as ggplot
#' @param facet_scales Either "fixed", "free", "free_x", "free_y", same as ggplot
#' @param facet_switch Either NULL, "x", "y", "both", same as ggplot
#' @param nrow Number of rows if facet_type is "wrap"
#'
#' @import ggplot2
#' @importFrom dplyr mutate bind_cols group_by summarize arrange
#' @importFrom purrr imap_dfr
#' @importFrom tidyr pivot_longer
#' @importFrom ggexp plot_facets theme_ggexp
#' @importFrom gtools mixedsort
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' sce = scater::mockSCE() %>% scater::logNormCounts() %>% scater::runPCA()
#' plot_reduced_dimensions(sce_list = list(sample_1 = sce, sample_2 = sce),
#'                         features = c("Gene_0001", "Gene_0002", "Gene_0003"),
#'                         facet_columns = ".sample",
#'                         facet_rows = ".feature",
#'                         facet_switch = "y")
plot_reduced_dimensions = function(sce_list,
                                   type = "PCA",
                                   assay = "logcounts",
                                   alt_exp = NULL,
                                   features,
                                   alpha = 0.3,
                                   pt_size = 1,
                                   facet_rows = c(),
                                   facet_columns = c(".feature"),
                                   facet_type = "grid",
                                   facet_scales = "free",
                                   facet_switch = NULL,
                                   nrow = 2) {
  if (is.null(names(sce_list)))
    names(sce_list) = paste0("sample_", 1:length(sce_list))

  data = imap_dfr(
    sce_list,
    ~ get_cell_features(.x, c(features, facet_rows, facet_columns), assay, alt_exp) %>%
      mutate(., .sample = .y) %>%
      bind_cols(get_reduced_dims(.x, type))
  ) %>%
    pivot_longer(
      cols = intersect(features, colnames(.)),
      names_to = ".feature",
      values_to = "value"
    )

  data$.feature = factor(data$.feature, levels = features)

  if (is.numeric(data$value)) {
    min_max = data %>%
      group_by(.dots = c(facet_rows, facet_columns)) %>%
      summarize(min = min(value, na.rm = TRUE),
                max = max(value, na.rm = TRUE)) %>%
      mutate(value = paste0(round(min, 2), "-", round(max, 2)))

    data = data %>%
      group_by(.dots = c(".feature")) %>%
      mutate(value = (value - min(value, na.rm = TRUE)) / (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))) %>%
      arrange(value)

  } else {
    data$value = factor(as.character(data$value),
                        levels = c("NA", mixedsort(as.character(
                          unique(data$value[data$value != "NA"])
                        ), na.last = FALSE)))
    data = arrange(data, value)
  }

  plot = ggplot(data, aes_string(
    x = paste0((type), "_1"),
    y = paste0((type), "_2"),
    color = "value"
  )) +
    geom_point(alpha = alpha,
               size = pt_size) +
    theme_ggexp() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  if (is.numeric(data$value))
    plot = plot + geom_text(
      data = min_max,
      aes(x = Inf, y = Inf, label = value),
      hjust = 1,
      vjust = 1.2,
      size = 3,
      inherit.aes = FALSE
    )

  plot = plot_facets(plot,
                     facet_rows,
                     facet_columns,
                     facet_type,
                     facet_scales,
                     facet_switch,
                     nrow)

  if (!is.numeric(data$value)) {
    plot = plot +
      guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5)))
  } else {
    plot = plot + scale_color_gradient(
      low = "#E8E8E8",
      high = "firebrick",
      breaks = c(0, 1),
      labels = c("min", "max"),
      limits = c(0, 1)
    )
  }

  return(plot + theme(legend.title = element_blank()))
}

get_reduced_dims = function(sce, type) {
  reduced_dims = data.frame(
    dim1 = SingleCellExperiment::reducedDims(sce)@listData[[type]][, 1],
    dim2 = SingleCellExperiment::reducedDims(sce)@listData[[type]][, 2]
  )

  colnames(reduced_dims) = paste0(type, "_", c(1, 2))

  reduced_dims$barcode = SummarizedExperiment::colData(sce)$Barcode

  return(reduced_dims)
}