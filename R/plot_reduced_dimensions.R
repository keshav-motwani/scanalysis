#' Plot reduced dimensional plot with multiple features
#'
#' Features are min-max normalized per feature, and the range of each feature is annotated per facet to consolidate multiple features into one color scale.
#'
#' If multiple SingleCellExperiments are provided in the sce_list, and you want to facet by this, you can add ".sample" to one of the faceting variables, as this is implicitly added into the data frame being plotted.
#'
#' In almost all cases, you would want to facet by feature, so be sure to also include ".feature" in either facet_columns or facet_rows
#'
#' @param sce_list list of SingleCellExperiment objects to plot
#' @param features features to plot - can be from reducedDims, colData, or assay data, but note that all must be either numeric or categorical for one plot
#' @param type name of reducedDim attribute to plot
#' @param label_value boolean to annotate text label for the value - only works if all features are discrete
#' @param alpha alpha for points
#' @param point_size size of points
#' @param text_size size of font for text annotation
#' @param lower_quantile quantile which should be used to determine the lower limit of the color bar
#' @param upper_quantile quantile which should be used to determine the upper limit of the color bar
#' @param min_value minimum feature value, below which to set to this value
#' @param facet_rows variables from colData to facet on, can also include ".sample" or ".feature" as described below
#' @param facet_columns variables from colData to facet on, can also include ".sample" or ".feature" as described below
#' @param facet_type either "wrap" or "grid", same as ggplot
#' @param assay assay to obtain data from (ex: counts, logcounts)
#' @param alt_exp alternate experiment to obtain data from
#' @param ... other params passed into either facet_wrap or facet_grid, depending on facet_type parameter
#'
#' @import ggplot2
#' @importFrom dplyr mutate bind_cols group_by summarize arrange
#' @importFrom purrr imap_dfr
#' @importFrom tidyr pivot_longer
#' @importFrom ggexp plot_facets theme_ggexp
#' @importFrom ggrepel geom_label_repel
#' @importFrom gtools mixedsort
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' library(scanalysis)
#'
#' sce = scater::mockSCE() %>%
#'     scater::logNormCounts() %>%
#'     scater::runPCA()
#'
#' plot_reduced_dimensions(sce_list = list(sample_1 = sce, sample_2 = sce),
#'                         features = c("Gene_0001", "Gene_0002", "Gene_0003"),
#'                         facet_columns = ".sample",
#'                         facet_rows = ".feature",
#'                         switch = "y")
plot_reduced_dimensions = function(sce_list,
                                   type = "PCA",
                                   assay = "logcounts",
                                   alt_exp = NULL,
                                   features,
                                   label = NULL,
                                   alpha = 1,
                                   point_size = 0.05,
                                   text_size = 3,
                                   lower_quantile = 0,
                                   upper_quantile = 1,
                                   min_value = NULL,
                                   facet_rows = c(),
                                   facet_columns = c(),
                                   facet_type = "grid",
                                   ...) {
  if (is.null(names(sce_list)))
    names(sce_list) = paste0("sample_", 1:length(sce_list))

  data = imap_dfr(
    sce_list,
    ~ get_cell_features(
      .x,
      c(features, facet_rows, facet_columns, label),
      assay,
      alt_exp
    ) %>%
      mutate(., .sample = .y) %>%
      bind_cols(.get_reduced_dims(.x, type))
  ) %>%
    pivot_longer(
      cols = intersect(features, colnames(.)),
      names_to = ".feature",
      values_to = "value"
    )

  data$.sample = factor(data$.sample, levels = names(sce_list))
  data$.feature = factor(data$.feature, levels = features)

  if (is.numeric(data$value)) {
    min_max = data %>%
      group_by(.dots = c(facet_rows, facet_columns)) %>%
      summarize(
        min = quantile(value, lower_quantile, na.rm = TRUE),
        max = quantile(value, upper_quantile, na.rm = TRUE)
      ) %>%
      mutate(value = paste0(round(min, 2), "-", round(max, 2)))

    if (!is.null(min_value)) {
      data$value[data$value < min_value] = min_value
    }

    data = data %>%
      group_by(.dots = c(".feature")) %>%
      mutate(value = (value - quantile(value, lower_quantile, na.rm = TRUE)) / (
        quantile(value, upper_quantile, na.rm = TRUE) - quantile(value, lower_quantile, na.rm = TRUE)
      ))

    data$value[data$value > 1] =  1
    data$value[data$value < 0] =  0

  } else {
    data$value = factor(as.character(data$value),
                        levels = c("NA", mixedsort(as.character(
                          unique(data$value[data$value != "NA"])
                        ), na.last = FALSE)))
  }

  data = arrange(data,!is.na(value), value)

  plot = ggplot(data, aes_string(
    x = paste0((type), "_1"),
    y = paste0((type), "_2"),
    color = "value"
  )) +
    geom_point(alpha = alpha,
               size = point_size) +
    theme_ggexp() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  if (is.numeric(data$value)) {
    plot = plot + geom_text(
      data = min_max,
      aes(x = Inf, y = Inf, label = value),
      hjust = 1,
      vjust = 1.2,
      size = 3,
      inherit.aes = FALSE,
      show.legend = FALSE
    )
  }

  if (!is.null(label)) {
    if (label %in% features) {
      color = "value"
      label = "value"
    } else {
      color = NULL
    }

    annotations = data %>%
      group_by(.dots = c(facet_rows, facet_columns, label)) %>%
      summarize(x = median(!!as.name(paste0((
        type
      ), "_1"))),
      y = median(!!as.name(paste0((
        type
      ), "_2"))))

    plot = plot + geom_label_repel(
      data = annotations,
      aes_string(
        x = "x",
        y = "y",
        label = label,
        color = color
      ),
      label.padding = unit(0.1, "lines"),
      alpha = 1,
      fill = "white",
      size = text_size
    )
  }

  plot = plot_facets(plot,
                     facet_rows,
                     facet_columns,
                     facet_type,
                     ...)

  if (!is.numeric(data$value)) {
    plot = plot +
      guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5)))
  } else {
    plot = plot +
      scale_color_gradient(
        low = "#E8E8E8",
        high = "firebrick",
        breaks = c(0, 1),
        labels = c("min", "max"),
        limits = c(0, 1)
      )
  }

  return(plot + theme(legend.title = element_blank()))
}

#' Get reduced dimensions of object
#'
#' @param sce SingleCellExperiment object
#' @param type Name of reduction type in reducedDims
#'
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom SummarizedExperiment colData colData<-
#'
#' @return
#' @keywords internal
#'
#' @examples
#' NULL
.get_reduced_dims = function(sce, type) {
  reduced_dims = data.frame(dim1 = reducedDims(sce)@listData[[type]][, 1],
                            dim2 = reducedDims(sce)@listData[[type]][, 2])

  colnames(reduced_dims) = paste0(type, "_", c(1, 2))

  reduced_dims$barcode = colData(sce)$Barcode

  return(reduced_dims)
}
