#' Scatterplot of two features of interest from colData with annotated thresholds and counts based on filters
#'
#' @param sce_list List of SingleCellExperiment object
#' @param x Numeric column from colData that is in all objects in sce_list
#' @param y Numeric column from colData that is in all objects in sce_list
#' @param color Column from colData that is in all objects in sce_list
#' @param shape Column from colData that is in all objects in sce_list
#' @param x_filters List of filters for each SCE object in sce_list from scater::isOutlier, or a vector with attribute named thresholds that is a vector with min and max allowed values
#' @param y_filters List of filters for each SCE object in sce_list from scater::isOutlier, or a vector with attribute named thresholds that is a vector with min and max allowed values
#' @param x_log Boolean to use log x-axis
#' @param y_log Boolean to use log y-axis
#' @param text_size Font size for annotations
#' @param facet_rows Columns to facet on
#' @param facet_columns Columns to facet on
#' @param facet_type Either "wrap" or "grid", same as ggplot
#' @param ... params passed into either facet_wrap or facet_grid, depending on facet_type parameter
#'
#' @import ggplot2
#' @importFrom ggexp theme_ggexp plot_facets
#' @importFrom purrr imap pmap_dfr
#' @importFrom dplyr bind_rows
#' @importFrom scales trans_breaks trans_format math_format
#'
#' @return
#' @export
#'
#' @examples
#' NULL
plot_gex_bivariate_qc = function(sce_list,
                                 x,
                                 y,
                                 color = NULL,
                                 shape = NULL,
                                 x_filters = NULL,
                                 y_filters = NULL,
                                 x_log = TRUE,
                                 y_log = TRUE,
                                 text_size = 3,
                                 facet_rows = NULL,
                                 facet_columns = NULL,
                                 facet_type = "grid",
                                 ...) {

  if (is.null(names(sce_list))) {
    names(sce_list) = paste0("sample_", 1:length(sce_list))
  }

  if (is.null(x_filters)) {
    x_filters = rep(list(NULL), length(sce_list))
  } else {
    stopifnot(names(sce_list) == names(x_filters))
  }

  if (is.null(y_filters)) {
    y_filters = rep(list(NULL), length(sce_list))
  } else {
    stopifnot(names(sce_list) == names(y_filters))
  }

  data = imap(sce_list, .prepare_gex_data)

  counts = pmap_dfr(
    list(data, x_filters, y_filters, names(data)),
    ~ .prepare_gex_bivariate_counts(
      ..1,
      x,
      y,
      ..2,
      ..3,
      x_log,
      y_log,
      ..4,
      c(facet_rows, facet_columns)
    )
  )

  data = bind_rows(data)

  plot = data %>%
    ggplot(aes_string(
      x = x,
      y = y,
      color = color,
      shape = shape
    )) +
    geom_point(alpha = 0.5) +
    geom_density2d(color = "blue", alpha = 0.5) +
    theme_ggexp() +
    scale_color_viridis_c()

  xlim = c(min(data[, x, drop = TRUE]), max(data[, x, drop = TRUE]))
  ylim = c(min(data[, y, drop = TRUE]), max(data[, y, drop = TRUE]))

  if (x_log) {
    plot = plot +
      scale_x_log10(
        breaks = trans_breaks("log10", function(x)
          10 ^ x),
        labels = trans_format("log10", math_format(10 ^ .x)),
        limits = xlim
      )
  } else {
    plot = plot + xlim(xlim)
  }

  if (y_log) {
    plot = plot + scale_y_log10(
      breaks = trans_breaks("log10", function(x)
        10 ^ x),
      labels = trans_format("log10", math_format(10 ^ .x)),
      limits = ylim
    )
  }

  plot = plot +
    geom_label(
      data = counts %>% filter(count != 0),
      aes(
        label = count,
        x = x,
        y = y,
        vjust  = vjust,
        hjust = hjust
      ),
      color = "black",
      size = text_size,
      label.padding = unit(0.1, "lines"),
      alpha = 0.5
    ) +
    geom_vline(data = counts %>% filter(x1_anno == 1),
               aes(xintercept = x1),
               linetype = "dashed") +
    geom_vline(data = counts %>% filter(x2_anno == 1),
               aes(xintercept = x2),
               linetype = "dashed") +
    geom_hline(data = counts %>% filter(y1_anno == 1),
               aes(yintercept = y1),
               linetype = "dashed") +
    geom_hline(data = counts %>% filter(y2_anno == 1),
               aes(yintercept = y2),
               linetype = "dashed")

  plot = plot_facets(plot,
                     facet_rows,
                     facet_columns,
                     facet_type,
                     ...)

  return(plot)
}

#' Scatterplot of two features of interest from colData with annotated thresholds and counts based on filters
#'
#' @param sce_list List of SingleCellExperiment object
#' @param x_filters List of filters for each SCE object in sce_list from scater::isOutlier, or a vector with attribute named thresholds that is a vector with min and max allowed values
#' @param x Numeric column from colData that is in all objects in sce_list
#' @param y Discrete column from colData that is in all objects in sce_list to split histograms by
#' @param color Column from colData that is in all objects in sce_list
#' @param shape Column from colData that is in all objects in sce_list
#' @param x_log Boolean to use log x-axis
#' @param text_size Font size for annotations
#' @param facet_rows Columns to facet on
#' @param facet_columns Columns to facet on
#' @param facet_type Either "wrap" or "grid", same as ggplot
#' @param ... params passed into either facet_wrap or facet_grid, depending on facet_type parameter
#'
#' @import ggplot2
#' @importFrom ggexp theme_ggexp plot_facets
#' @importFrom purrr imap pmap_dfr
#' @importFrom dplyr bind_rows
#' @importFrom scales trans_breaks trans_format math_format
#' @importFrom ggridges geom_density_ridges2
#'
#' @return
#' @export
#'
#' @examples
#' NULL
plot_gex_univariate_qc = function(sce_list,
                                  x_filters = NULL,
                                  x,
                                  y = NULL,
                                  color = NULL,
                                  shape = NULL,
                                  x_log = TRUE,
                                  text_size = 3,
                                  facet_rows = NULL,
                                  facet_columns = NULL,
                                  facet_type = "wrap",
                                  ...) {
  data = imap(sce_list, .prepare_gex_data)

  counts = pmap_dfr(
    list(data, x_filters, names(data)),
    ~ .prepare_gex_univariate_counts(..1, x, ..2, x_log, ..3, c(facet_rows, facet_columns, y))
  )

  data = bind_rows(data)

  if (is.null(y)) {
    y = ".null"
    data$.null = "all"
    counts$.null = "all"
  }

  plot = ggplot(data, aes_string(x = x, y = y)) +
    geom_density_ridges2(
      aes_string(point_color = color, point_shape = shape),
      alpha = .2,
      point_alpha = 0.5,
      jittered_points = TRUE
    ) +
    theme_ggexp()

  if (x_log) {
    plot = plot +
      scale_x_log10(
        breaks = trans_breaks("log10", function(x)
          10 ^ x),
        labels = trans_format("log10", math_format(10 ^ .x)),
        limits = c(min(data[, x, drop = TRUE]), max(data[, x, drop = TRUE]))
      )
  }

  if (!is.null(x_filters)) {
    plot = plot +
      geom_label(
        data = counts %>% filter(count != 0),
        aes(
          label = count,
          x = x,
          hjust = hjust
        ),
        color = "black",
        size = text_size,
        label.padding = unit(0.1, "lines"),
        alpha = 0.5,
        vjust = -0.5
      ) +
      geom_vline(
        data = counts %>% filter(x1_anno == 1),
        aes(xintercept = x1),
        linetype = "dashed"
      ) +
      geom_vline(
        data = counts %>% filter(x2_anno == 1),
        aes(xintercept = x2),
        linetype = "dashed"
      )
  }

  plot = plot_facets(plot,
                     facet_rows,
                     facet_columns,
                     facet_type,
                     ...)

  if (y == ".null") {
    plot = plot + labs(y = NULL)
  }

  if (is.numeric(data[, color, drop = TRUE])) {
    plot = plot +
      scale_color_viridis_c(aesthetics = c("point_colour"),
                            guide = guide_colorbar(available_aes = c("point_colour")))
  } else {
    plot = plot +
      scale_color_viridis_d(aesthetics = c("point_colour"))
  }

  plot = plot + scale_y_discrete(expand = c(0, 0.005))

  return(plot)
}

#' Prepare data frame to plot for gene expression related QC plots
#'
#' @param sce SingleCellExperiment object
#' @param sample_name Name of sample
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom dplyr mutate
#'
#' @return
#' @keywords internal
#'
#' @examples
#' NULL
.prepare_gex_data = function(sce, sample_name) {
  colData(sce) %>%
    as.data.frame() %>%
    mutate(.sample = sample_name)
}

#' Prepare counts of cells that meet filter criteria to annotate based on a single filter
#'
#' @param data Result from .prepare_gex_data
#' @param x Feature from colData to plot distribution of
#' @param x_filter Filter for a sample from scater::isOutlier, or a vector with attribute named thresholds that is a vector with min and max allowed values
#' @param x_log Boolean to plot x on log-scale or not
#' @param sample_name Sample name
#' @param facets Features from colData to facet on
#'
#' @importFrom dplyr mutate group_by summarize left_join
#' @importFrom tidyr unnest
#'
#' @return
#' @keywords internal
#'
#' @examples
#' NULL
.prepare_gex_univariate_counts = function(data,
                                          x,
                                          x_filter,
                                          x_log,
                                          sample_name,
                                          facets) {
  x_scale_limits = .prepare_gex_scale_limits(data, x, x_filter)

  counts = data.frame(
    filter = paste0("filter_", 1:3),
    x1 = x_scale_limits[1:3],
    x1_anno = ifelse(is.null(x_filter), 0, list(c(0, 1, 1)))[[1]],
    x2 = x_scale_limits[2:4],
    x2_anno = ifelse(is.null(x_filter), 0, list(c(1, 1, 0)))[[1]],
    hjust = 0.5,
    vjust = 0.5
  )

  count_fn = function(x) {
    data.frame(
      filter = paste0("filter_", 1:3),
      count = apply(as.matrix(counts[, c("x1", "x2")]), 1,
                    function(row)
                      sum(x >= row["x1"] &
                            x < row["x2"]))
    )
  }

  counts = data %>%
    group_by(.dots = facets) %>%
    summarize(counts = list(count_fn(!!as.name(x)))) %>%
    unnest() %>%
    left_join(counts)

  if (x_log) {
    counts = counts %>%
      mutate(x = sqrt(ifelse(x1 == 0, 1, x1) * (x2)))
  } else {
    counts = counts %>%
      mutate(x = 0.5 * (x1 + x2))
  }

  counts$.sample = sample_name

  return(counts)
}

#' Prepare counts of cells that meet filter criteria to annotate based on two filters
#'
#' @param data Result from .prepare_gex_data
#' @param x Feature from colData to plot distribution of
#' @param x_filter Filter for a sample from scater::isOutlier, or a vector with attribute named thresholds that is a vector with min and max allowed values
#' @param y Feature from colData to plot distribution of
#' @param y_filter Filter for a sample from scater::isOutlier, or a vector with attribute named thresholds that is a vector with min and max allowed values
#' @param x_log Boolean to plot x on log-scale or not
#' @param y_log Boolean to plot y on log-scale or not
#' @param sample_name Sample name
#' @param facets Features from colData to facet on
#'
#' @importFrom dplyr mutate group_by summarize left_join
#' @importFrom tidyr unnest
#'
#' @return
#' @keywords internal
#'
#' @examples
#' NULL
.prepare_gex_bivariate_counts = function(data,
                                         x,
                                         y,
                                         x_filter,
                                         y_filter,
                                         x_log,
                                         y_log,
                                         sample_name,
                                         facets) {
  x_scale_limits = .prepare_gex_scale_limits(data, x, x_filter)
  y_scale_limits = .prepare_gex_scale_limits(data, y, y_filter)

  counts = data.frame(
    filter = paste0("filter_", 1:9),
    x1 = x_scale_limits[1:3],
    x1_anno = ifelse(is.null(x_filter), 0, list(c(0, 1, 1)))[[1]],
    x2 = x_scale_limits[2:4],
    x2_anno = ifelse(is.null(x_filter), 0, list(c(1, 1, 0)))[[1]],
    y1 = rep(y_scale_limits[1:3], each = 3),
    y1_anno = ifelse(is.null(y_filter), 0, list(rep(c(
      0, 1, 1
    ), each = 3)))[[1]],
    y2 = rep(y_scale_limits[2:4], each = 3),
    y2_anno = ifelse(is.null(y_filter), 0, list(rep(c(
      1, 1, 0
    ), each = 3)))[[1]],
    hjust = 0.5,
    vjust = 0.5
  )

  count_fn = function(x, y) {
    data.frame(
      filter = paste0("filter_", 1:9),
      count = apply(as.matrix(counts[, c("x1", "x2", "y1", "y2")]), 1,
                    function(row)
                      sum(x >= row["x1"] &
                            x < row["x2"] &
                            y > row["y1"] &
                            y <= row["y2"]))
    )
  }

  counts = data %>%
    group_by(.dots = facets) %>%
    summarize(counts = list(count_fn(!!as.name(x), !!as.name(y)))) %>%
    unnest() %>%
    left_join(counts)

  if (x_log) {
    counts = counts %>%
      mutate(x = sqrt(ifelse(x1 == 0, 1, x1) * (x2)))
  } else {
    counts = counts %>%
      mutate(x = 0.5 * (x1 + x2))
  }

  if (y_log) {
    counts = counts %>%
      mutate(y = sqrt(ifelse(y1 == 0, 1, y1) * (y2)))
  } else {
    counts = counts %>%
      mutate(x = 0.5 * (y1 + y2))
  }

  counts$.sample = sample_name

  return(counts)
}


#' Prepare the limits of plot scale
#'
#' @param data Result from .prepare_gex_data
#' @param var Feature from data that is being plotted
#' @param filter Filter for a sample from scater::isOutlier, or a vector with attribute named thresholds that is a vector with min and max allowed values
#'
#' @return
#' @keywords internal
#'
#' @examples
#' NULL
.prepare_gex_scale_limits = function(data, var, filter) {
  if (!is.null(filter)) {
    scale_limits = c(
      min =  min(data[, var, drop = TRUE]),
      lower = max(attributes(filter)$thresholds[1], min(data[, var, drop = TRUE])),
      higher = min(attributes(filter)$thresholds[2], max(data[, var, drop = TRUE])),
      max = max(data[, var, drop = TRUE])
    )
  } else {
    scale_limits = c(
      min =  min(data[, var, drop = TRUE]),
      lower = min(data[, var, drop = TRUE]),
      higher = max(data[, var, drop = TRUE]),
      max = max(data[, var, drop = TRUE])
    )
  }
  return(scale_limits)
}
