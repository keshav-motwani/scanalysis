#' Plot volcano plot with annotations
#'
#' @param data data frame containing dataset to use for plotting
#' @param fold_change string scalar indicating column containing fold change data (x-axis)
#' @param p_value string scalar indicating column containing p-value data (y-axis is -log10(p_value))
#' @param color string scalar indicating column for color
#' @param label string scalar indicating column for text annotation
#' @param annotations string vector indicating features to be annotated
#' @param annotations_if_threshold string vector indicating features to be annotated iff they meet the thresholds
#' @param n_annotate_top integer indicating the number of top DE genes to annotate
#' @param p_value_threshold float indicating maximum p-value
#' @param fold_change_threshold float indicating minimum absolute fold change threshold
#' @param p_value_rank_threshold integer indicating maximum p-value rank threshold
#' @param alpha numeric scalar for alpha of points
#' @param facet_rows string vector indicating columns for faceting by row
#' @param facet_columns string vector indicating columns for faceting by column
#' @param facet_type string scalar that is either "wrap" or "grid", corresponding to facet_wrap and facet_grid respectively
#' @param facet_switch string scalar that is either NULL, "both", "x", or "y", same as switch argument in facet calls
#' @param facet_scales string scalar that is either "fixed", "free_x", "free_y", or "free", same as scales argument in facet calls
#' @param nrow numeric scalar indicating the number of rows in plot, only applies if facet_type == "wrap"
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggexp plot_facets theme_ggexp
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' NULL
plot_volcano = function(data,
                        fold_change,
                        p_value,
                        color = NULL,
                        label = "gene",
                        annotations = c(),
                        annotations_if_threshold = c(),
                        n_annotate_top = 10,
                        p_value_threshold = 0.05,
                        fold_change_threshold = 1.5,
                        p_value_rank_threshold = 100,
                        alpha = 0.5,
                        facet_rows = c(),
                        facet_columns = c(),
                        facet_type = "grid",
                        facet_scales = "free",
                        facet_switch = NULL,
                        nrow = 3) {
  data[, paste0("log10_", p_value)] = -log10(data[, p_value])

  annotation_data = get_annotations(
    data,
    fold_change,
    p_value,
    label,
    n_annotate_top,
    annotations,
    annotations_if_threshold,
    p_value_threshold,
    fold_change_threshold,
    p_value_rank_threshold,
    facet_rows,
    facet_columns
  )

  annotation_data$color = ifelse(annotation_data[, fold_change, drop = TRUE] > 0, "firebrick", "cornflowerblue")

  plot = data %>%
    ggplot(aes_string(
      x = fold_change,
      y = paste0("log10_", p_value),
      color = color,
      label = label
    )) +
    geom_point(alpha = alpha, shape = 1) +
    geom_text_repel(data = annotation_data,
                    color = annotation_data$color,
                    size = 3) +
    geom_point(data = annotation_data, color = annotation_data$color)

  if (!is.null(fold_change_threshold)) {
    plot = plot +
      geom_vline(xintercept = fold_change_threshold, linetype = 2) +
      geom_vline(xintercept = -1 * fold_change_threshold,
                 linetype = 2)
  }

  if (!is.null(p_value_threshold)) {
    plot = plot +
      geom_hline(yintercept = -1 * log10(p_value_threshold),
                 linetype = 2)
  }

  plot = plot_facets(plot,
                     facet_rows,
                     facet_columns,
                     facet_type,
                     facet_scales,
                     facet_switch,
                     nrow) + theme_ggexp()
  plot = plot + labs(x = "log2(Fold Change)", y = "-log10(p-value)")

  return(plot)
}

get_annotations = function(data,
                           fold_change,
                           p_value,
                           label,
                           n_annotate_top,
                           annotations,
                           annotations_if_threshold,
                           p_value_threshold,
                           fold_change_threshold,
                           p_value_rank_threshold,
                           facet_rows,
                           facet_columns) {
  annotation_data = data %>%
    dplyr::filter(abs(!!as.name(fold_change)) > fold_change_threshold) %>%
    dplyr::group_by(.dots = c(facet_rows, facet_columns)) %>%
    dplyr::mutate(rank = rank(!!as.name(p_value))) %>%
    dplyr::top_n(p_value_rank_threshold,-rank) %>%
    dplyr::filter(
      !!as.name(label) %in% annotations_if_threshold,!!as.name(p_value) < p_value_threshold,
    ) %>%
    dplyr::bind_rows(data %>% dplyr::filter(!!as.name(label) %in% annotations)) %>%
    dplyr::bind_rows(
      get_top_annotations(
        data,
        n_annotate_top,
        fold_change,
        p_value,
        facet_rows,
        facet_columns
      )
    )

  annotation_data = annotation_data[!duplicated(annotation_data), ]

  annotation_data
}

get_top_annotations = function(data,
                               n_annotate_top,
                               fold_change,
                               p_value,
                               facet_rows,
                               facet_columns) {
  dplyr::arrange(data, !!as.name(p_value),-abs(!!as.name(fold_change))) %>%
    dplyr::mutate(rank = seq_len(nrow(.))) %>%
    dplyr::group_by(.dots = c(facet_rows, facet_columns),!!as.name(fold_change) > 0) %>%
    dplyr::top_n(round(n_annotate_top / 2), -rank)
}
