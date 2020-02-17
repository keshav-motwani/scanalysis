#' Plot volcano plot with annotations
#'
#' @param data differential expression results
#' @param fold_change column containing fold change data (x-axis)
#' @param p_value column containing p-value data (y-axis is -log10(p_value))
#' @param color column for color
#' @param label column for text annotation
#' @param annotations features to be annotated
#' @param annotations_if_threshold features to be annotated iff they meet the thresholds
#' @param n_annotate_top number of top DE genes to annotate
#' @param p_value_threshold maximum p-value
#' @param fold_change_threshold minimum absolute fold change threshold
#' @param p_value_rank_threshold maximum p-value rank threshold
#' @param alpha alpha of points
#' @param facet_rows columns for faceting by row
#' @param facet_columns columns for faceting by column
#' @param facet_type either "wrap" or "grid", corresponding to facet_wrap and facet_grid respectively
#' @param ... other params passed into either facet_wrap or facet_grid, depending on facet_type parameter
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
                        ...) {
  data[, paste0("log10_", p_value)] = -log10(data[, p_value])

  annotation_data = .prepare_volcano_annotations(
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
                     ...) +
    theme_ggexp()

  plot = plot + labs(x = "log(Fold Change)", y = "-log10(p-value)")

  return(plot)
}

#' Prepare annotations for volcano plot
#'
#' @param data differential expression results
#' @param fold_change column containing fold change data (x-axis)
#' @param p_value column containing p-value data (y-axis is -log10(p_value))
#' @param label column for text annotation
#' @param annotations features to be annotated
#' @param annotations_if_threshold features to be annotated iff they meet the thresholds
#' @param n_annotate_top number of top DE genes to annotate
#' @param p_value_threshold maximum p-value
#' @param fold_change_threshold minimum absolute fold change threshold
#' @param p_value_rank_threshold maximum p-value rank threshold
#' @param facet_rows columns for faceting by row
#' @param facet_columns columns for faceting by column
#'
#' @importFrom dplyr filter group_by mutate top_n filter bind_rows
#'
#' @return
#' @keywords internal
#'
#' @examples
#' NULL
.prepare_volcano_annotations = function(data,
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
  all_annotation_data = data %>%
    filter(abs(!!as.name(fold_change)) > fold_change_threshold,!!as.name(p_value) < p_value_threshold)

  requested_annotation_data = all_annotation_data %>%
    group_by(.dots = c(facet_rows, facet_columns)) %>%
    mutate(rank = rank(!!as.name(p_value))) %>%
    top_n(p_value_rank_threshold, -rank) %>%
    filter(!!as.name(label) %in% annotations_if_threshold) %>%
    bind_rows(all_annotation_data %>% filter(!!as.name(label) %in% annotations))

  top_annotation_data = .prepare_volcano_top_annotations(
    all_annotation_data,
    n_annotate_top,
    fold_change,
    p_value,
    facet_rows,
    facet_columns
  )

  annotation_data = bind_rows(top_annotation_data, requested_annotation_data)
  annotation_data$rank = NULL

  annotation_data = annotation_data[!duplicated(annotation_data),]

  return(annotation_data)
}

#' Prepare the top N annotations for volcano plot
#'
#' @param data differential expression results
#' @param n_annotate_top number of top DE genes to annotate
#' @param fold_change column containing fold change data (x-axis)
#' @param p_value column containing p-value data (y-axis is -log10(p_value))
#' @param facet_rows columns for faceting by row
#' @param facet_columns columns for faceting by column
#'
#' @importFrom dplyr arrange mutate group_by top_n
#'
#' @return
#' @keywords internal
#'
#' @examples
#' NULL
.prepare_volcano_top_annotations = function(data,
                                            n_annotate_top,
                                            fold_change,
                                            p_value,
                                            facet_rows,
                                            facet_columns) {
  arrange(data,!!as.name(p_value),-abs(!!as.name(fold_change))) %>%
    mutate(rank = seq_len(nrow(.))) %>%
    group_by(.dots = c(facet_rows, facet_columns),
             !!as.name(fold_change) > 0) %>%
    top_n(round(n_annotate_top / 2),-rank)
}
