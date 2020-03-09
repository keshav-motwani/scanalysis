#' Scatterplot of two features of interest from colData with annotated thresholds and counts based on filters
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param x_filters list of filters for each SCE object in sce_list from scater::isOutlier, or a vector with attribute named thresholds that is a vector with min and max allowed values
#' @param x numeric column from colData that is in all objects in sce_list
#' @param facet_rows columns to facet on
#' @param facet_columns columns to facet on
#' @param ... other params passed into plot_gex_univariate_qc
#'
#' @import ggplot2
#' @importFrom purrr map map2 imap_dfr map_lgl
#' @importFrom dplyr summarise select group_by arrange summarize desc ungroup distinct
#' @importFrom SingleCellExperiment colData
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
                                      facet_rows = NULL,
                                      facet_columns = NULL,
                                      ...) {

  chain_count_annotated = all(map_lgl(sce_list, ~ "count_TRA_TRB_IGL_IGK_IGH" %in% colnames(colData(.x))))

  if (!chain_count_annotated) {
    stop("Chain count annotations not available, please run annotate_chain_count on each object in sce_list first")
  }

  plot = plot_gex_univariate_qc(sce_list,
                                x_filters = x_filters,
                                x = x,
                                y = "count_TRA_TRB_IGL_IGK_IGH",
                                facet_rows = facet_rows,
                                facet_columns = facet_columns,
                                ...)

  if (!is.null(vdj_filters)) {

    order = imap_dfr(vdj_filters, ~ attributes(.x)$allowed_values %>%
                       mutate(.sample = .y)) %>%
      group_by(.sample) %>%
      arrange(desc(allowed), count_TRA_TRB_IGL_IGK_IGH)

    yintercepts = summarize(order, yintercept = sum(allowed) + 0.9)

    plot$data$count_TRA_TRB_IGL_IGK_IGH = factor(plot$data$count_TRA_TRB_IGL_IGK_IGH,
                                                 levels = distinct(ungroup(order), count_TRA_TRB_IGL_IGK_IGH, allowed)$count_TRA_TRB_IGL_IGK_IGH)

    plot = plot + geom_hline(data = yintercepts,
                             aes(yintercept = yintercept),
                             linetype = "dashed")

  }

  return(plot)
}