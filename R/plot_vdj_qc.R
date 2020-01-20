#' Plot number of VDJ chains versus colData variable of interest with annotated filters
#'
#' @param sce SingleCellExperiment object
#' @param y Numeric column from colData
#' @param vdj_doublet_filter Filter returned from filter_vdj_doublets
#' @param y_filter Filter result from scater::isOutlier or from scanalysis::filter_gex_outliers
#' @param y_log Boolean to use log y-axis
#'
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom ggexp theme_ggexp
#' @importFrom dplyr mutate as_tibble group_by tally
#'
#' @return
#' @export
#'
#' @examples
#' NULL
plot_vdj_qc = function(sce, y, vdj_doublet_filter, y_filter, y_log = TRUE, y_label = y) {
  data = as_tibble(colData(sce)) %>%
    mutate(chain_count = vdj_doublet_filter)
  plot = ggplot(data, aes_string(x = "count_TRA_TRB_IGL_IGH", y = y, color = "chain_count")) +
    geom_quasirandom(size = 0.5) +
    theme_ggexp() +
    scale_color_manual(values = c("gray", "black"), guide = "none") +
    geom_hline(aes(yintercept = attributes(y_filter)$threshold["lower"])) +
    geom_hline(aes(yintercept = attributes(y_filter)$threshold["higher"]))
  counts = data %>%
    group_by(.dots = "count_TRA_TRB_IGL_IGH") %>%
    tally()
  plot = plot + geom_text(
    data = counts,
    aes_string(label = "n",
               x = "count_TRA_TRB_IGL_IGH",
               y = min(data[, y])),
    hjust = 0.5,
    vjust = 1,
    size = 2,
    color = "black",
    angle = 0
  ) + labs(x = "Count: TRA_TRB_IGL_IGH", y = y_label)
  if (y_log) {
    plot = plot + scale_y_log10()
  }
  return(plot)
}
