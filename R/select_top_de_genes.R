#' Select the top DE genes, ranked on either fold change or -log10(p-value)
#'
#' @param data differential expression results
#' @param n number of top DE genes to select
#' @param fold_change column containing fold change data (x-axis)
#' @param p_value column containing p-value data (y-axis is -log10(p_value))
#' @param groups column containing groups to get top DE genes per group (for example, cluster column)
#' @param rank_by whether to rank by fold_change or p_value first
#' @param only_pos only select features with positive fold change
#'
#' @importFrom dplyr arrange mutate group_by top_n
#'
#' @return
#' @export
#'
#' @examples
#' NULL
select_top_de_genes = function(data,
                               n = 10,
                               fold_change = "fold_change",
                               p_value = "p_value",
                               groups = c(),
                               rank_by = "fold_change",
                               only_pos = FALSE) {
  if (!only_pos) {
    data$abs_fold_change = abs(data[, fold_change, drop = TRUE])
  } else {
    data$abs_fold_change = data[, fold_change, drop = TRUE]
  }

  data$fold_change = data[, fold_change, drop = TRUE]
  data$p_value = data[, p_value, drop = TRUE]

  if (rank_by == "fold_change") {
    result = arrange(data,-abs_fold_change, p_value)
  } else {
    result = arrange(data, p_value,-abs_fold_change)
  }

  result = result %>%
    mutate(rank = seq_len(nrow(.)))

  if (!only_pos) {
    result = result %>%
      group_by(.dots = groups,
               fold_change > 0) %>%
      top_n(round(n / 2), -rank)
  } else {
    result = result %>%
      group_by(.dots = groups) %>%
      top_n(n, -rank)
  }

  return(result)
}
