#' Compute matrix of evenness profiles per group based on clonotype frequency distributions
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param group_by column names from colData to use to construct groups for clonotype counting
#' @param min_alpha minimum Renyi entropy alpha
#' @param max_alpha maximum Renyi entropy alpha
#' @param alpha_steps number of steps in between min_alpha and max_alpha
#'
#' @importFrom purrr map_lgl map_dbl map2_dfr map
#' @importFrom dplyr as_tibble mutate group_by group_keys group_split
#' @importFrom tidyr unite
#'
#' @return
#' @export
#'
#' @examples
#' NULL
compute_evenness_profile_matrix = function(sce_list,
                                           group_by = c(),
                                           min_alpha = 0,
                                           max_alpha = 10,
                                           alpha_steps = 51) {
  if (is.null(names(sce_list)))
    names(sce_list) = paste0("sample_", 1:length(sce_list))

  stopifnot(all(map_lgl(
    sce_list, ~ "clonotype" %in% colnames(colData(.x))
  )))

  data = map2_dfr(
    sce_list,
    names(sce_list),
    ~ as_tibble(colData(.x)) %>%
      mutate(.sample = .y)
  )

  data = data %>%
    group_by(.dots = group_by)

  if (ncol(group_keys(data)) > 0) {
    groups = unite(group_keys(data), col = group,!!group_by)$group
  } else {
    groups = ""
  }

  count_vectors = data %>%
    group_split() %>%
    map( ~ dplyr::pull(.x, clonotype)) %>%
    map(table)

  matrix = t(compute_evenness_profile(count_vectors, min_alpha, max_alpha, alpha_steps))
  rownames(matrix) = groups

  attributes(matrix)$group_annotations = as.data.frame(group_keys(data))
  attributes(matrix)$group_annotations$total_reads = map_dbl(count_vectors, sum)
  attributes(matrix)$group_annotations$total_reads = map_dbl(count_vectors, sum)

  return(matrix)
}

#' Compute tidied data frame of evenness profiles per group based on clonotype frequency distributions
#'
#' @param sce_list Named list of SingleCellExperiment objects
#' @param group_by Character vector of column names from colData to use to construct groups for clonotype counting
#' @param min_alpha Minimum Renyi entropy alpha
#' @param max_alpha Maximum Renyi entropy alpha
#' @param alpha_steps Number of steps in between min_alpha and max_alpha
#'
#' @return
#' @export
#'
#' @examples
#' NULL
compute_evenness_profile_long = function(sce_list,
                                         group_by = c(),
                                         min_alpha = 0,
                                         max_alpha = 10,
                                         alpha_steps = 51) {
  evenness_profile_matrix = compute_evenness_profile_matrix(sce_list, group_by, min_alpha, max_alpha, alpha_steps)
  group_annotations = attributes(evenness_profile_matrix)$group_annotations

  evenness_profile = dplyr::as_tibble(evenness_profile_matrix) %>%
    dplyr::mutate(group = rownames(evenness_profile_matrix)) %>%
    dplyr::bind_cols(group_annotations) %>%
    tidyr::pivot_longer(
      cols = -c("group", colnames(group_annotations)),
      names_to = "alpha",
      values_to = "evenness"
    ) %>%
    dplyr::mutate(alpha = as.numeric(alpha))

  return(evenness_profile)
}

shannon_entropy = function(x) {
  p = x[!(x == 0)] / sum(x[!(x == 0)])
  - sum(p * log(p))
}

min_entropy = function(pop) {
  f = pop / sum(pop)
  - log(max(f))
}

renyi_entropy = function(pop, alpha = 1) {
  f = pop / sum(pop)
  if (abs(alpha - 1) < 10 ^ (-10)) {
    return(shannon_entropy(pop))
  } else if (alpha == Inf) {
    return(min_entropy(pop))
  } else if (abs(alpha - 1) > 10 ^ (-10)) {
    return((1 / (1 - alpha)) * log(sum(f ^ alpha)))
  }
}

compute_evenness_profile = function(count_vectors,
                                    min_alpha,
                                    max_alpha,
                                    alpha_steps) {
  alphas = seq(min_alpha, max_alpha, length.out = alpha_steps)
  matrix = sapply(1:length(count_vectors),
                  function(x)
                    sapply(1:length(alphas),
                           function(i)
                             exp(
                               renyi_entropy(count_vectors[[x]] / sum(count_vectors[[x]]), alpha = alphas[i])
                             )) / length(count_vectors[[x]]))
  rownames(matrix) = alphas
  return(matrix)
}
