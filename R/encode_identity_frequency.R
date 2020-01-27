#' Encode frequencies of (combinations) of values in columns from colData into a matrix with rows as groups and columns as features
#'
#' @param sce_list List of SCE objects, all containing relevant attributes in colData as columns
#' @param attributes Column names from colData which to encode
#' @param group_by Column names from colData which to group by (compute frequencies per group)
#' @param normalize Normalization method - either "none" or "relative_frequency"
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom purrr map2_dfr
#' @importFrom dplyr mutate
#'
#' @return
#' @export
#'
#' @examples
#' NULL
encode_cell_identity_frequency_matrix = function(sce_list,
                                                 attributes,
                                                 group_by = c(),
                                                 normalize = "none") {
  if (is.null(names(sce_list)))
    names(sce_list) = paste0("sample_", 1:length(sce_list))

  data = map2_dfr(
    sce_list,
    names(sce_list),
    ~ as.data.frame(colData(.x)) %>%
      mutate(.sample = .y)
  )

  group_by = unique(c(group_by, ".sample"))

  return(.encode_identity_frequency_matrix(data, attributes, group_by, normalize))
}

#' Encode frequencies of (combinations) of values in columns from colData into a long data frame
#'
#' @param sce_list List of SCE objects, all containing relevant attributes in colData as columns
#' @param attributes Column names from colData which to encode
#' @param group_by Column names from colData which to group by (compute frequencies per group)
#' @param normalize Normalization method - either "none" or "relative_frequency"
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom purrr imap_dfr
#' @importFrom dplyr mutate
#'
#' @return
#' @export
#'
#' @examples
#' NULL
encode_cell_identity_frequency_long = function(sce_list,
                                               attributes,
                                               group_by = c(),
                                               normalize = "none") {
  if (is.null(names(sce_list)))
    names(sce_list) = paste0("sample_", 1:length(sce_list))

  data = imap_dfr(sce_list,
                  ~ as.data.frame(colData(.x)) %>%
                    mutate(.sample = .y))

  group_by = unique(c(group_by, ".sample"))

  return(convert_identity_frequency_matrix_to_long(data, attributes, group_by, normalize))
}

#' Encode frequencies of (combinations) of values in columns from colData into a matrix with rows as groups and columns as features
#'
#' @param sce_list List of SCE objects, all containing relevant attributes in colData(sce)$vdj as columns
#' @param attributes Column names from unnested colData(sce)$vdj which to encode
#' @param group_by Column names from colData which to group by (compute frequencies per group)
#' @param normalize Normalization method - either "none" or "relative_frequency"
#'
#' @importFrom dplyr mutate
#' @importFrom purrr imap_dfr
#'
#' @return
#' @export
#'
#' @examples
#' NULL
encode_vdj_identity_frequency_matrix = function(sce_list,
                                                attributes,
                                                group_by = c(),
                                                normalize = "none") {
  data = imap_dfr(sce_list,
                  ~ unnest_vdj(.x) %>%
                    mutate(.sample = .y))

  group_by = unique(c(group_by, ".sample"))

  return(.encode_identity_frequency_matrix(data, attributes, group_by, normalize))
}

#' Encode frequencies of (combinations) of values in columns from colData into long data frame
#'
#' @param sce_list List of SCE objects, all containing relevant attributes in colData(sce)$vdj as columns
#' @param attributes Column names from unnested colData(sce)$vdj which to encode
#' @param group_by Column names from colData which to group by (compute frequencies per group)
#' @param normalize Normalization method - either "none" or "relative_frequency"
#'
#' @importFrom dplyr mutate
#' @importFrom purrr imap_dfr
#'
#' @return
#' @export
#'
#' @examples
#' NULL
encode_vdj_identity_frequency_long = function(sce_list,
                                              attributes,
                                              group_by = c(),
                                              normalize = "none") {
  data = imap_dfr(sce_list,
                  ~ unnest_vdj(.x) %>%
                    mutate(.sample = .y))

  return(convert_identity_frequency_matrix_to_long(data, attributes, group_by, normalize))
}


#' Helper function to encode identity frequency matrix given data frame
#'
#' @param data data frame to encode
#' @param attributes columns to use to define features
#' @param group_by columns to use to define groups
#' @param normalize normalization method
#'
#' @importFrom dplyr group_by group_keys group_split pull bind_rows
#' @importFrom tidyr unite separate
#' @importFrom purrr map
#' @importFrom plyr ldply
#'
#' @return
#'
#' @examples
#' NULL
.encode_identity_frequency_matrix = function(data, attributes, group_by, normalize) {
  data[, attributes][is.na(data[, attributes])] = "NA"

  data = data %>%
    group_by(.dots = group_by) %>%
    unite(identity, !!attributes, remove = FALSE, sep = "///")

  row_annotations = group_keys(data)

  if (ncol(row_annotations) > 0) {
    groups = unite(row_annotations, col = group, !!group_by)$group
  } else {
    groups = ""
    row_annotations = data.frame(annotation = c("1"))
  }

  count_vectors = data %>%
    group_split() %>%
    map(~ pull(.x, identity)) %>%
    map(~ table(.x, useNA = "no"))
  names(count_vectors) = groups

  data = ldply(count_vectors, bind_rows)
  data$.id = NULL
  data[, "NA"] = NULL

  matrix = as.matrix(data)
  matrix[is.na(matrix)] = 0
  if (normalize == "relative_frequency") {
    diag = diag(1 / rowSums(matrix))
    if (nrow(matrix) == 1)
      diag = matrix(1 / rowSums(matrix))
    matrix = diag %*% matrix
  }
  rownames(matrix) = groups

  column_annotations = data.frame(feature = colnames(matrix))

  attributes(matrix)$group_annotations = row_annotations
  attributes(matrix)$feature_annotations = separate(
    column_annotations,
    feature,
    into = attributes,
    sep = "///",
    remove = FALSE
  )

  return(matrix)
}

#' Convert an identity frequency matrix
#'
#' @param data data frame to encode
#' @param attributes columns to use to define features
#' @param group_by columns to use to define groups
#' @param normalize normalization method
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr bind_cols left_join
#'
#' @return
#' @export
#'
#' @examples
#' NULL
convert_identity_frequency_matrix_to_long = function(data, attributes, group_by, normalize) {
  matrix = .encode_identity_frequency_matrix(data, attributes, group_by, normalize)
  data = bind_cols(
    data.frame(attributes(matrix)$group_annotations),
    as.data.frame(matrix),
    data.frame(group = rownames(matrix))
  ) %>%
    pivot_longer(cols = colnames(matrix), names_to = "feature") %>%
    left_join(attributes(matrix)$feature_annotations, by = "feature")
  return(data)
}
