plot_feature_heatmap = function(sce_list,
                                features,
                                cell_annotations,
                                feature_annotations,
                                split_cells,
                                split_features,
                                assay = "logcounts",
                                alt_exp = NULL,
                                ...) {

  data = do.call(
    cbind,
    map(sce_list,
        ~ as.matrix(get_assay_data(.x, assay, alt_exp)[features, ]))
  )

  if (is.character(feature_annotations)) {
    feature_annotations = .get_feature_annotations(sce_list, feature_annotations, alt_exp)
  }

  if (is.character(split_features)) {
    split_features = .get_feature_annotations(sce_list, split_features, alt_exp)
  }

  cell_annotations = .get_cell_annotations(sce_list, cell_annotations)

  split_cells = .get_cell_annotations(sce_list, split_cells)

  max = apply(data, 1, function(x) 1/quantile(x, 0.99))

  normalized_data = diag(max) %*% data

  heatmap = plot_heatmap(normalized_data,
                         row_annotations = feature_annotations,
                         column_annotations = cell_annotations,
                         split_rows = split_features,
                         split_columns = split_cells,
                         ...)

  return(heatmap)
}

#' Get cell annotations for feature_heatmap
#'
#' @param sce_list
#' @param cell_features
#'
#' @importFrom purrr map
#' @importFrom SingleCellExperiment colData
#'
#' @return
#' @export
#'
#' @examples
.get_cell_annotations = function(sce_list, cell_features) {
  if (!is.null(cell_features)) {
    do.call(
      cbind,
      map(sce_list,
          ~ as.data.frame(colData(.x)[, cell_features, drop = FALSE]))
    )
  }
}

.get_feature_annotations = function(sce_list, features, alt_exp) {
  as.data.frame(get_row_data(sce_list[[1]], alt_exp)[, features, drop = FALSE])
}