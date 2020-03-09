#' Plot heatmap of features across cells with annotations
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param features features to use
#' @param cell_annotations columns in colData to annotate on heatmap
#' @param feature_annotations data frame of annotations related to each feature in `features`
#' @param split_cells columns in colData to split columns (cells) by
#' @param split_features data frame of groups to split rows (features) by
#' @param assay assay in SingleCellExperiment to pull data from
#' @param alt_exp altExp in SingleCellExperiment to pull data from
#' @param subsample_cells whether to subsample cells to smallest group (if split_cells specified)
#' @param show_cell_dend whether to show dendrogram for cells
#' @param show_feature_dend whether to show dendrogram for features
#' @param ... params passed to `ggexp::plot_heatmap`
#'
#' @importFrom ggexp plot_heatmap
#'
#' @return
#' @export
#'
#' @examples
#' NULL
plot_feature_heatmap = function(sce_list,
                                features,
                                cell_annotations,
                                feature_annotations,
                                split_cells,
                                split_features,
                                assay = "logcounts",
                                alt_exp = NULL,
                                subsample_cells = FALSE,
                                show_cell_dend = FALSE,
                                show_feature_dend = TRUE,
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

  max = apply(data, 1, function(x) 1/max(x))

  normalized_data = diag(max) %*% data

  rownames(normalized_data) = features

  normalized_data[normalized_data > 1] = 1
  normalized_data[normalized_data < 0] = 0

  if (subsample_cells && !is.null(split_cells)) {
    groups = apply(split_cells, 1, paste)
    min = min(table(groups))
    barcodes = data.frame(group = groups, barcode = colnames(normalized_data)) %>%
      dplyr::group_by(group) %>%
      dplyr::sample_n(min) %>%
      pull(barcode) %>%
      as.character()
    normalized_data = normalized_data[, barcodes, drop = FALSE]
    if (!is.null(cell_annotations)) {
      cell_annotations = cell_annotations[barcodes, , drop = FALSE]
    }
    if (!is.null(split_cells)) {
      split_cells = split_cells[barcodes, , drop = FALSE]
    }
  }

  heatmap = plot_heatmap(normalized_data,
                         row_annotations = feature_annotations,
                         column_annotations = cell_annotations,
                         split_rows = split_features,
                         split_columns = split_cells,
                         show_row_dend = show_feature_dend,
                         show_column_dend = show_cell_dend,
                         ...)

  return(heatmap)
}

#' Get cell annotations for feature_heatmap
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param cell_features columns from colData
#'
#' @importFrom purrr map
#' @importFrom SingleCellExperiment colData
#'
#' @return
#'
#' @examples
#' NULL
.get_cell_annotations = function(sce_list, cell_features) {
  if (!is.null(cell_features)) {
    do.call(
      cbind,
      map(sce_list,
          ~ as.data.frame(colData(.x)[, cell_features, drop = FALSE]))
    )
  }
}

#' Get feature annotations for feature_heatmap
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param features columns from rowData
#' @param alt_exp altExp to plot
#'
#' @importFrom SingleCellExperiment colData
#'
#' @return
#'
#' @examples
#' NULL
.get_feature_annotations = function(sce_list, features, alt_exp) {
  as.data.frame(get_row_data(sce_list[[1]], alt_exp)[, features, drop = FALSE])
}