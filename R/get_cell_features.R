#' Get feature from assay data, colData, or reducedDims at once from main experiment or alternate experiments
#'
#' @param sce SingleCellExperiment object
#' @param features Names of features desired. It will first check row names of the assay data, then column names of the colData, then column names of each of the reducedDims slots. Only those found will be returned
#' @param assay Assay to use (counts, logcounts, etc.)
#' @param alt_exp Name of the altExp to use (if any)
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr bind_cols
#'
#' @return data.frame with rows as cells and columns as features
#'
#' @examples
#' NULL
get_cell_features = function(sce, features, assay, alt_exp = NULL) {
  assay_data = .get_cell_features_assay(sce, features, assay, alt_exp)
  coldata_data = .get_cell_features_coldata(sce, features)
  reduced_dimensions_data = .get_cell_features_reduced_dimensions(sce, features)
  result = bind_cols(assay_data, coldata_data, reduced_dimensions_data)
  result$barcode = colData(sce)$Barcode
  return(result)
}

#' Get feature from assay (from alternate experiment)
#'
#' @param sce SingleCellExperiment object
#' @param features Row names from assay object of relevant alternate experiment
#' @param assay Assay to use (counts, logcounts, etc.)
#' @param alt_exp Name of the altExp to use (if any)
#'
#' @return data.frame with rows as cells and columns as features
#' @keywords internal
#'
#' @examples
#' NULL
.get_cell_features_assay = function(sce, features, assay, alt_exp) {
  rownames_intersection = intersect(features, rownames(get_assay_data(sce, assay, alt_exp)))
  matrix = t(as.matrix(get_assay_data(sce, assay, alt_exp)[rownames_intersection, , drop = FALSE]))
  feature_values = data.frame(matrix)
  colnames(feature_values) = rownames_intersection
  return(feature_values)
}

#' Get feature from column metadata
#'
#' @param sce SingleCellExperiment object
#' @param features Column names from column metadata
#'
#' @return data.frame with rows as cells and columns as features
#' @keywords internal
#'
#' @examples
#' NULL
.get_cell_features_coldata = function(sce, features) {
  colnames_intersection = intersect(features, colnames(SummarizedExperiment::colData(sce)))
  feature_values = as.data.frame(SummarizedExperiment::colData(sce)[, colnames_intersection, drop = FALSE])
  return(feature_values)
}

#' Get feature from reduced dimensional representation
#'
#' @param sce SingleCellExperiment object
#' @param features Names of dimensions to retreive in format of [red_dim_name]_[dimension_number] (e.g. PCA_1)
#'
#' @importFrom purrr map_dfc
#' @importFrom SingleCellExperiment reducedDims
#'
#' @return data.frame with rows as cells and columns as features
#' @keywords internal
#'
#' @examples
#' NULL
.get_cell_features_reduced_dimensions = function(sce, features) {
  if (length(reducedDims(sce)) > 0) {
    reduced_dimensions_names = names(reducedDims(sce))
    reduced_dimensions = map_dfc(reduced_dimensions_names, ~ {
      data = as.data.frame(reducedDim(sce, .x))
      colnames(data) = paste0(.x, "_", 1:ncol(data))
      data
    })
    colnames_intersection = intersect(features, colnames(reduced_dimensions))
    feature_values = as.data.frame(reduced_dimensions[, colnames_intersection, drop = FALSE])
  } else {
    feature_values = NULL
  }
  return(feature_values)
}
