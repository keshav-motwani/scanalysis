#' Get assay data from either the main experiment or altExps
#'
#' @param sce SingleCellExperiment object
#' @param assay assay to get data from (counts, logcounts, etc.)
#' @param alt_exp alternate experiment to get assay data from
#'
#' @importFrom SingleCellExperiment altExp
#' @importFrom SummarizedExperiment assay
#'
#' @return
#' @export
#'
#' @examples
#' NULL
get_assay_data = function(sce, assay, alt_exp = NULL) {
  if (is.null(alt_exp) ||
      alt_exp == metadata(sce)$default_assay ||
      (is.null(metadata(sce)$default_assay) && alt_exp == "RNA")) {
    data = sce
  } else {
    data = altExp(sce, alt_exp)
  }
  return(assay(data, assay))
}

#' Get rowData from either the main experiment or altExps
#'
#' @param sce SingleCellExperiment object
#' @param alt_exp alternate experiment to get rowData from
#'
#' @importFrom SingleCellExperiment altExp
#' @importFrom SummarizedExperiment assay
#'
#' @return
#' @export
#'
#' @examples
#' NULL
get_row_data = function(sce, assay, alt_exp = NULL) {
  if (is.null(alt_exp) ||
      alt_exp == metadata(sce)$default_assay ||
      (is.null(metadata(sce)$default_assay) && alt_exp == "RNA")) {
    data = sce
  } else {
    data = altExp(sce, alt_exp)
  }
  return(rowData(data))
}