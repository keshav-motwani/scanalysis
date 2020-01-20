#' Get assay data from either the main experiment or altExps
#'
#' @param sce SingleCellExperiment object
#' @param assay Assay to get data from (counts, logcounts, etc.)
#' @param alt_exp Alternate experiment to get assay data from
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
  if (!is.null(alt_exp)) {
    data = altExp(sce, alt_exp)
  } else {
    data = sce
  }
  return(assay(data, assay))
}