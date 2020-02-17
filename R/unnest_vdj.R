#' Get VDJ data with cell-level metadata
#'
#' @param sce SingleCellExperiment object with a column named vdj containing contig_annotations in colData
#' @param include_all whether to include all metadata, including metadata about each cell (from all columns in colData)
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom dplyr as_tibble
#'
#' @return DataFrame with row for each vdj sequence containing cell-level metadata
#' @export
#'
#' @examples
#' NULL
unnest_vdj = function(sce, include_all = TRUE) {
  if ("vdj" %in% colnames(colData(sce))) {
    if (class(sce$vdj) != "CompressedSplitDFrameList") {
      warning("having to convert vdj into CompressedSplitDFrameList, this may take time")
      sce$vdj = as(sce$vdj, "CompressedSplitDFrameList")
    }
    vdj = unlist(sce$vdj)
    lengths = lengths(sce$vdj)
    indices = rep(1:ncol(sce), lengths)
    if (include_all) {
      result = cbind(colData(sce)[indices, , drop = FALSE], vdj)
    } else {
      result = vdj
      result$Barcode = colData(sce)[indices, , drop = TRUE]$Barcode
    }
    result$vdj = NULL
    return(as.data.frame(result))
  }
}