#' Get filter for barcodes based on the number of chains for TRA, TRB, IGL, and IGH present
#'
#' Wrapper around scater::isOutlier, arguments type and nmads refer to the same thing. For mitochondrial reads, uses non-log transformed values for deciding cutoffs.
#'
#' @param sce SingleCellExperiment object
#' @param min_tra Range of unique TRA sequences allowed per barcode. Barcodes with less than the first number and greater than the second number will be excluded from filter.
#' @param min_trb Range of unique TRB sequences allowed per barcode. Barcodes with less than the first number and greater than the second number will be excluded from filter.
#' @param min_igl Range of unique IGL sequences allowed per barcode. Barcodes with less than the first number and greater than the second number will be excluded from filter.
#' @param min_igh Range of unique IGH sequences allowed per barcode. Barcodes with less than the first number and greater than the second number will be excluded from filter.
#' @param tcr_and_bcr_allowed Boolean specifying if reads to both TRA/TRB locus and IGL/IGH locus are allowed on a single barcode. If FALSE, barcodes with reads to both will be excluded from filter.
#' @param only_filter_if_any_vdj_hits Boolean specifying if barcodes should only be filtered if there are any TRA/TRB/IGL/IGH reads at all. For diverse populations with non B/T cells, should be set to TRUE in order to avoid meaningless filtering on non B/T cells.
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom scater isOutlier
#' @importFrom Matrix rowSums
#'
#' @return Boolean filter with TRUE for cells kept and FALSE for cells outside the allowed ranges.
#' @export
#'
#' @examples
#' NULL
filter_chain_count = function(sce,
                                  tra_range = c(0, 2),
                                  trb_range = c(0, 1),
                                  igl_range = c(0, 2),
                                  igh_range = c(0, 1),
                                  tcr_and_bcr_allowed = TRUE, only_filter_if_any_vdj_hits = TRUE) {

  stopifnot("count_TRA" %in% colnames(colData(sce)))

  min_tra_filter = colData(sce)[, "count_TRA"] >= tra_range[1]
  max_tra_filter = colData(sce)[, "count_TRA"] <= tra_range[2]

  min_trb_filter = colData(sce)[, "count_TRB"] >= trb_range[1]
  max_trb_filter = colData(sce)[, "count_TRB"] <= trb_range[2]

  min_igl_filter = colData(sce)[, "count_IGL"] >= igl_range[1]
  max_igl_filter = colData(sce)[, "count_IGL"] <= igl_range[2]

  min_igh_filter = colData(sce)[, "count_IGH"] >= igh_range[1]
  max_igh_filter = colData(sce)[, "count_IGH"] <= igh_range[2]

  filter_list = list(min_tra_filter,
                     max_tra_filter,
                     min_trb_filter,
                     max_trb_filter,
                     min_igl_filter,
                     max_igl_filter,
                     min_igh_filter,
                     max_igh_filter)

  if (!tcr_and_bcr_allowed) {
    tr_and_ig_filter = !((colData(sce)$count_TRA + colData(sce)$count_TRB) >= 1 &
                           (colData(sce)$count_IGL + colData(sce)$count_IGH) >= 1)
    filter_list = c(filter_list, list(tr_and_ig_filter))
  }

  result = Reduce("&", filter_list)

  if (only_filter_if_any_vdj_hits) {
    result = result | rowSums(as.matrix(colData(sce)[, c("count_TRA", "count_TRB", "count_IGL", "count_IGH")])) < 1
  }

  return(result)
}