#' Get filter for total umis per barcode based on log transformed values outside (either above, below, or both based on type parameter) nmads median absolute deviations from the median
#'
#' Wrapper around scater::isOutlier, arguments type and nmads refer to the same thing. For total_umi, uses log transformed values for deciding cutoffs.
#'
#' @param sce SingleCellExperiment object
#' @param nmads number of median absolute deviations from the median for cutoff
#' @param type either "both", "lower", or "higher", referring to which side to filter on
#' @param column column in colData containing total UMI count
#'
#' @import SingleCellExperiment
#' @importFrom scater isOutlier
#'
#' @return Boolean filter with TRUE for cells and FALSE for outliers. Contains attribute with number of median absolute deviations specified.
#' @export
#'
#' @examples
#' NULL
filter_total_umi = function(sce, nmads, type, column = "total_umi") {
  not_outlier = !isOutlier(colData(sce)[, column],
                           nmads = nmads,
                           type = type,
                           log = TRUE)
  attributes(not_outlier)$nmads = nmads
  return(not_outlier)
}

#' Get filter for number of genes expressed per barcode based on log transformed values
#'
#' Wrapper around scater::isOutlier, arguments type and nmads refer to the same thing. For number of genes expressed, uses log transformed values for deciding cutoffs.
#'
#' @param sce SingleCellExperiment object
#' @param nmads Number of median absolute deviations from the median for cutoff
#' @param type Either "both", "lower", or "higher", referring to which side to filter on
#' @param column column in colData containing number of genes expressed
#'
#' @import SingleCellExperiment
#' @importFrom scater isOutlier
#'
#' @return Boolean filter with TRUE for cells and FALSE for outliers. Contains attribute with number of median absolute deviations specified.
#' @export
#'
#' @examples
#' NULL
filter_n_genes_expr = function(sce, nmads, type, column = "n_genes_expr") {
  not_outlier = !isOutlier(colData(sce)[, column],
                           nmads = nmads,
                           type = type,
                           log = TRUE)
  attributes(not_outlier)$nmads = nmads
  return(not_outlier)
}

#' Get filter for percentage of mitochondrial reads expressed per barcode
#'
#' Wrapper around scater::isOutlier, arguments type and nmads refer to the same thing. For mitochondrial reads, uses non-log transformed values for deciding cutoffs.
#'
#' Ensure that \code{\link{add_pct_gene}} has been run before this or that a column named pct_mito is present in colData(sce)
#'
#' @param sce SingleCellExperiment object
#' @param nmads number of median absolute deviations from the median
#' @param type either "both", "lower", or "higher", referring to which side to filter on
#' @param column column in colData containing percentage of mitochondrial reads
#'
#' @import SingleCellExperiment
#' @importFrom scater isOutlier
#'
#' @return Boolean filter with TRUE for cells and FALSE for outliers. Contains attribute with number of median absolute deviations specified.
#' @export
#'
#' @examples
#' NULL
filter_pct_mito = function(sce, nmads, type, column = "pct_mito") {
  not_outlier = !isOutlier(colData(sce)[, column],
                           nmads = nmads,
                           type = type)
  attributes(not_outlier)$nmads = nmads
  return(not_outlier)
}
