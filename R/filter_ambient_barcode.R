#' Get filter to remove ambient RNA barcodes based on the combination of multiple methods:
#'
#' inflection - all barcodes below the inflection point in the UMI rank vs total UMI plot are called ambient
#' knee - all barcodes below the knee in the UMI rank vs total UMI plot are called ambient
#' cellranger (v2) - all barcodes below 10% of the 99th percentile of total UMIs in the top n_cells barcodes are called ambient
#' empty_drops - uses the algorithm from DropletUtils::emptyDrops to compare each barcode to the profile of barcodes below
#'
#' @param sce SingleCellExperiment object
#' @param n_cells Number of input cells (used in cellranger cutoff computation)
#' @param fdr FDR threshold for keeping a cell for emptyDrops
#' @param lower_umi_limit Number of total UMIs below which a barcode is assumed to contain ambient RNA, used in emptyDrops method only
#'
#' @return List for each method in containing filter for method and attributes with details specific to each method.
#' @export
#'
#' @examples
#' NULL
filter_ambient_barcode = function(sce, n_cells, fdr = 0.01, lower_umi_limit = 100) {

  filters = list()
  filters$inflection = filter_ambient_barcode_inflection(sce)
  filters$empty_drops = filter_ambient_barcode_empty_drops(sce, fdr, lower_umi_limit)
  filters$knee = filter_ambient_barcode_knee(sce)
  filters$cellranger = filter_ambient_barcode_cellranger(sce, n_cells)

  return(filters)
}

#' Get filter for cells based on inflection point
#'
#' @param sce SingleCellExperiment object
#'
#' @importFrom DropletUtils barcodeRanks
#'
#' @return Boolean filter with TRUE for cells and FALSE for ambient RNA. Attribute named value containing total UMI at inflection point.
#' @export
#'
#' @examples
#' NULL
filter_ambient_barcode_inflection = function(sce) {
  bcrank = barcodeRanks(counts(sce))
  inflection = bcrank@metadata$inflection
  filter = bcrank$total > inflection
  attributes(filter)$value = inflection
  return(filter)
}

#' Get filter for cells based on emptyDrops algorithm
#'
#' @param sce SingleCellExperiment object
#' @param fdr Upper limit for FDR
#' @param lower_limit_umi Lower limit for UMI, below which all barcodes are assumed to be ambient
#'
#' @importFrom DropletUtils emptyDrops
#' @importFrom dplyr as_tibble mutate
#'
#' @return Boolean filter with TRUE for cells and FALSE for ambient RNA. Attribute named data containing output of emptyDrops function
#' @export
#'
#' @examples
#' NULL
filter_ambient_barcode_empty_drops = function(sce, fdr, lower_umi_limit) {
  out = emptyDrops(counts(sce), lower_umi_limit) %>%
    as_tibble() %>%
    mutate(cell = FDR <= fdr) %>%
    mutate(cell = ifelse(cell == TRUE, "cell", "ambient")) %>%
    mutate(cell = ifelse(is.na(cell), "< min UMI", cell))
  filter = out$cell == "cell"
  attributes(filter)$data = out
  return(filter)
}

#' Get filter for cells based on knee point
#'
#' @param sce SingleCellExperiment object
#'
#' @importFrom DropletUtils barcodeRanks
#'
#' @return Boolean filter with TRUE for cells and FALSE for ambient RNA. Attribute named value containing total UMI at knee point.
#' @export
#'
#' @examples
#' NULL
filter_ambient_barcode_knee = function(sce) {
  bcrank = barcodeRanks(counts(sce))
  knee = bcrank@metadata$knee
  filter = bcrank$total > knee
  attributes(filter)$value = knee
  return(filter)
}

#' Get filter for cells based on inflection point
#'
#' @param sce SingleCellExperiment object
#' @param n_cells Number of cells expected
#'
#' @importFrom DropletUtils barcodeRanks
#'
#' @return Boolean filter with TRUE for cells and FALSE for ambient RNA. Attribute named value containing total UMI at the cutoff point.
#' @export
#'
#' @examples
#' NULL
filter_ambient_barcode_cellranger = function(sce, n_cells) {
  bcrank = barcodeRanks(counts(sce))
  cutoff = quantile(sort(bcrank$total, decreasing = TRUE)[1:n_cells], 0.99, na.rm = TRUE)/10
  filter = bcrank$total > cutoff
  attributes(filter)$value = cutoff
  return(filter)
}
