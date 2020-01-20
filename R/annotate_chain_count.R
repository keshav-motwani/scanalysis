#' Annotate number of TRA, TRB, IGL, and IGH reads per barcode
#'
#' Columns count_TRA, count_TRB, count_IGL, count_IGH, and count_TRA_TRB_IGL_IGH (underscore delimited concatenation of previous 4 columns) are added to colData(sce)
#'
#' @param sce SingleCellExperiment object with a column named vdj containing contig_annotations in colData
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom dplyr as_tibble count
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#'
#' @return SingleCellExperiment object with new count columns in colData
#' @export
#'
#' @examples
#' NULL
annotate_chain_count = function(sce) {
  if ("vdj" %in% colnames(colData(sce))) {

    chain_count = unlist(sce$vdj) %>%
      as_tibble() %>%
      count(barcode, chain) %>%
      pivot_wider(id_cols = barcode, values_from = n, names_from = chain)

    chain_count = column_to_rownames(chain_count, "barcode")[colData(sce)$Barcode, ]

    chain_types = c("TRA", "TRB", "IGL", "IGH")
    for (type in chain_types) {
      if (!type %in% colnames(chain_count)) {
        chain_count[, type] = 0
      }
      chain_count[, type][is.na(chain_count[, type])] = 0
    }

    colnames(chain_count) = paste0("count_", colnames(chain_count))
    chain_count = chain_count[, paste0("count_", chain_types)]
    chain_count$count_TRA_TRB_IGL_IGH = paste0(chain_count$count_TRA, "_", chain_count$count_TRB, "_", chain_count$count_IGL, "_", chain_count$count_IGH)

    colData(sce) = cbind(colData(sce), chain_count)
  }
  return(sce)
}
