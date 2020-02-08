#' Read 10X gene expression and VDJ data into a SingleCellExperiment object
#'
#' VDJ data is stored in colData(sce)$vdj for now. This makes it easy to manipulate using the tidyverse for assignment of clonotypes, etc.
#'
#' @param gene_expr_path Path to gene expression data - can be full path to the cellranger count output for the sample, or the output path of cellranger count if sample name is specified
#' @param vdj_path Path to vdj data - can be full path to the cellranger vdj output for the sample, or the output path of cellranger vdj if sample name is specified
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @importMethodsFrom S4Vectors split
#'
#' @return
#' @export
#'
#' @examples
#' NULL
read_10x = function(gene_expr_path = getwd(),
                    vdj_path = NA) {

  sce = .read_10x_gene_expr(gene_expr_path)

  if (!is.na(vdj_path)) {
    vdj = .read_10x_vdj(vdj_path)
    colData(sce)$vdj = split(vdj, factor(vdj$barcode, colData(sce)$Barcode))
  }

  return(sce)
}


#' Read gene expression matrix from 10X data
#'
#' @param path Path to gene expression data - can be full path to the cellranger count output for the sample, or the output path of cellranger count if sample name is specified
#'
#' @importFrom DropletUtils read10xCounts
#' @importFrom scater uniquifyFeatureNames
#' @importFrom SingleCellExperiment rowData splitAltExps
#'
#' @return SingleCellExperiment object
#' @keywords internal
#'
#' @examples
#' NULL
.read_10x_gene_expr = function(path) {

  sce = read10xCounts(path, col.names = TRUE) %>%
    splitAltExps(ifelse(grepl("^HTO_", rownames(.)), "HTO", "gene")) %>%
    splitAltExps(ifelse(rowData(.)$Type ==  "Antibody Capture", "ADT", "gene"))

  rownames(sce) = uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)

  metadata(sce)$default_assay = "RNA"

  return(sce)
}


#' Read VDJ information from 10X data
#'
#' @param path Path to vdj data - can be full path to the cellranger vdj output for the sample, or the output path of cellranger vdj if sample name is specified
#'
#' @importFrom readr read_csv
#' @importFrom dplyr filter
#' @importClassesFrom S4Vectors DataFrame
#'
#' @return DataFrame with vdj information
#' @keywords internal
#'
#' @examples
#' NULL
.read_10x_vdj = function(path) {

  vdj = read_csv(file.path(path, "all_contig_annotations.csv")) %>%
    filter(
      high_confidence,
      is_cell,
      productive == "True" | productive,
      chain %in% c("TRA", "TRB", "IGL", "IGK", "IGH")
    ) %>%
    as("DataFrame")

  return(vdj)
}
