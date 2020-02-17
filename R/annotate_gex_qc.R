#' Add number of total umis per barcode to colData(sce)$total_umi
#'
#' @param sce SingleCellExperiment object
#'
#' @importFrom SummarizedExperiment assays colData<-
#' @importFrom Matrix colSums
#'
#' @return
#' @export
#'
#' @examples
#' NULL
annotate_total_umi_count = function(sce, name = "total_umi", alt_exp = NULL) {
  data = get_assay_data(sce, "counts", alt_exp)
  total_umi = colSums(data)
  colData(sce)[, name] = total_umi
  return(sce)
}

#' Add number of genes expressed per barcode to colData(sce)$n_genes_expr
#'
#' @param sce SingleCellExperiment object
#'
#' @importFrom SummarizedExperiment assays colData<-
#'
#' @return
#' @export
#'
#' @examples
#' NULL
annotate_n_genes_expr = function(sce, name = "n_genes_expr", alt_exp = NULL) {
  data = get_assay_data(sce, "counts", alt_exp)
  n_genes_expr = diff(data@p)
  colData(sce)[, name] = n_genes_expr
  return(sce)
}

#' Add percent of gene set defined by a regular expression pattern to colData(sce)
#'
#' @param sce SingleCellExperiment object
#' @param pattern Regular expression for gene names to sum up expression
#' @param name Column name to save data in colData(sce)
#'
#' @importFrom SummarizedExperiment assays colData<-
#' @importFrom Matrix colSums
#'
#' @return
#' @export
#'
#' @examples
#' NULL
annotate_pct_gene_set = function(sce, pattern, name, alt_exp = NULL) {
  genes = grep(pattern, rownames(sce), value = TRUE)

  data = get_assay_data(sce, "counts", alt_exp)
  total_umi = colSums(data)
  data = data[genes,]
  gene_set_umi = colSums(data)

  pct_gene_set = gene_set_umi / total_umi * 100

  colData(sce)[, name] = pct_gene_set

  return(sce)
}

#' Add number of cells with nonzero expression for each gene to rowData(sce)$num_cells_expr
#'
#' @param sce SingleCellExperiment object
#'
#' @importFrom SummarizedExperiment assays rowData<-
#' @importFrom SingleCellExperiment altExp altExp<-
#' @importFrom Matrix rowSums
#'
#' @return
#' @export
#'
#' @examples
#' NULL
annotate_n_cells_expr = function(sce, name = "n_cells_expr", alt_exp =  NULL) {
  data = get_assay_data(sce, "counts", alt_exp)
  num_cells = rowSums(data != 0)
  if (is.null(alt_exp)) {
    rowData(sce)[, name] = num_cells
  } else {
    rowData(altExp(sce, alt_exp))[, name] = num_cells
  }
  return(sce)
}

#' Add percentage of total reads that each gene takes up across the whole dataset to rowData(sce)$pct_reads
#'
#' @param sce SingleCellExperiment object
#'
#' @importFrom SummarizedExperiment assays rowData<-
#' @importFrom Matrix rowSums
#'
#' @return
#' @export
#'
#' @examples
#' NULL
annotate_pct_total_reads = function(sce, name = "pct_reads", alt_exp = NULL) {
  data = get_assay_data(sce, "counts", alt_exp)
  reads = rowSums(data)
  pct_reads = reads / sum(reads) * 100
  if (is.null(alt_exp)) {
    rowData(sce)[, name] = pct_reads
  } else {
    rowData(altExp(sce, alt_exp))[, name] = pct_reads
  }
  return(sce)
}

#' Add average of ambient feature reads to rowData
#'
#' @param sce SingleCellExperiment object
#'
#' @importFrom SummarizedExperiment assays rowData<- colData<-
#' @importFrom SingleCellExperiment altExp
#' @importFrom Matrix rowMeans Matrix
#' @importFrom qlcMatrix corSparse
#'
#' @return
#' @export
#'
#' @examples
#' NULL
annotate_ambient_profile = function(sce,
                                    ambient_indices,
                                    name = "ambient_profile",
                                    alt_exp = NULL) {
  data = get_assay_data(sce, "counts", alt_exp)[, ambient_indices]
  ambient = rowMeans(data)
  if (is.null(alt_exp)) {
    rowData(sce)[, name] = ambient
  } else {
    rowData(altExp(sce, alt_exp))[, name] = ambient
  }
  pearson_corr = corSparse(data, Matrix(matrix(ambient), sparse = TRUE))
  colData(sce)[, paste0(name, "_corr")] = pearson_corr
  return(sce)
}
