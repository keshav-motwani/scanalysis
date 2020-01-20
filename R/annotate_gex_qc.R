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
annotate_total_umi_count = function(sce) {
  data = assays(sce)$counts
  total_umi = colSums(data)
  colData(sce)$total_umi = total_umi
  return(sce)
}

#' Add number of genes expressed per barcode to colData(sce)$n_genes_expr
#'
#' @param sce
#'
#' @importFrom SummarizedExperiment assays colData<-
#'
#' @return
#' @export
#'
#' @examples
#' NULL
annotate_n_genes_expr = function(sce) {
  data = assays(sce)$counts
  n_genes_expr = diff(data@p)
  colData(sce)$n_genes_expr = n_genes_expr
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
annotate_pct_gene_set = function(sce, pattern, name) {
  genes = grep(pattern, rownames(sce), value = TRUE)

  data = assays(sce)$counts
  total_umi = colSums(data)
  data = data[genes, ]
  gene_set_umi = colSums(data)

  pct_gene_set = gene_set_umi/total_umi * 100

  colData(sce)[, name] = pct_gene_set

  return(sce)
}

#' Add number of cells with nonzero expression for each gene to rowData(sce)$num_cells_expr
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
annotate_num_cells_expr = function(sce) {
  data = assays(sce)$counts
  num_cells = rowSums(data != 0)
  rowData(sce)$num_cells_expr = num_cells
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
annotate_pct_total_reads = function(sce) {
  data = assays(sce)$counts
  reads = rowSums(data)
  pct_reads = reads/sum(reads) * 100
  rowData(sce)$pct_reads = pct_reads
  return(sce)
}
