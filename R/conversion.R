#' Convert Seurat object to SingleCellExperiment and retain multi-modal data
#'
#' @param seurat Seurat object
#'
#' @importFrom S4Vectors metadata<-
#' @importFrom SingleCellExperiment altExps
#' @importFrom Seurat DefaultAssay as.SingleCellExperiment VariableFeatures
#'
#' @return
#' @export
#'
#' @examples
#' NULL
seurat_to_sce = function(seurat) {
  if ("SCT" %in% names(seurat@assays)) {
    default_assay = "SCT"
  } else {
    default_assay = "RNA"
  }

  result = seurat_assay_to_sce(seurat, default_assay)

  for (assay in setdiff(names(seurat@assays), default_assay)) {
    value = seurat_assay_to_sce(seurat, assay)
    altExps(result) = c(altExps(result), list(value))
    names(altExps(result))[length(names(altExps(result)))] = assay
  }

  metadata(result)$default_assay = default_assay

  return(result)
}

#' Convert Seurat assay to SingleCellExperiment
#'
#' @param seurat Seurat object
#' @param assay Assay to convert
#'
#' @importFrom S4Vectors metadata<-
#' @importFrom Seurat as.SingleCellExperiment VariableFeatures
#'
#' @return
#' @export
#'
#' @examples
#' NULL
seurat_assay_to_sce = function(seurat, assay) {

  result = as.SingleCellExperiment(seurat, assay = assay)

  if (!is.null(seurat@assays[[assay]]@misc)) metadata(result) = seurat@assays[[assay]]@misc
  metadata(result)$scaled = seurat@assays[[assay]]@scale.data
  metadata(result)$variable_features = VariableFeatures(seurat, assay = assay)

  return(result)
}

#' Convert SingleCellExperiment object to Seurat and retain multi-modal data
#'
#' @param sce SingleCellExperiment object
#'
#' @importFrom S4Vectors metadata
#' @importFrom SingleCellExperiment altExps
#' @importFrom Seurat DefaultAssay<-
#' @importFrom purrr map
#'
#' @return
#' @export
#'
#' @examples
#' NULL
sce_to_seurat = function(sce) {

  main = sce_to_assay(sce, return_assay = FALSE)

  alt_exps = as.list(altExps(sce))
  alt_exps = map(alt_exps, sce_to_assay)

  assays = c(list(main@assays[[1]]), alt_exps)

  if (is.null(metadata(sce)$default_assay)) {
    default_assay = "RNA"
  } else {
    default_assay = metadata(sce)$default_assay
  }

  names(assays) = c(default_assay, names(alt_exps))

  main@assays = assays
  DefaultAssay(main) = default_assay

  return(main)
}

#' Convert SingleCellExperiment object to Seurat assay
#'
#' @param sce SingleCellExperiment object
#'
#' @importFrom S4Vectors metadata
#' @import SingleCellExperiment
#' @importFrom Seurat DefaultAssay as.Seurat VariableFeatures VariableFeatures<-
#' @importFrom purrr map
#'
#' @return
#' @export
#'
#' @examples
#' NULL
sce_to_assay = function(sce, return_assay = TRUE) {
  if ("logcounts" %in% names(assays(sce))) {
    data = "logcounts"
  } else {
    data = NULL
  }
  if ("counts" %in% names(assays(sce))) {
    counts = "counts"
  } else {
    counts = NULL
  }
  seurat = as.Seurat(sce, counts = counts, data = data)

  if (!is.null(metadata(sce)$scaled)) {
    seurat@assays[[1]]@scale.data = metadata(sce)$scaled
  }

  if (return_assay) {
    seurat = seurat@assays[[1]]
    seurat@misc = metadata(sce)
  } else {
    seurat@assays[[1]]@misc = metadata(sce)
  }

  VariableFeatures(seurat) = metadata(sce)[["variable_features"]]

  return(seurat)
}
