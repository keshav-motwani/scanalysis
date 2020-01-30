#' Get modified version of PBMC 10K data from 10X
#'
#' @importFrom BiocFileCache BiocFileCache
#' @importFrom Matrix colSums
#' @importFrom SummarizedExperiment assay
#' @importFrom purrr map
#'
#' @return
#' @export
#'
#' @examples
#' data = get_pbmc_10k_data()
get_pbmc_10k_data = function() {
  bfc = BiocFileCache("raw_data", ask = FALSE)
  link = "http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_raw_feature_bc_matrix.tar.gz"
  raw.path = bfcrpath(bfc, link)
  untar(raw.path, exdir = file.path("raw_data/pbmc10k"))

  pbmc_10k = read_10x("raw_data/pbmc10k/raw_feature_bc_matrix/")

  subset_data = function(indices) {
    sce = pbmc_10k[, indices]
    sce$timepoint = sample(paste0("Week ", 0:2), length(indices), replace = TRUE)
    sce[, colSums(assay(sce)) > 0]
  }

  dataset = map(split(1:ncol(pbmc_10k), paste0("patient_", 1:4)), subset_data)

  return(dataset)
}