#' Get modified version of PBMC 10k data from 10X
#'
#' @importFrom BiocFileCache BiocFileCache bfcrpath
#' @importFrom Matrix colSums
#' @importFrom SummarizedExperiment assay
#' @importFrom purrr map
#'
#' @return
#' @export
#'
#' @examples
#' data = get_pbmc_10k_data()
get_multi_sample_pbmc_10k = function() {
  bfc = BiocFileCache("example_data", ask = FALSE)

  link = "http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_raw_feature_bc_matrix.tar.gz"
  raw.path = bfcrpath(bfc, link)

  untar(raw.path, exdir = file.path("example_data/multi_sample_pbmc_10k"))

  pbmc_10k = read_10x("example_data/multi_sample_pbmc_10k/raw_feature_bc_matrix/")

  subset_data = function(indices) {
    sce = pbmc_10k[, indices]
    sce$timepoint = sample(paste0("Week ", 0:2), length(indices), replace = TRUE)
    sce[, colSums(assay(sce)) > 0]
  }

  dataset = map(split(1:ncol(pbmc_10k), paste0("patient_", 1:4)), subset_data)

  return(dataset)
}

#' Get PBMC 5k data using v3 chemistry
#'
#' @importFrom BiocFileCache BiocFileCache bfcrpath
#' @importFrom Matrix colSums
#' @importFrom SummarizedExperiment assay
#' @importFrom purrr map
#'
#' @return
#' @export
#'
#' @examples
#' data = get_pbmc_5k_v3()
get_pbmc_5k_v3 = function() {
  bfc = BiocFileCache("example_data", ask = FALSE)

  link = "http://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_filtered_feature_bc_matrix.tar.gz"
  raw.path = bfcrpath(bfc, link)

  untar(raw.path, exdir = file.path("example_data/pbmc_5k_v3"))

  pbmc_5k_v3 = read_10x("example_data/pbmc_5k_v3/filtered_feature_bc_matrix/")

  return(pbmc_5k_v3)
}

#' Get PBMC 5k data using v3 chemistry
#'
#' @importFrom BiocFileCache BiocFileCache bfcrpath
#' @importFrom Matrix colSums
#' @importFrom SummarizedExperiment assay
#' @importFrom purrr map
#'
#' @return
#' @export
#'
#' @examples
#' data = get_pbmc_5k_nextgem()
get_pbmc_5k_nextgem = function() {
  bfc = BiocFileCache("example_data", ask = FALSE)

  link = "http://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3_nextgem/5k_pbmc_v3_nextgem_filtered_feature_bc_matrix.tar.gz"
  raw.path = bfcrpath(bfc, link)

  untar(raw.path, exdir = file.path("example_data/pbmc_5k_nextgem"))

  pbmc_5k_nextgem = read_10x("example_data/pbmc_5k_nextgem/filtered_feature_bc_matrix/")

  return(pbmc_5k_nextgem)
}