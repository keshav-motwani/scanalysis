#' Annotate CDR3 sequences stored in colData(sce)$vdj based on a reference dataset (for example, Ag-specific data)
#'
#' @param sce SingleCellExperiment object
#' @param reference reference dataset for use in annotating query sequences
#' @param reference_cdr3_column column in reference which contains the CDR3 region
#' @param reference_annotation_column column in reference to annotate onto query dataset
#' @param chain_match chain(s) to match (for example, c("TRA", "TRB"))
#' @param max_dist maximum distance to allow between query and reference sequence
#' @param method distance metric- one of "osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex"
#'
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom dplyr filter group_by count mutate select left_join
#' @importFrom fuzzyjoin stringdist_left_join
#'
#' @return SingleCellExperiment object with additional column with the same name as reference_annotation_column
#' @export
#'
#' @examples
#' NULL
annotate_cdr3 = function(sce, reference, reference_cdr3_column, reference_annotation_column, chain_match = "TRB", max_dist = 1, method = "hamming") {
  vdj_data_all = unnest_vdj(sce, include_all = FALSE)
  vdj_data = vdj_data_all %>%
    filter(chain == chain_match)
  reference = reference[, c(reference_cdr3_column, reference_annotation_column), drop = FALSE]
  reference = reference[!duplicated(reference), ]
  result = stringdist_left_join(vdj_data, reference, by = c(cdr3 = reference_cdr3_column), max_dist = max_dist, method = method)[, c("Barcode", "cdr3", reference_cdr3_column, reference_annotation_column)]
  result = result %>%
    filter(!is.na(!!as.name(reference_annotation_column))) %>%
    group_by(Barcode) %>%
    count(!!as.name(reference_annotation_column)) %>%
    mutate(the_rank  = rank(-n, ties.method = "random")) %>%
    filter(the_rank == 1) %>%
    select(-the_rank)
  vdj_data = left_join(vdj_data_all, result, by = "Barcode")
  vdj_data = as(vdj_data, "DataFrame")
  colData(sce)$vdj = split(vdj_data, factor(vdj_data$Barcode, colData(sce)$Barcode))
  return(sce)
}
