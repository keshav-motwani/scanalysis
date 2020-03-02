#' Annotate number of occurrences of a given clonotype in a sample
#'
#' @param sce SingleCellExperiment object
#' @param clonotype_name name of column in colData containing clonotype - can be generated using `assign_clonotypes`.
#' @param group_by columns to group repertoires by - only counts clonotypes within groups
#'
#' @importFrom dplyr as_tibble filter group_by arrange left_join desc
#' @importFrom SummarizedExperiment colData<- colData
#'
#' @return SingleCellExperiment object with added columns to colData including "clonotype", "clonotype_count", and "clonotype_count_2:10"
#' @export
#'
#' @examples
#' NULL
annotate_clonotype_count = function(sce,
                                    clonotype_name = "clonotype",
                                    group_by = c()) {

  vdj = as.data.frame(colData(sce)[, c(clonotype_name, group_by), drop = FALSE])

  counts = vdj %>%
    group_by(.dots = c(group_by, clonotype_name)) %>%
    count(name = "clonotype_count") %>%
    arrange(desc(clonotype_count)) %>%
    filter(!is.na(!!as.name(clonotype_name)))

  vdj = left_join(vdj, counts, by = c(clonotype_name, group_by))

  for (i in 2:max(counts$clonotype_count)) {
    vdj[, paste0("clonotype_count_", i)] = sprintf("< %s counts", i)
    vdj[which(vdj$clonotype_count >= i), paste0("clonotype_count_", i)] = vdj[which(vdj$clonotype_count >= i), "clonotype", drop = TRUE]
  }

  colnames(vdj) = gsub("clonotype", clonotype_name, colnames(vdj))

  colData(sce) = cbind(colData(sce), vdj)

  return(sce)
}
