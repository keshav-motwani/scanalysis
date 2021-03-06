#' Function to use for assigning clonotypes based on custom definition
#'
#' @param sce SingleCellExperiment object
#' @param clonotype_fields fields to include in clonotype definition (e.g., c("v_gene", "cdr3", "j_gene", etc.))
#' @param chains chains which to include in clonotype definition(e.g., c("TRA", "TRB"))
#' @param tra/trb/igl/igk/igh_range vector of length 2 with minimum and maximum number of unique chains (inclusive)
#'
#' @importFrom dplyr as_tibble filter group_by summarise_all select left_join
#' @importFrom tidyr unite
#' @importFrom tibble column_to_rownames
#' @importFrom SummarizedExperiment colData<- colData
#'
#' @return SingleCellExperiment object with added columns to colData including "clonotype", "clonotype_count", and "clonotype_count_2:10"
#' @export
#'
#' @examples
#' NULL
assign_clonotypes = function(sce,
                             clonotype_fields = c("v_gene", "cdr3", "j_gene"),
                             chains = c("TRA", "TRB"),
                             tra_range = c(0, 2),
                             trb_range = c(0, 1),
                             igl_range = c(0, 2),
                             igk_range = c(0, 2),
                             igh_range = c(0, 1),
                             name = "clonotype") {

  col_data = colData(sce)

  indices = col_data$count_TRA >= tra_range[1] &
    col_data$count_TRB >= trb_range[1] &
    col_data$count_IGL >= igl_range[1] &
    col_data$count_IGK >= igk_range[1] &
    col_data$count_IGH >= igh_range[1] &
    col_data$count_TRA <= tra_range[2] &
    col_data$count_TRB <= trb_range[2] &
    col_data$count_IGL <= igl_range[2] &
    col_data$count_IGK <= igk_range[2] &
    col_data$count_IGH <= igh_range[2]

  vdj = unnest_vdj(sce[, indices], FALSE) %>%
    as_tibble() %>%
    filter(chain %in% chains) %>%
    group_by(Barcode) %>%
    summarise_all(function(x)
      paste(sort(unique(x)), collapse = "_")) %>%
    unite(
      data = .,
      col = clonotype, !!clonotype_fields,
      sep = "///",
      remove = TRUE,
      na.rm = TRUE
    ) %>%
    select(Barcode, clonotype)

  vdj = column_to_rownames(vdj, "Barcode")[colData(sce)$Barcode, , drop = FALSE]

  colnames(vdj) = gsub("clonotype", name, colnames(vdj))

  colData(sce) = cbind(colData(sce), vdj)

  return(sce)
}
