#' QC plot of UMI rank vs total number of UMIs per barcode
#'
#' @param sce_list list of SingleCellExperiment objects
#' @param ambient_rna_filters list of results from filter_ambient_rna function, one for each SCE object in sce_list - must be the same names and order
#' @param color column from colData that is in all objects in sce_list
#' @param shape column from colData that is in all objects in sce_list
#' @param text_size font size for annotations
#' @param facet_rows columns to facet on
#' @param facet_columns columns to facet on
#' @param facet_type either "wrap" or "grid", same as ggplot
#' @param ... other params passed into either facet_wrap or facet_grid, depending on facet_type parameter
#'
#' @import ggplot2
#' @importFrom DropletUtils barcodeRanks
#' @importFrom purrr map_chr pmap_dfr
#' @importFrom tibble enframe
#' @importFrom dplyr arrange filter
#' @importFrom ggexp theme_ggexp
#' @importFrom scales trans_breaks trans_format math_format
#' @importFrom ggexp plot_facets theme_ggexp
#'
#' @return ggplot with plot of UMI rank vs total number of UMIs per barcode
#' @export
#'
#' @examples
#' NULL
plot_barcode_qc = function(sce_list,
                           ambient_rna_filters,
                           color = "empty_drops",
                           shape = NULL,
                           facet_rows = NULL,
                           facet_columns = NULL,
                           facet_type = "grid",
                           point_size = 1,
                           text_size = 2,
                           alpha = 0.5,
                           ...) {

  stopifnot(names(sce_list) == names(ambient_rna_filters))

  data = pmap_dfr(list(sce_list, names(sce_list), ambient_rna_filters),
                  ~ .prepare_barcode_plot_data(..1, ..2, ..3, c(color, shape, facet_rows, facet_columns)))

  hlines = pmap_dfr(list(sce_list, names(sce_list), ambient_rna_filters),
                    .prepare_barcode_cutoffs)

  counts = pmap_dfr(
    list(sce_list, names(sce_list), ambient_rna_filters),
    ~ .prepare_barcode_counts(..1, ..2, ..3, setdiff(
      c(facet_rows, facet_columns), ".sample"
    ))
  )

  rank_vs_umi = ggplot(mapping = aes(rank, total_umi)) +
    geom_point(
      data = data,
      aes_string(color = color, shape = shape),
      size = point_size,
      alpha = alpha
    ) +
    geom_hline(data = hlines, aes_string(yintercept = "y")) +
    scale_x_log10(
      breaks = trans_breaks("log10", function(x)
        10 ^ x),
      labels = trans_format("log10", math_format(10 ^ .x))
    ) +
    scale_y_log10(
      breaks = trans_breaks("log10", function(x)
        10 ^ x),
      labels = trans_format("log10", math_format(10 ^ .x))
    ) +
    theme_ggexp() +
    labs(y = "Total UMI",
         x = "UMI Rank") +
    geom_text(
      data = hlines,
      aes(
        x = Inf,
        y = y,
        label = paste0(" ", color)
      ),
      hjust = 1,
      vjust = -0.5,
      size = 2
    ) +
    geom_text(
      data = counts,
      aes(x = 0, y = 0, label = label),
      hjust = 0,
      vjust = 0,
      size = text_size
    )

  rank_vs_umi = plot_facets(rank_vs_umi,
                            facet_rows,
                            facet_columns,
                            facet_type,
                            ...)

  return(rank_vs_umi)
}

#' Prepare data for plotting barcode QC plot
#'
#' @param sce SingleCellExperiment object
#' @param sample_name name of sample
#' @param ambient_rna_filter result from filter_ambient_rna function
#'
#' @importFrom DropletUtils barcodeRanks
#' @importFrom dplyr filter arrange
#' @importFrom SingleCellExperiment colData
#'
#' @return
#' @keywords internal
#'
#' @examples
#' NULL
.prepare_barcode_plot_data = function(sce, sample_name, ambient_rna_filter, metadata) {

  bcrank = barcodeRanks(counts(sce))

  rank_total =
    data.frame(
      rank = bcrank$rank,
      total_umi = bcrank$total,
      empty_drops = attributes(ambient_rna_filter$empty_drops)$data$cell,
      .sample = sample_name
    )

  rank_total = cbind(rank_total, colData(sce)[, intersect(metadata, colnames(colData(sce))), drop = FALSE])
  rank_total = as.data.frame(rank_total[, unique(colnames(rank_total)), drop = FALSE])

  return(rank_total)
}

#' Prepare cutoff values for each method
#'
#' @param sce SingleCellExperiment object
#' @param sample_name name of sample
#' @param ambient_rna_filter result from filter_ambient_rna function
#'
#' @return
#' @keywords internal
#'
#' @examples
#' NULL
.prepare_barcode_cutoffs = function(sce, sample_name, ambient_rna_filter) {

  hlines = data.frame(
    y = c(
      attributes(ambient_rna_filter$knee)$value,
      attributes(ambient_rna_filter$inflection)$value,
      attributes(ambient_rna_filter$cellranger)$value
    ),
    color = c("knee", "inflection", "cellranger"),
    .sample = sample_name
  )

  return(hlines)

}


#' Prepare count annotations for number of cells passing each filtering method
#'
#' @param sce SingleCellExperiment object
#' @param sample_name name of sample
#' @param ambient_rna_filter result from filter_ambient_rna function
#' @param facets columns to group by in plot
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom dplyr summarize group_by mutate
#'
#' @return
#' @keywords internal
#'
#' @examples
#' NULL
.prepare_barcode_counts = function(sce,
                                   sample_name,
                                   ambient_rna_filter,
                                   facets) {

  as.data.frame(cbind(data.frame(ambient_rna_filter), colData(sce)[, facets, drop = FALSE])) %>%
    group_by(.dots = facets) %>%
    summarize(
      label = paste0(
        " knee: ",
        sum(knee, na.rm = TRUE),
        "\n",
        " inflection: ",
        sum(inflection, na.rm = TRUE),
        "\n",
        " cellranger: ",
        sum(cellranger, na.rm = TRUE),
        "\n",
        " empty_drops: ",
        sum(empty_drops, na.rm = TRUE),
        "\n",
        " total: ",
        length(empty_drops),
        "\n"
      )
    ) %>%
    mutate(.sample = sample_name)

}
