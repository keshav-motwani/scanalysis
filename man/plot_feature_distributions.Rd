% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_feature_distributions.R
\name{plot_feature_distributions}
\alias{plot_feature_distributions}
\alias{plot_features}
\title{Plot feature distributions from SingleCellExperiment objects}
\usage{
plot_feature_distributions(
  sce_list,
  features,
  x,
  color = x,
  fill = NULL,
  facet_rows = c(".feature"),
  facet_columns = c(".sample"),
  assay = "logcounts",
  alt_exp = NULL,
  ...
)

plot_features(
  sce_list,
  features,
  x,
  color = x,
  fill = NULL,
  facet_rows = c(".feature"),
  facet_columns = c(".sample"),
  assay = "logcounts",
  alt_exp = NULL,
  ...
)
}
\arguments{
\item{sce_list}{list of SingleCellExperiment objects}

\item{features}{numeric features to plot - can be from reducedDims, colData, or assay data}

\item{x}{colData variable to plot on x-axis}

\item{color}{colData variable to color points by}

\item{fill}{colData variable to fill by}

\item{facet_rows}{variables from colData to facet on, can also include ".sample" or ".feature" as described below}

\item{facet_columns}{variables from colData to facet on, can also include ".sample" or ".feature" as described below}

\item{assay}{assay to obtain data from (counts, logcounts, etc)}

\item{alt_exp}{alternate experiment to obtain data from}

\item{...}{other parameters passed to ggexp::plot_distributions}
}
\value{
ggplot object
}
\description{
Plot feature distributions from SingleCellExperiment objects
}
\examples{
NULL
}
