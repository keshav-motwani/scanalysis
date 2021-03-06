% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/encode_identity_frequency.R
\name{encode_cell_identity_frequency_matrix}
\alias{encode_cell_identity_frequency_matrix}
\title{Encode frequencies of (combinations) of values in columns from colData into a matrix with rows as groups and columns as features}
\usage{
encode_cell_identity_frequency_matrix(
  sce_list,
  attributes,
  group_by = c(),
  normalize = "none"
)
}
\arguments{
\item{sce_list}{Llist of SCE objects, all containing relevant attributes in colData as columns}

\item{attributes}{column names from colData which to encode}

\item{group_by}{column names from colData which to group by (compute frequencies per group)}

\item{normalize}{normalization method - either "none" or "relative_frequency"}
}
\value{

}
\description{
Encode frequencies of (combinations) of values in columns from colData into a matrix with rows as groups and columns as features
}
\examples{
NULL
}
