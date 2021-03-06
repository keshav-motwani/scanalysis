% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_evenness_profile.R
\name{compute_evenness_profile_matrix}
\alias{compute_evenness_profile_matrix}
\title{Compute matrix of evenness profiles per group based on clonotype frequency distributions}
\usage{
compute_evenness_profile_matrix(
  sce_list,
  group_by = c(),
  min_alpha = 0,
  max_alpha = 10,
  alpha_steps = 51
)
}
\arguments{
\item{sce_list}{list of SingleCellExperiment objects}

\item{group_by}{column names from colData to use to construct groups for clonotype counting}

\item{min_alpha}{minimum Renyi entropy alpha}

\item{max_alpha}{maximum Renyi entropy alpha}

\item{alpha_steps}{number of steps in between min_alpha and max_alpha}
}
\value{

}
\description{
Compute matrix of evenness profiles per group based on clonotype frequency distributions
}
\examples{
NULL
}
