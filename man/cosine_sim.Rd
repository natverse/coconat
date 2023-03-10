% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cosine.R
\name{cosine_sim}
\alias{cosine_sim}
\title{Efficient calculation of cosine similarity for sparse or dense matrices}
\usage{
cosine_sim(x, sparse = FALSE, transpose = FALSE)
}
\arguments{
\item{x}{A data matrix suitable for clustering}

\item{sparse}{Whether to return a sparse or dense matrix (default dense)}

\item{transpose}{When \code{F} (the default) calculates the cosine distance
between columns. When \code{T} calculates the distance between rows.}
}
\value{
A square matrix, dense unless \code{sparse=TRUE} and \code{x} is
  sparse.
}
\description{
Efficient calculation of cosine similarity for sparse or dense matrices
}
\details{
This is much faster than e.g. \code{lsa::cosine} with the added
  benefit that it works for sparse input matrices
}
\examples{
da2ds15=readRDS(system.file('sampledata/da2ds15.rds', package = 'coconat'))
am=partner_summary2adjacency_matrix(da2ds15, inputcol = 'partner', outputcol = 'bodyid')
cosine_sim(am)

\dontrun{
fam_pnkc2=flywire_adjacency_matrix2(
  inputids = 'class:ALPN_R', outputids = 'class:Kenyon_Cell_R',
  sparse=T , threshold = 2)
kckc.cos=cosine_sim(fam_pnkc2)
pnpn.cos=cosine_sim(fam_pnkc2, transpose=T)
}
}
