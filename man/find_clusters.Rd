% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_clusters.R
\name{find_clusters}
\alias{find_clusters}
\title{Find cell clusters}
\usage{
find_clusters(
  sce,
  use_dimred = "PCA",
  seed = 100,
  snn_k = 10,
  ncores = 1,
  method = c("walktrap", "spinglass"),
  steps = 4,
  spins = 25,
  min_member = 20,
  prefix = NULL,
  plot = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{sce}{A SingleCellExperiment object containing expression values, usually counts.}

\item{use_dimred}{A string specifying whether existing values in \code{reducedDims(sce)} should be used.}

\item{seed}{Random seed.}

\item{snn_k}{The number of nearest neighbors to consider during graph construction.}

\item{ncores}{Number of cores.}

\item{method}{"walktrap" or "spinglass" for finding communities in graphs via short random walks or a spin-glass model
and simulated annealing.}

\item{steps}{The length of the random walks to perform.}

\item{spins}{Integer constant, the number of spins to use. This is the upper
limit for the number of communities. It is not a problem to supply a
(reasonably) big number here, in which case some spin states will be
unpopulated.}

\item{min_member}{Minimal number of cluster members.}

\item{prefix}{Prefix for file name for the QC metrics histograms.}

\item{plot}{TRUE/FALSE for whether to plot the QC metrics histograms.}

\item{verbose}{TRUE/FALSE for specifying if diagnostics should be printed to screen.}
}
\value{
A SingleCellExperiment object with cell cluster information.
}
\description{
Find cell clusters using \pkg{igraph}.
}
