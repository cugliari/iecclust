\name{epclust-package}
\alias{epclust-package}
\alias{epclust}
\docType{package}

\title{
	\packageTitle{epclust}
}

\description{
	\packageDescription{epclust}
}

\details{
  Non-R-base dependencies:
  \itemize{
    \item cluster: for PAM algorithm
    \item bigmemory: to share (big) matrices between processes
    \item wavelets: to compute curves contributions using DWT
    \item Rwave: to compute the CWT
  }

  Suggested packages:
  \itemize{
    \item synchronicity: to compute synchrones in // (concurrent writes)
    \item devtools,testthat,roxygen2: development environment (doc, reload, test...)
    \item MASS: generate multivariate gaussian samples in tests
    \item wmtsa: generate sample curves for \code{claws} examples
    \item DBI: for the example with series in an SQLite DB
    \item digest: to compare \code{claws} examples results
  }

  The package vignette was generated with Jupyter, outside R packaging flow.
}

\author{
	\packageAuthor{epclust}

	Maintainer: \packageMaintainer{epclust}
}
