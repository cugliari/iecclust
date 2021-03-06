% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{clustering}
\alias{clustering}
\alias{clusteringTask1}
\alias{clusteringTask2}
\alias{clusteringTask1}
\alias{clusteringTask2}
\title{Two-stage clustering, within one task (see \code{claws()})}
\usage{
clusteringTask1(indices, getContribs, K1, algoClust1, nb_items_clust,
  ncores_clust = 3, verbose = FALSE)

clusteringTask2(indices, getSeries, K2, algoClust2, nb_series_per_chunk,
  smooth_lvl, nvoice, nbytes, endian, ncores_clust = 3, verbose = FALSE)
}
\arguments{
\item{indices}{Indices of the series to consider}

\item{getContribs}{Function to retrieve contributions from initial series indices:
\code{getContribs(indices)} outputs a contributions matrix, in columns}

\item{K1}{Number of clusters to be found after stage 1 (K1 << N)}

\item{algoClust1}{Clustering algorithm for stage 1. A function which takes (data, K)
as argument where data is a matrix in columns and K the desired number of clusters,
and outputs K medoids ranks. Default: PAM.}

\item{nb_items_clust}{Number of items in 1st clustering algorithm input}

\item{ncores_clust}{Number of parallel clusterings in one task}

\item{verbose}{FALSE: nothing printed; TRUE: some execution traces}

\item{getSeries}{Function to retrieve series (argument: 'inds', integer vector),
as columns of a matrix}

\item{K2}{Number of clusters to be found after stage 2 (K2 << K1)}

\item{algoClust2}{Clustering algorithm for stage 2. A function which takes (dists, K)
as argument where dists is a matrix of distances and K the desired number of
clusters, and outputs K medoids ranks. Default: PAM.}

\item{nb_series_per_chunk}{Number of series to retrieve in one batch}

\item{smooth_lvl}{Smoothing level: odd integer, 1 == no smoothing.}

\item{nvoice}{Number of voices within each octave for CWT computations}

\item{nbytes}{4 or 8 bytes to (de)serialize a floating-point number}

\item{endian}{Endianness for (de)serialization: "little" or "big"}
}
\value{
The indices of the computed (resp. K1 and K2) medoids.
}
\description{
\code{clusteringTask1()} runs one full stage-1 task, which consists in iterated
clustering on nb_curves / ntasks energy contributions, computed through
discrete wavelets coefficients.
\code{clusteringTask2()} runs a full stage-2 task, which consists in WER distances
computations between medoids (indices) output from stage 1, before applying
the second clustering algorithm on the distances matrix.
}
