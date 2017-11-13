#' @include de_serialize.R
#' @include clustering.R
#' @include main.R
#'
#' @useDynLib epclust
#'
#' @importFrom Rwave cwt
#' @importFrom cluster pam
#' @importFrom parallel makeCluster clusterExport parLapply stopCluster
#' @importFrom wavethresh wd filter.select
#' @importFrom stats spline
#' @importFrom methods is
#' @importFrom bigmemory big.matrix as.big.matrix is.big.matrix
#' @importFrom utils head tail
NULL
