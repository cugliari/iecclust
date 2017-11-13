#' Two-stage clustering, within one task (see \code{claws()})
#'
#' \code{clusteringTask1()} runs one full stage-1 task, which consists in iterated
#' clustering on nb_curves / ntasks energy contributions, computed through
#' discrete wavelets coefficients.
#' \code{clusteringTask2()} runs a full stage-2 task, which consists in WER distances
#' computations between medoids (indices) output from stage 1, before applying
#' the second clustering algorithm on the distances matrix.
#'
#' @param getContribs Function to retrieve contributions from initial series indices:
#'   \code{getContribs(indices)} outputs a contributions matrix, in columns
#' @inheritParams claws
#' @inheritParams computeSynchrones
#' @inheritParams computeWerDists
#'
#' @return The indices of the computed (resp. K1 and K2) medoids.
#'
#' @name clustering
#' @rdname clustering
#' @aliases clusteringTask1 clusteringTask2
NULL

#' @rdname clustering
#' @export
clusteringTask1 <- function(indices, getContribs, K1, algoClust1, nb_items_clust,
	ncores_clust=3, verbose=FALSE)
{
	if (verbose)
		cat(paste("*** Clustering task 1 on ",length(indices)," series [start]\n", sep=""))

	if (length(indices) <= K1)
		return (indices)

	parll <- (ncores_clust > 1)
	if (parll)
	{
		# outfile=="" to see stderr/stdout on terminal
		cl <-
			if (verbose)
				parallel::makeCluster(ncores_clust, outfile = "")
			else
				parallel::makeCluster(ncores_clust)
		parallel::clusterExport(cl, c("getContribs","K1","verbose"), envir=environment())
	}
	# Iterate clustering algorithm 1 until K1 medoids are found
	while (length(indices) > K1)
	{
		# Balance tasks by splitting the indices set - as evenly as possible
		indices_workers <- .splitIndices(indices, nb_items_clust, min_size=K1+1)
		indices <-
			if (parll)
			{
				unlist( parallel::parLapply(cl, indices_workers, function(inds) {
					require("epclust", quietly=TRUE)
					inds[ algoClust1(getContribs(inds), K1) ]
				}) )
			}
			else
			{
				unlist( lapply(indices_workers, function(inds)
					inds[ algoClust1(getContribs(inds), K1) ]
				) )
			}
		if (verbose)
		{
			cat(paste("*** Clustering task 1 on ",length(indices)," medoids [iter]\n", sep=""))
		}
	}
	if (parll)
		parallel::stopCluster(cl)

	indices #medoids
}

#' @rdname clustering
#' @export
clusteringTask2 <- function(indices, getSeries, K2, algoClust2, nb_series_per_chunk,
	smooth_lvl, nvoice, nbytes, endian, ncores_clust=3, verbose=FALSE)
{
	if (verbose)
		cat(paste("*** Clustering task 2 on ",length(indices)," medoids\n", sep=""))

	if (length(indices) <= K2)
		return (indices)

	# A) Compute the WER distances (Wavelets Extended coefficient of deteRmination)
	distances <- computeWerDists(indices, getSeries, nb_series_per_chunk,
		smooth_lvl, nvoice, nbytes, endian, ncores_clust, verbose)

	# B) Apply clustering algorithm 2 on the WER distances matrix
	if (verbose)
		cat(paste("*** algoClust2() on ",nrow(distances)," items\n", sep=""))
	indices[ algoClust2(distances,K2) ]
}
