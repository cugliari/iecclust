#' computeSynchrones
#'
#' Compute the synchrones curves (sums of clusters elements) from a matrix of medoids,
#' using euclidian distance.
#'
#' @param medoids matrix of K medoids curves in columns
#' @param nb_curves How many series? (this is known, at this stage)
#' @inheritParams claws
#' @inheritParams computeWerDists
#'
#' @return A matrix of K synchrones in columns (same length as the series)
#'
#' @export
computeSynchrones <- function(medoids, getSeries, nb_curves,
	nb_series_per_chunk, ncores=3, verbose=FALSE)
{
	# Synchrones computation is embarassingly parallel: compute it by chunks of series
	computeSynchronesChunk <- function(indices)
	{
		# Obtain a chunk of reference series
		series_chunk <- getSeries(indices)
		nb_series_chunk <- ncol(series_chunk)

		# Get medoids indices for this chunk of series
		mi <- assignMedoids(series_chunk, medoids[,])

		# Update synchrones using mi above, grouping it by values of mi (in 1...K)
		# to avoid too many lock/unlock
		for (i in seq_len(K))
		{
			# lock / unlock required because several writes at the same time
			if (parll)
				synchronicity::lock(m)
			synchrones[,i] <- synchrones[,i] + rowSums(as.matrix(series_chunk[,mi==i]))
			if (parll)
				synchronicity::unlock(m)
		}
		NULL
	}

	K <- ncol(medoids)
	L <- nrow(medoids)
	# Use bigmemory (shared==TRUE by default) + synchronicity to fill synchrones in //
	synchrones <- bigmemory::big.matrix(nrow=L, ncol=K, type="double", init=0.)
	# NOTE: synchronicity is only for Linux & MacOS; on Windows: run sequentially
	parll <- (ncores > 1 && requireNamespace("synchronicity",quietly=TRUE)
		&& Sys.info()['sysname'] != "Windows")

	if (parll)
		m <- synchronicity::boost.mutex() #for lock/unlock, see computeSynchronesChunk

	if (verbose)
		cat(paste("--- Compute ",K," synchrones with ",nb_curves," series\n", sep=""))

	# Balance tasks by splitting 1:nb_curves into groups of size <= nb_series_per_chunk
	indices_workers <- .splitIndices(seq_len(nb_curves), nb_series_per_chunk)
	ignored <-
		if (parll)
		{
			parallel::mclapply(indices_workers,
				function(inds) computeSynchronesChunk(inds), mc.cores=ncores)
		}
		else
			lapply(indices_workers, computeSynchronesChunk)

	return (synchrones[,])
}
