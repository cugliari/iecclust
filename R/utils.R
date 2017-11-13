# Check integer arguments with functional conditions
.toInteger <- function(x, condition)
{
	errWarn <- function(ignored)
		paste("Cannot convert argument' ",substitute(x),"' to integer", sep="")
	if (!is.integer(x))
		tryCatch({x <- as.integer(x)[1]; if (is.na(x)) stop()},
			warning=errWarn, error=errWarn)
	if (!condition(x))
	{
		stop(paste("Argument '",substitute(x),
			"' does not verify condition ",body(condition), sep=""))
	}
	x
}

# Check logical arguments
.toLogical <- function(x)
{
	errWarn <- function(ignored)
		paste("Cannot convert argument' ",substitute(x),"' to logical", sep="")
	if (!is.logical(x))
		tryCatch({x <- as.logical(x)[1]; if (is.na(x)) stop()},
			warning=errWarn, error=errWarn)
	x
}

#' curvesToContribs
#'
#' Compute the discrete wavelet coefficients for each series, and aggregate them in
#' energy contribution across scales as described in https://arxiv.org/abs/1101.4744v2
#'
#' @param curves [big.]matrix of series (in columns), of size L x n
#' @param wav_filt Wavelet transform filter, as a vector c(Family,FilterNumber)
#' @inheritParams claws
#'
#' @return A matrix of size log(L) x n containing contributions in columns
#'
#' @export
curvesToContribs <- function(curves, wav_filt, contrib_type)
{
	series <- as.matrix(curves)
	L <- nrow(series)
	D <- ceiling( log2(L) )
	# Series are interpolated to all have length 2^D
	nb_sample_points <- 2^D
	apply(series, 2, function(x) {
		interpolated_curve <- spline(1:L, x, n=nb_sample_points)$y
		W <- wavethresh::wd(interpolated_curve, wav_filt[2], wav_filt[1])$D
		# Compute the sum of squared discrete wavelet coefficients, for each scale
		nrj <- sapply( 1:D, function(i) ( sqrt( sum(W[(2^D-(2^i-1)):(2^D-2^(i-1))]^2) ) ) )
		if (contrib_type!="absolute")
			nrj <- nrj / sum(nrj)
		if (contrib_type=="logit")
			nrj <- - log(1 - nrj)
		unname( nrj )
	})
}

# Helper function to divide indices into balanced sets.
# Ensure that all indices sets have at least min_size elements.
.splitIndices <- function(indices, nb_per_set, min_size=1)
{
	L <- length(indices)
	nb_workers <- floor( L / nb_per_set )
	rem <- L %% nb_per_set
	if (nb_workers == 0 || (nb_workers==1 && rem==0))
	{
		# L <= nb_per_set, simple case
		return (list(indices))
	}

	indices_workers <- lapply( seq_len(nb_workers), function(i)
		indices[(nb_per_set*(i-1)+1):(nb_per_set*i)] )

	rem <- L %% nb_per_set #number of remaining unassigned items
	if (rem == 0)
		return (indices_workers)

	rem <- (L-rem+1):L
	# If remainder is smaller than min_size, feed it with indices from other sets
	# until either its size exceed min_size (success) or other sets' size
	# get lower min_size (failure).
	while (length(rem) < min_size)
	{
		index <- length(rem) %% nb_workers + 1
		if (length(indices_workers[[index]]) <= min_size)
		{
			stop("Impossible to split indices properly for clustering.
				Try increasing nb_items_clust or decreasing K1")
		}
		rem <- c(rem, tail(indices_workers[[index]],1))
		indices_workers[[index]] <- head( indices_workers[[index]], -1)
	}
	return ( c(indices_workers, list(rem) ) )
}

#' assignMedoids
#'
#' Find the closest medoid for each curve in input
#'
#' @param curves (Chunk) of series whose medoids indices must be found
#' @param medoids Matrix of medoids (in columns)
#'
#' @return The vector of integer assignments
#' @export
assignMedoids <- function(curves, medoids)
{
	nb_series <- ncol(curves)
	mi <- rep(NA,nb_series)
	for (i in seq_len(nb_series))
		mi[i] <- which.min( colSums( sweep(medoids, 1, curves[,i], '-')^2 ) )
	mi
}

#' filterMA
#'
#' Filter [time-]series by replacing all values by the moving average of values
#' centered around current one. Border values are averaged with available data.
#'
#' @param M_ A real matrix of size LxD
#' @param w_ The (odd) number of values to average
#'
#' @return The filtered matrix (in columns), of same size as the input
#' @export
filterMA <- function(M_, w_)
	.Call("filterMA", M_, w_, PACKAGE="epclust")

#' cleanBin
#'
#' Remove binary files to re-generate them at next run of \code{claws()}.
#' To be run in the folder where computations occurred (or no effect).
#'
#' @export
cleanBin <- function()
{
	bin_files <- list.files(pattern="*.epclust.bin", all.files=TRUE)
	for (file in bin_files)
		unlink(file)
}
