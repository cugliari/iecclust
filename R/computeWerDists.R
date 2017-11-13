#' computeWerDists
#'
#' Compute the WER distances between the series at specified indices, which are
#' obtaind by \code{getSeries(indices)}
#'
#' @param indices Indices of the series to consider
#' @param getSeries Function to retrieve series (argument: 'inds', integer vector),
#'   as columns of a matrix
#' @param ncores Number of cores for parallel runs
#' @inheritParams claws
#'
#' @return A distances matrix of size K x K where K == length(indices)
#'
#' @export
computeWerDists <- function(indices, getSeries, nb_series_per_chunk, smooth_lvl=3,
	nvoice=4, nbytes=4, endian=.Platform$endian, ncores=3, verbose=FALSE)
{
	n <- length(indices)
	L <- length(getSeries(1)) #TODO: not very neat way to get L
	noctave <- ceiling(log2(L)) #min power of 2 to cover serie range
	# Since a CWT contains noctave*nvoice complex series, we deduce the number of CWT to
	# retrieve/put in one chunk.
	nb_cwt_per_chunk <- max(1, floor(nb_series_per_chunk / (nvoice*noctave*2)))

	# Initialize result as a square big.matrix of size 'number of medoids'
	Xwer_dist <- bigmemory::big.matrix(nrow=n, ncol=n, type="double")

	shift <- 1 #roughly equivalent to s0 in biwavelets & cie. TODO: as arg?

	cwt_file <- tempfile(pattern="epclust_cwt.bin_")
	# Compute the getSeries(indices) CWT, and store the results in the binary file
	computeSaveCWT <- function(inds)
	{
		if (verbose)
			cat("   Compute save CWT on ",length(inds)," indices\n", sep="")

		# Obtain CWT as big vectors of real part + imaginary part (concatenate)
		ts_cwt <- sapply(inds, function(i) {
			ts <- scale(ts(getSeries(i)), center=TRUE, scale=FALSE)
			ts_cwt <- Rwave::cwt(ts, noctave+ceiling(shift/nvoice), nvoice,
				w0=2*pi, twoD=TRUE, plot=FALSE)
			ts_cwt <- ts_cwt[,(1+shift):(noctave*nvoice+shift)]
			c( as.double(Re(ts_cwt)),as.double(Im(ts_cwt)) )
		})

		# Serialization
		binarize(ts_cwt, cwt_file, nb_cwt_per_chunk, ",", nbytes, endian)
	}

	# Function to retrieve a synchrone CWT from (binary) file
	getCWT <- function(index, L)
	{
		flat_cwt <- getDataInFile(index, cwt_file, nbytes, endian)
		cwt_length <- length(flat_cwt) / 2
		re_part <- as.matrix(flat_cwt[1:cwt_length], nrow=L)
		im_part <- as.matrix(flat_cwt[(cwt_length+1):(2*cwt_length)], nrow=L)
		re_part + 1i * im_part
	}

	# Compute distances between columns i and j for j>i
	computeDistances <- function(i)
	{
		if (parll)
		{
			# parallel workers start with an empty environment
			require("epclust", quietly=TRUE)
			Xwer_dist <- bigmemory::attach.big.matrix(Xwer_dist_desc)
		}

		if (verbose)
			cat(paste("   Distances from ",i," to ",i+1,"...",n,"\n", sep=""))

		# Get CWT of column i, and run computations for columns j>i
		cwt_i <- getCWT(i, L)
		WX  <- filterMA(Mod(cwt_i * Conj(cwt_i)), smooth_lvl)

		for (j in (i+1):n)
		{
			cwt_j <- getCWT(j, L)

			# Compute the ratio of integrals formula 5.6 for WER^2
			# in https://arxiv.org/abs/1101.4744v2 paragraph 5.3
		
			
		# TODO:	
#		[X,period,scale,coix] =
#wavelet(x(:,2),dt,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother);%#ok
#[Y,period,scale,coiy] = wavelet(y(:,2),dt,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother);
#
#%Smooth X and Y before truncating!  (minimize coi)
#%sinv=1./(scale');
#%
#%
#%sX=smoothwavelet(sinv(:,ones(1,nx)).*(abs(X).^2),dt,period,Args.Dj,scale);
#%sY=smoothwavelet(sinv(:,ones(1,ny)).*(abs(Y).^2),dt,period,Args.Dj,scale);
#%
#%Wxy=X.*conj(Y);
#%
#%% ----------------------- Wavelet coherence ---------------------------------
#%sWxy=smoothwavelet(sinv(:,ones(1,n)).*Wxy,dt,period,Args.Dj,scale);
#%Rsq=abs(sWxy).^2./(sX.*sY);
	
			
			
			
			
			
			
			num <- filterMA(Mod(cwt_i * Conj(cwt_j)), smooth_lvl)
			WY <- filterMA(Mod(cwt_j * Conj(cwt_j)), smooth_lvl)
			wer2 <- sum(colSums(num)^2) / sum(colSums(WX) * colSums(WY))

			Xwer_dist[i,j] <- sqrt(L * ncol(cwt_i) * (1 - wer2))
			Xwer_dist[j,i] <- Xwer_dist[i,j]
		}
		Xwer_dist[i,i] <- 0.
	}

	if (verbose)
		cat(paste("--- Precompute and serialize synchrones CWT\n", sep=""))

	# Split indices by packets of length at most nb_cwt_per_chunk
	indices_cwt <- .splitIndices(indices, nb_cwt_per_chunk)
	# NOTE: next loop could potentially be run in //. Indices would be permuted (by
	# serialization order), and synchronicity would be required because of concurrent
	# writes. Probably not worth the effort - but possible.
	for (inds in indices_cwt)
		computeSaveCWT(inds)

	parll <- (ncores > 1)
	if (parll)
	{
		# outfile=="" to see stderr/stdout on terminal
		cl <-
			if (verbose)
				parallel::makeCluster(ncores, outfile="")
			else
				parallel::makeCluster(ncores)
		Xwer_dist_desc <- bigmemory::describe(Xwer_dist)
		parallel::clusterExport(cl, envir=environment(),
			varlist=c("parll","n","L","Xwer_dist_desc","getCWT","verbose"))
	}

	if (verbose)
		cat(paste("--- Compute WER distances\n", sep=""))

	ignored <-
		if (parll)
			parallel::clusterApplyLB(cl, seq_len(n-1), computeDistances)
		else
			lapply(seq_len(n-1), computeDistances)
	Xwer_dist[n,n] <- 0.

	if (parll)
		parallel::stopCluster(cl)

	unlink(cwt_file) #remove binary file

	Xwer_dist[,] #~small matrix K1 x K1
}
