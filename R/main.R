#' CLAWS: CLustering with wAvelets and Wer distanceS
#'
#' Cluster electricity power curves (or any series of similar nature) by applying a
#' two stage procedure in parallel (see details).
#' Input series must be sampled on the same time grid, no missing values.
#'
#' Summary of the function execution flow:
#' \enumerate{
#'   \item Compute and serialize all contributions, obtained through discrete wavelet
#'     decomposition (see Antoniadis & al. [2013])
#'   \item Divide series into \code{ntasks} groups to process in parallel. In each task:
#'   \enumerate{
#'     \item iterate the first clustering algorithm on its aggregated outputs,
#'       on inputs of size \code{nb_items_clust}\cr
#'         -> K1 medoids indices
#'     \item optionally, if WER=="mix":\cr
#'       a. compute WER distances (K1xK1) between medoids\cr
#'       b. apply the 2nd clustering algorithm\cr
#'          -> K2 medoids indices
#'   }
#'   \item Launch a final task on the aggregated outputs of all previous tasks:
#'     ntasks*K1 if WER=="end", ntasks*K2 otherwise
#'   \item Compute synchrones (sum of series within each final group)
#' }
#' 
#' The main argument -- \code{series} -- has a quite misleading name, since it can be
#' either a [big.]matrix, a CSV file, a connection or a user function to retrieve series.
#' When \code{series} is given as a function it must take a single argument,
#' 'indices': integer vector equal to the indices of the curves to retrieve;
#' see SQLite example.
#' WARNING: the return value must be a matrix (in columns), or NULL if no matches.
#' 
#' Note: Since we don't make assumptions on initial data, there is a possibility that
#' even when serialized, contributions do not fit in RAM. For example,
#' 30e6 series of length 100,000 would lead to a +4Go contribution matrix. Therefore,
#' it's safer to place these in (binary) files; that's what we do.
#'
#' @param series Access to the N (time-)series, which can be of one of the four
#'   following types:
#'   \itemize{
#'     \item [big.]matrix: each column contains the (time-ordered) values of one
#'       time-serie
#'     \item connection: any R connection object providing lines as described above
#'     \item character: name of a CSV file containing series in rows (no header)
#'     \item function: a custom way to retrieve the curves; it has only one argument:
#'       the indices of the series to be retrieved. See SQLite example
#'   }
#' @param K1 Number of clusters to be found after stage 1 (K1 << N)
#' @param K2 Number of clusters to be found after stage 2 (K2 << K1)
#' @param nb_series_per_chunk Number of series to retrieve in one batch
#' @param nb_items_clust Number of items in 1st clustering algorithm input
#' @param algoClust1 Clustering algorithm for stage 1. A function which takes (data, K)
#'   as argument where data is a matrix in columns and K the desired number of clusters,
#'   and outputs K medoids ranks. Default: PAM.
#' @param algoClust2 Clustering algorithm for stage 2. A function which takes (dists, K)
#'   as argument where dists is a matrix of distances and K the desired number of
#'   clusters, and outputs K medoids ranks. Default: PAM.
#' @param wav_filt Wavelet transform filter, as a string "Family:FilterNumber"; see
#'   ?wavethresh::wd
#' @param contrib_type Type of contribution: "relative", "logit" or "absolute" (any
#'   prefix)
#' @param WER "end" to apply stage 2 after stage 1 has fully iterated, or "mix" to apply
#'   stage 2 at the end of each task
#' @param smooth_lvl Smoothing level: odd integer, 1 == no smoothing.
#' @param nvoice Number of voices within each octave for CWT computations
#' @param random TRUE (default) for random chunks repartition
#' @param ntasks Number of tasks (parallel iterations to obtain K1 [if WER=="end"]
#'   or K2 [if WER=="mix"] medoids); default: 1.\cr
#'   Note: ntasks << N (number of series), so that N is "roughly divisible" by ntasks
#' @param ncores_tasks Number of parallel tasks ('1' == sequential tasks)
#' @param ncores_clust Number of parallel clusterings in one task
#' @param sep Separator in CSV input file (if any provided)
#' @param nbytes 4 or 8 bytes to (de)serialize a floating-point number
#' @param endian Endianness for (de)serialization: "little" or "big"
#' @param verbose FALSE: nothing printed; TRUE: some execution traces
#'
#' @return A list:
#' \itemize{
#'   \item medoids: matrix of the final K2 medoids curves
#'   \item ranks: corresponding indices in the dataset
#'   \item synchrones: sum of series within each final group
#' }
#'
#' @references Clustering functional data using Wavelets [2013];
#'   A. Antoniadis, X. Brossat, J. Cugliari & J.-M. Poggi.
#'   Inter. J. of Wavelets, Multiresolution and Information Procesing,
#'   vol. 11, No 1, pp.1-30. doi:10.1142/S0219691313500033
#'
#' @examples
#' \dontrun{
#' # WER distances computations are too long for CRAN (for now)
#' # Note: on this small example, sequential run is faster
#'
#' # Random series around cos(x,2x,3x)/sin(x,2x,3x)
#' x <- seq(0,50,0.05)
#' L <- length(x) #1001
#' ref_series <- matrix( c(cos(x),cos(2*x),cos(3*x),sin(x),sin(2*x),sin(3*x)), ncol=6 )
#' library(wmtsa)
#' series <- do.call( cbind, lapply( 1:6, function(i)
#'   do.call(cbind, wmtsa::wavBootstrap(ref_series[,i], n.realization=40)) ) )
#' # Mix series so that all groups are evenly spread
#' permut <- (0:239)%%6 * 40 + (0:239)%/%6 + 1
#' series = series[,permut]
#' #dim(series) #c(240,1001)
#' res_ascii <- claws(series, K1=30, K2=6, nb_series_per_chunk=500,
#'   nb_items_clust=100, random=FALSE, verbose=TRUE, ncores_clust=1)
#'
#' # Same example, from CSV file
#' csv_file <- tempfile(pattern="epclust_series.csv_")
#' write.table(t(series), csv_file, sep=",", row.names=FALSE, col.names=FALSE)
#' res_csv <- claws(csv_file, 30, 6, 500, 100, random=FALSE, ncores_clust=1)
#'
#' # Same example, from binary file
#' bin_file <- tempfile(pattern="epclust_series.bin_")
#' nbytes <- 8
#' endian <- "little"
#' binarize(csv_file, bin_file, 500, ",", nbytes, endian)
#' getSeries <- function(indices) getDataInFile(indices, bin_file, nbytes, endian)
#' res_bin <- claws(getSeries, 30, 6, 500, 100, random=FALSE, ncores_clust=1)
#' unlink(csv_file)
#' unlink(bin_file)
#'
#' # Same example, from SQLite database
#' library(DBI)
#' series_db <- dbConnect(RSQLite::SQLite(), "file::memory:")
#' # Prepare data.frame in DB-format
#' n <- ncol(series)
#' times_values <- data.frame(
#'   id = rep(1:n,each=L),
#'   time = rep( as.POSIXct(1800*(1:L),"GMT",origin="2001-01-01"), n ),
#'   value = as.double(series) )
#' dbWriteTable(series_db, "times_values", times_values)
#' # Fill associative array, map index to identifier
#' indexToID_inDB <- as.character(
#'   dbGetQuery(series_db, 'SELECT DISTINCT id FROM times_values')[,"id"] )
#' serie_length <- as.integer( dbGetQuery(series_db,
#'   paste("SELECT COUNT(*) FROM times_values WHERE id == ",indexToID_inDB[1],sep="")) )
#' getSeries <- function(indices) {
#'   indices = indices[ indices <= length(indexToID_inDB) ]
#'   if (length(indices) == 0)
#'     return (NULL)
#'   request <- "SELECT id,value FROM times_values WHERE id in ("
#'   for (i in seq_along(indices)) {
#'     request <- paste(request, indexToID_inDB[ indices[i] ],  sep="")
#'     if (i < length(indices))
#'       request <- paste(request, ",", sep="")
#'   }
#'   request <- paste(request, ")", sep="")
#'   df_series <- dbGetQuery(series_db, request)
#'   matrix(df_series[,"value"], nrow=serie_length)
#' }
#' res_db <- claws(getSeries, 30, 6, 500, 100, random=FALSE, ncores_clust=1)
#' dbDisconnect(series_db)
#'
#' # All results should be equal:
#' all(res_ascii$ranks == res_csv$ranks
#'   & res_ascii$ranks == res_bin$ranks
#'   & res_ascii$ranks == res_db$ranks)
#' }
#' @export
claws <- function(series, K1, K2, nb_series_per_chunk, nb_items_clust=5*K1,
	algoClust1=function(data,K) cluster::pam(t(data),K,diss=FALSE,pamonce=1)$id.med,
	algoClust2=function(dists,K) cluster::pam(dists,K,diss=TRUE,pamonce=1)$id.med,
	wav_filt="Coiflets:1", contrib_type="absolute", WER="end", smooth_lvl=3, nvoice=4,
	random=TRUE, ntasks=1, ncores_tasks=1, ncores_clust=3, sep=",", nbytes=4,
	endian=.Platform$endian, verbose=FALSE)
{
	# Check/transform arguments
	if (!is.matrix(series) && !bigmemory::is.big.matrix(series)
		&& !is.function(series)
		&& !methods::is(series,"connection") && !is.character(series))
	{
		stop("'series': [big]matrix, function, file or valid connection (no NA)")
	}
	K1 <- .toInteger(K1, function(x) x>=2)
	K2 <- .toInteger(K2, function(x) x>=2)
	nb_series_per_chunk <- .toInteger(nb_series_per_chunk, function(x) x>=1)
	nb_items_clust <- .toInteger(nb_items_clust, function(x) x>K1)
	random <- .toLogical(random)
	wav_filt <- strsplit(wav_filt, ':')[[1]]
	wav_family <- wav_filt[1]
	wav_number <- wav_filt[2]
	tryCatch({ignored <- wavethresh::filter.select(wav_number, wav_family)},
		error=function(e) stop("Invalid wavelet filter; see ?wavethresh::filter.select") )
	ctypes <- c("relative","absolute","logit")
	contrib_type <- ctypes[ pmatch(contrib_type,ctypes) ]
	if (is.na(contrib_type))
		stop("'contrib_type' in {'relative','absolute','logit'}")
	if (WER!="end" && WER!="mix")
		stop("'WER': in {'end','mix'}")
	random <- .toLogical(random)
	ntasks <- .toInteger(ntasks, function(x) x>=1)
	ncores_tasks <- .toInteger(ncores_tasks, function(x) x>=1)
	ncores_clust <- .toInteger(ncores_clust, function(x) x>=1)
	if (!is.character(sep))
		stop("'sep': character")
	nbytes <- .toInteger(nbytes, function(x) x==4 || x==8)
	verbose <- .toLogical(verbose)

	# Binarize series if it is not a function; the aim is to always use a function,
	# to uniformize treatments. An equally good alternative would be to use a file-backed
	# bigmemory::big.matrix, but it would break the "all-is-function" pattern.
	if (!is.function(series))
	{
		if (verbose)
			cat("...Serialize time-series (or retrieve past binary file)\n")
		series_file <- ".series.epclust.bin"
		if (!file.exists(series_file))
			binarize(series, series_file, nb_series_per_chunk, sep, nbytes, endian)
		getSeries <- function(inds) getDataInFile(inds, series_file, nbytes, endian)
	}
	else
		getSeries <- series

	# Serialize all computed wavelets contributions into a file
	contribs_file <- ".contribs.epclust.bin"
	if (verbose)
		cat("...Compute contributions and serialize them (or retrieve past binary file)\n")
	if (!file.exists(contribs_file))
	{
		nb_curves <- binarizeTransform(getSeries,
			function(curves) curvesToContribs(curves, wav_filt, contrib_type),
			contribs_file, nb_series_per_chunk, nbytes, endian)
	}
	else
	{
		# TODO: duplicate from getDataInFile() in de_serialize.R
		contribs_size <- file.info(contribs_file)$size #number of bytes in the file
		contrib_length <- readBin(contribs_file, "integer", n=1, size=8, endian=endian)
		nb_curves <- (contribs_size-8) / (nbytes*contrib_length)
	}
	getContribs <- function(indices) getDataInFile(indices, contribs_file, nbytes, endian)

	# A few sanity checks: do not continue if too few data available.
	if (nb_curves < K2)
		stop("Not enough data: less series than final number of clusters")
	nb_series_per_task <- round(nb_curves / ntasks)
	if (nb_series_per_task < K2)
		stop("Too many tasks: less series in one task than final number of clusters")

	# Generate a random permutation of 1:N (if random==TRUE);
	# otherwise just use arrival (storage) order.
	indices_all <- if (random) sample(nb_curves) else seq_len(nb_curves)
	# Split (all) indices into ntasks groups of ~same size
	indices_tasks <- lapply(seq_len(ntasks), function(i) {
		upper_bound <- ifelse( i<ntasks, min(nb_series_per_task*i,nb_curves), nb_curves )
		indices_all[((i-1)*nb_series_per_task+1):upper_bound]
	})

	parll <- (ncores_tasks > 1)
	if (parll && ntasks>1)
	{
		# Initialize parallel runs: outfile="" allow to output verbose traces in the console
		# under Linux. All necessary variables are passed to the workers.
		cl <-
			if (verbose)
				parallel::makeCluster(ncores_tasks, outfile="")
			else
				parallel::makeCluster(ncores_tasks)
		varlist <- c("ncores_clust","verbose", #task 1 & 2
			"K1","getContribs","algoClust1","nb_items_clust") #task 1
		if (WER=="mix")
		{
			# Add variables for task 2
			varlist <- c(varlist, "K2","getSeries","algoClust2","nb_series_per_chunk",
				"smooth_lvl","nvoice","nbytes","endian")
		}
		parallel::clusterExport(cl, varlist, envir <- environment())
	}

	# This function achieves one complete clustering task, divided in stage 1 + stage 2.
	# stage 1: n indices  --> clusteringTask1(...) --> K1 medoids (indices)
	# stage 2: K1 indices --> K1xK1 WER distances --> clusteringTask2(...) --> K2 medoids,
	# where n == N / ntasks, N being the total number of curves.
	runTwoStepClustering <- function(inds)
	{
		# When running in parallel, the environment is blank: we need to load the required
		# packages, and pass useful variables.
		if (parll && ntasks>1)
			require("epclust", quietly=TRUE)
		indices_medoids <- clusteringTask1(inds, getContribs, K1, algoClust1,
			nb_items_clust, ncores_clust, verbose)
		if (WER=="mix")
		{
			indices_medoids <- clusteringTask2(indices_medoids, getSeries, K2, algoClust2,
				nb_series_per_chunk,smooth_lvl,nvoice,nbytes,endian,ncores_clust,verbose)
		}
		indices_medoids
	}

	if (verbose)
	{
		message <- paste("...Run ",ntasks," x stage 1", sep="")
		if (WER=="mix")
			message <- paste(message," + stage 2", sep="")
		cat(paste(message,"\n", sep=""))
	}

	# As explained above, we obtain after all runs ntasks*[K1 or K2] medoids indices,
	# depending whether WER=="end" or "mix", respectively.
	indices_medoids_all <-
		if (parll && ntasks>1)
			unlist( parallel::parLapply(cl, indices_tasks, runTwoStepClustering) )
		else
			unlist( lapply(indices_tasks, runTwoStepClustering) )

	if (parll && ntasks>1)
		parallel::stopCluster(cl)

	# For the last stage, ncores_tasks*(ncores_clusts+1) cores should be available:
	#  - ntasks for level 1 parallelism
	#  - ntasks*ncores_clust for level 2 parallelism,
	# but since an extension MPI <--> tasks / OpenMP <--> sub-tasks is on the way,
	# it's better to just re-use ncores_clust
	ncores_last_stage <- ncores_clust

	# Run last clustering tasks to obtain only K2 medoids indices
	if (verbose)
		cat("...Run final // stage 1 + stage 2\n")
	indices_medoids <- clusteringTask1(indices_medoids_all, getContribs, K1, algoClust1,
		nb_items_clust, ncores_tasks*ncores_clust, verbose)

	indices_medoids <- clusteringTask2(indices_medoids, getSeries, K2, algoClust2,
		nb_series_per_chunk,smooth_lvl,nvoice,nbytes,endian,ncores_last_stage,verbose)

	# Compute synchrones, that is to say the cumulated power consumptions for each of the
	# K2 final groups.
	medoids <- getSeries(indices_medoids)
	synchrones <- computeSynchrones(medoids, getSeries, nb_curves, nb_series_per_chunk,
		ncores_last_stage, verbose)

	# NOTE: no need to use big.matrix here, because only K2 << K1 << N remaining curves
	list("medoids"=medoids, "ranks"=indices_medoids, "synchrones"=synchrones)
}
