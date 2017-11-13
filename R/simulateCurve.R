Ref_Centroids <<- NULL
Tot_Nb_Curves <<- NA

#' initCentroids
#'
#' Must be called before getDetRandCurve()
#'
#' @param N Number of curves
#' @param D Number of sample points
#' @param K Number of clusters
#'
#' @export
initClustersParams <- function(N, D, K)
{
	# Generate K centroids
	Ref_Centroids <<- sapply(1:K, function(k) cumsum(rnorm(D)))
	Tot_Nb_Curves <<- N
}

#' getDetRandCurve
#'
#' TODO: improve this (using known clusters centers...)
#'
#' @param i Some index
#'
#' @examples
#' \dontrun{
#' initCentroids(N=250000, D=15000, K=15)
#' res_ascii <- claws(getDetRandCurve, K1=50, K2=15, nb_series_per_chunk=500,
#'   nb_items_clust=100, random=FALSE, verbose=TRUE, ncores_clust=3)
#' }
#' @export
getDetRandCurve <- function(indices)
{
	sapply(indices, function(i) {
		if (i > Tot_Nb_Curves)
			return (NULL)
		set.seed(i)
		j <- sample(ncol(Ref_Centroids), 1)
		Ref_Centroids[,j] + rnorm(nrow(Ref_Centroids))
	})
}
