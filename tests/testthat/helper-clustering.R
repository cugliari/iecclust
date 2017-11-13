# Compute the sum of (normalized) sum of squares of closest distances to a medoid.
computeDistortion <- function(series, medoids)
{
	n <- ncol(series)
	L <- nrow(series)
	distortion <- 0.
	for (i in seq_len(n))
		distortion <- distortion + min( colSums( sweep(medoids,1,series[,i],'-')^2 ) / L )

	sqrt( distortion / n )
}
