context("computeWerDists")

test_that("computeWerDists output correct results",
{
	nbytes <- 8
	endian <- "little"

	# On two identical series
	serie <- rnorm(212, sd=5)
	series <- cbind(serie, serie)
	getSeries <- function(indices) as.matrix(series[,indices])
	dists <- computeWerDists(1:2, getSeries, 50, 3, 4, nbytes, endian,
		verbose=TRUE)
	expect_equal(dists, matrix(0.,nrow=2,ncol=2))

	# On two constant series
	# TODO: ref results. Ask Jairo to check function.
})
