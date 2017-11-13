context("assignMedoids")

test_that("assignMedoids behave as expected",
{
	# Generate a gaussian mixture
	n <- 999
	L <- 7
	medoids <- cbind( rep(0,L), rep(-5,L), rep(5,L) )
	# short series...
	require("MASS", quietly=TRUE)
	series <- t( rbind( MASS::mvrnorm(n/3, medoids[,1], diag(L)),
		MASS::mvrnorm(n/3, medoids[,2], diag(L)),
		MASS::mvrnorm(n/3, medoids[,3], diag(L)) ) )

	# With high probability, medoids indices should resemble 1,1,1,...,2,2,2,...,3,3,3,...
	mi <- assignMedoids(series, medoids)
	mi_ref <- rep(1:3, each=n/3)
	expect_lt( mean(mi != mi_ref), 0.01 )
})
