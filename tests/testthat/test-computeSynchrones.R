context("computeSynchrones")

test_that("computeSynchrones behave as expected",
{
	# Generate 300 sinusoïdal series of 3 kinds: all series of indices == 0 mod 3 are the same
	# (plus noise), all series of indices == 1 mod 3 are the same (plus noise) ...
	n <- 300
	x <- seq(0,9.5,0.1)
	L <- length(x) #96 1/4h
	K <- 3
	s1 <- cos(x)
	s2 <- sin(x)
	s3 <- c( s1[1:(L%/%2)] , s2[(L%/%2+1):L] )
	#sum((s1-s2)^2) == 96
	#sum((s1-s3)^2) == 58
	#sum((s2-s3)^2) == 38
	s <- list(s1, s2, s3)
	series <- matrix(nrow=L, ncol=n)
	for (i in seq_len(n))
		series[,i] <- s[[I(i,K)]] + rnorm(L,sd=0.01)

	getSeries <- function(indices) {
		indices <- indices[indices <= n]
		if (length(indices)>0) as.matrix(series[,indices]) else NULL
	}

	synchrones <- computeSynchrones(cbind(s1,s2,s3),getSeries,n,100,verbose=TRUE)

	expect_equal(dim(synchrones), c(L,K))
	for (i in 1:K)
	{
		# Synchrones are (for each medoid) sums of closest curves.
		# Here, we expect exactly 100 curves of each kind to be assigned respectively to
		# synchrone 1, 2 and 3 => division by 100 should be very close to the ref curve
		expect_equal(synchrones[,i]/100, s[[i]], tolerance=0.01)
	}
})
