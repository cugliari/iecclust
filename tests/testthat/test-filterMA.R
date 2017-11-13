context("filterMA")

test_that("[time-]serie filtering behave as expected",
{
	# Mean of 3 values
	M <- matrix(runif(1000,min=-7,max=7), ncol=10)
	ref_fM <- stats::filter(M, rep(1/3,3), circular=FALSE)
	fM <- epclust:::filterMA(M, 3)

	# Expect an agreement on all inner values
	expect_equal(dim(fM), c(100,10))
	expect_equal(fM[2:99,], ref_fM[2:99,])

	# Border values should be averages of 2 values
	expect_equal(fM[1,], colMeans(M[1:2,]))
	expect_equal(fM[100,], colMeans(M[99:100,]))

	# Mean of 5 values
	ref_fM <- stats::filter(M, rep(1/5,5), circular=FALSE)
	fM <- epclust:::filterMA(M, 5)

	# Expect an agreement on all inner values
	expect_equal(dim(fM), c(100,10))
	expect_equal(fM[3:98,], ref_fM[3:98,])

	# Border values should be averages of 3 or 4 values
	expect_equal(fM[1,], colMeans(M[1:3,]))
	expect_equal(fM[2,], colMeans(M[1:4,]))
	expect_equal(fM[99,], colMeans(M[97:100,]))
	expect_equal(fM[100,], colMeans(M[98:100,]))
})
