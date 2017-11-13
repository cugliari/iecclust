context("utils functions")

test_that("Helper function to split indices work properly",
{
	indices <- 1:400

	# bigger nb_per_set than length(indices)
	expect_equal(epclust:::.splitIndices(indices,500), list(indices))

	# nb_per_set == length(indices)
	expect_equal(epclust:::.splitIndices(indices,400), list(indices))

	# length(indices) %% nb_per_set == 0
	expect_equal(epclust:::.splitIndices(indices,200),
		c( list(indices[1:200]), list(indices[201:400]) ))
	expect_equal(epclust:::.splitIndices(indices,100),
		c( list(indices[1:100]), list(indices[101:200]),
			list(indices[201:300]), list(indices[301:400]) ))

	# length(indices) / nb_per_set == 1, length(indices) %% nb_per_set == 100
	expect_equal(epclust:::.splitIndices(indices,300,min_size=1),
		list(1:300, 301:400))
	split_inds <- epclust:::.splitIndices(indices,300,min_size=200)
	expect_equal(length(unique(unlist(split_inds))), 400)
	expect_equal(length(split_inds), 2)
	expect_equal(length(split_inds[[1]]), 200)
	expect_equal(length(split_inds[[2]]), 200)
	expect_error(epclust:::.splitIndices(indices,300,min_size=300), "Impossible to split*")

	# length(indices) / nb_per_set == 2, length(indices) %% nb_per_set == 42
	expect_equal(epclust:::.splitIndices(indices,179,min_size=1),
		list(1:179, 180:358, 359:400))
	split_inds <-epclust:::.splitIndices(indices,179,min_size=60)
	expect_equal(length(unique(unlist(split_inds))), 400)
	expect_equal(length(split_inds), 3)
	expect_equal(length(split_inds[[1]]), 170)
	expect_equal(length(split_inds[[2]]), 170)
	expect_equal(length(split_inds[[3]]), 60)
	expect_error(epclust:::.splitIndices(indices,179,min_size=150), "Impossible to split*")
})

test_that("curvesToContribs output correct results",
{
	L <- 220 #extended to 256, log2(256) == 8

	# Zero serie
	expect_equal(curvesToContribs(rep(0,L), "d8", "absolute"), as.matrix(rep(0,8)))
	expect_equal(curvesToContribs(rep(0,L), "haar", "absolute"), as.matrix(rep(0,8)))

	# Constant serie
	expect_equal(curvesToContribs(rep(5,L), "haar", "absolute"), as.matrix(rep(0,8)))
	expect_equal(curvesToContribs(rep(10,L), "haar", "absolute"), as.matrix(rep(0,8)))
#	expect_equal(curvesToContribs(rep(5,L), "d8", ctype), rep(0,8))
#TODO:
})
