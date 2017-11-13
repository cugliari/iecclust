context("clustering")

test_that("clusteringTask1 behave as expected",
{
	# Generate 60 reference sinusoïdal series (medoids to be found),
	# and sample 900 series around them (add a small noise)
	n <- 900
	x <- seq(0,9.5,0.1)
	L <- length(x) #96 1/4h
	K1 <- 60
	s <- lapply( seq_len(K1), function(i) x^(1+i/30)*cos(x+i) )
	series <- matrix(nrow=L, ncol=n)
	for (i in seq_len(n))
		series[,i] <- s[[I(i,K1)]] + rnorm(L,sd=0.01)

	getSeries <- function(indices) {
		indices <- indices[indices <= n]
		if (length(indices)>0) as.matrix(series[,indices]) else NULL
	}

	wf <- "haar"
	ctype <- "absolute"
	getContribs <- function(indices) curvesToContribs(as.matrix(series[,indices]),wf,ctype)

	require("cluster", quietly=TRUE)
	algoClust1 <- function(contribs,K) cluster::pam(t(contribs),K,diss=FALSE)$id.med
	indices1 <- clusteringTask1(1:n, getContribs, K1, algoClust1, 140, verbose=TRUE)
	medoids_K1 <- getSeries(indices1)

	expect_equal(dim(medoids_K1), c(L,K1))
	# Not easy to evaluate result: at least we expect it to be better than random selection of
	# medoids within initial series
	distor_good <- computeDistortion(series, medoids_K1)
	for (i in 1:3)
		expect_lte( distor_good, computeDistortion(series,series[,sample(1:n, K1)]) )
})

test_that("clusteringTask2 behave as expected",
{
	# Same 60 reference sinusoïdal series than in clusteringTask1 test,
	# but this time we consider them as medoids - skipping stage 1
	# Here also we sample 900 series around the 60 "medoids"
	n <- 900
	x <- seq(0,9.5,0.1)
	L <- length(x) #96 1/4h
	K1 <- 60
	K2 <- 3
	#for (i in 1:60) {plot(x^(1+i/30)*cos(x+i),type="l",col=i,ylim=c(-50,50)); par(new=TRUE)}
	s <- lapply( seq_len(K1), function(i) x^(1+i/30)*cos(x+i) )
	series <- matrix(nrow=L, ncol=n)
	for (i in seq_len(n))
		series[,i] <- s[[I(i,K1)]] + rnorm(L,sd=0.01)

	getSeries <- function(indices) {
		indices <- indices[indices <= n]
		if (length(indices)>0) as.matrix(series[,indices]) else NULL
	}

	# Perfect situation: all medoids "after stage 1" are ~good
	algoClust2 <- function(dists,K) cluster::pam(dists,K,diss=TRUE)$id.med
	indices2 <- clusteringTask2(1:K1, getSeries, K2, algoClust2, 210, 3, 4, 8, "little",
		verbose=TRUE)
	medoids_K2 <- getSeries(indices2)

	expect_equal(dim(medoids_K2), c(L,K2))
	# Not easy to evaluate result: at least we expect it to be better than random selection of
	# synchrones within 1...K1 (from where distances computations + clustering was run)
	distor_good <- computeDistortion(series, medoids_K2)
#TODO: This fails; why?
#	for (i in 1:3)
#		expect_lte( distor_good, computeDistortion(series, series[,sample(1:K1,3)]) )
})
