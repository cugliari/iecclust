context("de_serialize")

test_that("serialization + getDataInFile retrieve original data / from matrix",
{
	data_bin_file <- ".epclust_test_m.bin"
	unlink(data_bin_file)

	# Dataset 200 cols / 30 rows
	data_ascii <- matrix(runif(200*30,-10,10),nrow=30)
	nbytes <- 4 #lead to a precision of 1e-7 / 1e-8
	endian <- "little"

	# Simulate serialization in one single call
	binarize(data_ascii, data_bin_file, 500, ",", nbytes, endian)
	expect_equal(file.info(data_bin_file)$size, length(data_ascii)*nbytes+8)
	for (indices in list(c(1,3,5), 3:13, c(5,20,50), c(75,130:135), 196:200))
	{
		data_lines <- getDataInFile(indices, data_bin_file, nbytes, endian)
		expect_equal(data_lines, data_ascii[,indices], tolerance=1e-6)
	}
	unlink(data_bin_file)

	# Serialization in several calls (last call complete, next call NULL)
	for (i in 1:20)
		binarize(data_ascii[,((i-1)*10+1):(i*10)], data_bin_file, 20, ",", nbytes, endian)
	expect_equal(file.info(data_bin_file)$size, length(data_ascii)*nbytes+8)
	for (indices in list(c(1,3,5), 3:13, c(5,20,50), c(75,130:135), 196:200))
	{
		data_lines <- getDataInFile(indices, data_bin_file, nbytes, endian)
		expect_equal(data_lines, data_ascii[,indices], tolerance=1e-6)
	}
	unlink(data_bin_file)
})

test_that("serialization + transform + getDataInFile retrieve original transforms",
{
	data_bin_file <- ".epclust_test_t.bin"
	unlink(data_bin_file)

	# Dataset 200 cols / 30 rows
	data_ascii <- matrix(runif(200*30,-10,10),nrow=30)
	nbytes <- 8
	endian <- "little"

	binarize(data_ascii, data_bin_file, 500, ",", nbytes, endian)
	# Serialize transformation (just compute range) into a new binary file
	trans_bin_file <- ".epclust_test_t_trans.bin"
	unlink(trans_bin_file)
	getSeries <- function(inds) getDataInFile(inds, data_bin_file, nbytes, endian)
	binarizeTransform(getSeries, function(series) apply(series, 2, range),
		trans_bin_file, 250, nbytes, endian)
	unlink(data_bin_file)
	expect_equal(file.info(trans_bin_file)$size, 2*ncol(data_ascii)*nbytes+8)
	for (indices in list(c(1,3,5), 3:13, c(5,20,50), c(75,130:135), 196:200))
	{
		trans_cols <- getDataInFile(indices, trans_bin_file, nbytes, endian)
		expect_equal(trans_cols, apply(data_ascii[,indices],2,range), tolerance=1e-6)
	}
	unlink(trans_bin_file)
})

test_that("serialization + getDataInFile retrieve original data / from connection",
{
	data_bin_file <- ".epclust_test_c.bin"
	unlink(data_bin_file)

	# Dataset 300 cols / 50 rows
	data_csv <- system.file("testdata","de_serialize.csv",package="epclust")
	nbytes <- 8
	endian <- "big"
	data_ascii <- unname( t( as.matrix(read.table(data_csv,sep=";",header=FALSE)) ) ) #ref

	# Simulate serialization in one single call
	binarize(data_csv, data_bin_file, 350, ";", nbytes, endian)
	expect_equal(file.info(data_bin_file)$size, 300*50*8+8)
	for (indices in list(c(1,3,5), 3:13, c(5,20,50), c(75,130:135), 196:200))
	{
		data_cols <- getDataInFile(indices,data_bin_file,nbytes,endian)
		expect_equal(data_cols, data_ascii[,indices])
	}
	unlink(data_bin_file)

	# Serialization in several calls / chunks of 29 --> 29*10 + 10, incomplete last
	data_con <- file(data_csv, "r")
	binarize(data_con, data_bin_file, 29, ";", nbytes, endian)
	expect_equal(file.info(data_bin_file)$size, 300*50*8+8)
	for (indices in list(c(1,3,5), 3:13, c(5,20,50), c(75,130:135), 196:200))
	{
		data_cols <- getDataInFile(indices,data_bin_file,nbytes,endian)
		expect_equal(data_cols, data_ascii[,indices])
	}
	unlink(data_bin_file)
	#close(data_con) --> done in binarize()
})
