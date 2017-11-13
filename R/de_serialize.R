#' (De)Serialization of a [big]matrix or data stream
#'
#' \code{binarize()} serializes a matrix or CSV file with minimal overhead, into a
#' binary file. \code{getDataInFile()} achieves the inverse task: she retrieves (ASCII)
#' data rows from indices in the binary file. Finally, \code{binarizeTransform()}
#' serialize transformations of all data chunks. To use it a data-retrieval function
#' must be provided -- thus \code{binarize} will most likely be used first
#' (and then a function defined to seek in generated binary file)
#'
#' @param data_ascii Matrix (by columns) or CSV file or connection (by rows)
#' @param data_bin_file Name of binary file on output of \code{binarize()}
#'   or input of \code{getDataInFile()}
#' @param nb_per_chunk Number of lines to process in one batch
#' @param getData Function to retrieve data chunks
#' @param transform Transformation function to apply on data chunks
#' @param indices Indices of the lines to retrieve
#' @inheritParams claws
#'
#' @return For \code{getDataInFile()}, a matrix with columns corresponding to the
#'   requested indices. \code{binarizeTransform()} returns the number of processed lines.
#'   \code{binarize()} is designed to serialize in several calls, thus returns nothing.
#'
#' @name de_serialize
#' @rdname de_serialize
#' @aliases binarize binarizeTransform getDataInFile
NULL

#' @rdname de_serialize
#' @export
binarize <- function(data_ascii, data_bin_file, nb_per_chunk,
	sep=",", nbytes=4, endian=.Platform$endian)
{
	# data_ascii can be of two types: [big.]matrix, or connection
	if (is.character(data_ascii))
		data_ascii <- file(data_ascii, open="r")
	else if (methods::is(data_ascii,"connection") && !isOpen(data_ascii))
		open(data_ascii)
	is_matrix <- !methods::is(data_ascii,"connection")

	# At first call, the length of a stored row is written. So it's important to determine
	# if the serialization process already started.
	first_write <- (!file.exists(data_bin_file) || file.info(data_bin_file)$size == 0)

	# Open the binary file for writing (or 'append' if already exists)
	data_bin <- file(data_bin_file, open=ifelse(first_write,"wb","ab"))

	if (first_write)
	{
		# Write data length on first call: number of items always on 8 bytes
		writeBin(0L, data_bin, size=8, endian=endian)
		if (is_matrix)
			data_length <- nrow(data_ascii)
		else #connection
		{
			# Read the first line to know data length, and write it then
			data_line <- scan(data_ascii, double(), sep=sep, nlines=1, quiet=TRUE)
			writeBin(data_line, data_bin, size=nbytes, endian=endian)
			data_length <- length(data_line)
		}
	}

	if (is_matrix)
	{
		# Data is processed by chunks; although this may not be so useful for (normal) matrix
		# input, it could for a file-backed big.matrix. It's easier to follow a unified pattern.
		index <- 1
	}
	repeat
	{
		if (is_matrix)
		{
			data_chunk <-
				if (index <= ncol(data_ascii))
					as.double(data_ascii[,index:min(ncol(data_ascii),index+nb_per_chunk-1)])
				else
					double(0)
			index <- index + nb_per_chunk
		}
		else #connection
			data_chunk <- scan(data_ascii, double(), sep=sep, nlines=nb_per_chunk, quiet=TRUE)

		# Data size is unknown in the case of a connection
		if (length(data_chunk)==0)
			break

		# Write this chunk of data to the binary file
		writeBin(data_chunk, data_bin, size=nbytes, endian=endian)
	}

	if (first_write)
	{
		# Write data_length, == (file_size-1) / (nbytes*nbWritten) at offset 0 in data_bin
		ignored <- seek(data_bin, 0)
		writeBin(data_length, data_bin, size=8, endian=endian)
	}
	close(data_bin)

	if ( ! is_matrix )
		close(data_ascii)
}

#' @rdname de_serialize
#' @export
binarizeTransform <- function(getData, transform, data_bin_file, nb_per_chunk,
	nbytes=4, endian=.Platform$endian)
{
	nb_items <- 0 #side-effect: store the number of transformed items
	index <- 1
	repeat
	{
		# Retrieve a chunk of data in a binary file (generally obtained by binarize())
		data_chunk <- getData((index-1)+seq_len(nb_per_chunk))
		if (is.null(data_chunk))
			break

		# Apply transformation on the current chunk (by columns)
		transformed_chunk <- transform(data_chunk)

		# Save the result in binary format
		binarize(transformed_chunk, data_bin_file, nb_per_chunk, ",", nbytes, endian)

		index <- index + nb_per_chunk
		nb_items <- nb_items + ncol(data_chunk)
	}
	nb_items #number of transformed items
}

#' @rdname de_serialize
#' @export
getDataInFile <- function(indices, data_bin_file, nbytes=4, endian=.Platform$endian)
{
	data_bin <- file(data_bin_file, "rb") #source binary file

	data_size <- file.info(data_bin_file)$size #number of bytes in the file
	# data_length: length of a vector in the binary file (first element, 8 bytes)
	data_length <- readBin(data_bin, "integer", n=1, size=8, endian=endian)

	# Seek all 'indices' columns in the binary file, using data_length and nbytes
	# to compute the offset ( index i at 8 + i*data_length*nbytes )
	data_ascii <- do.call( cbind, lapply( indices, function(i) {
		offset <- 8+(i-1)*data_length*nbytes
		if (offset >= data_size)
			return (NULL)
		ignored <- seek(data_bin, offset) #position cursor at computed offset
		readBin(data_bin, "double", n=data_length, size=nbytes, endian=endian)
	} ) )
	close(data_bin)

	data_ascii #retrieved data, in columns
}
