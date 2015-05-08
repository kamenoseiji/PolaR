#-------- Extract number of channel
GetChNum <- function(fname){
	file_ptr <- file(fname, "rb")								# Access to the file
	header <- readBin(file_ptr, what=integer(), size=4, n=32)	# Read header information
	close(file_ptr)
	return(header[16])
}

#-------- Read PolariS Power Spectrum (A file)
readPolariS <- function(fname){
	chnum <- GetChNum(fname)
	head_size <- 128						# 128-byte Header
	file_size <- file.info(fname)$size		# File size [byte]
	file_ptr <- file(fname, "rb")			# Access to the file
	header <- readBin(file_ptr, what=integer(), size=4, n=32)	# Read header information
	ipnum <- (file_size - head_size) / (4* chnum)	# 1 record = sizeof(float) x chnum
	temparray <- array(readBin(file_ptr, what=numeric(), n=chnum* ipnum, size=4), dim=c(chnum, ipnum))
	close(file_ptr)
	return(temparray)
}

#-------- Read PolariS Bit Distribution (P file)
readBitDist <- function(fname){
	head_size <- 128						# 128-byte Header
	file_size <- file.info(fname)$size		# File size [byte]
	file_ptr <- file(fname, "rb")			# Access to the file
	header   <- readBin(file_ptr, what=integer(), size=4, n=32)	# Read header information
	levelnum <- 2^header[11]
	ipnum <- (file_size - head_size) / (4* levelnum)	# 1 record = sizeof(int) x levelnum
	temparray <- array(readBin(file_ptr, what=integer(), n=levelnum* ipnum, size=4), dim=c(levelnum, ipnum))
	close(file_ptr)
	return(temparray)
}

#-------- Read PolariS Cross Power Spectrum (C file)
readPolariS_X <- function(fname){
	head_size <- 128						# 128-byte Header
	file_size <- file.info(fname)$size		# File size [byte]
	file_ptr <- file(fname, "rb")			# Access to the file
	header   <- readBin(file_ptr, what=integer(), size=4, n=32)	# Read header information
	chnum <- header[16]						# Number of channels in the header
	ipnum <- (file_size - head_size) / (8* chnum)	# 1 record = sizeof(float) x chnum x (re, im)
	tempvec <- readBin(file_ptr, what=numeric(), n=2*chnum* ipnum, size=4)
	close(file_ptr)
	even_index <- (0:(chnum* ipnum-1))*2; odd_index <- even_index + 1; even_index <- odd_index + 1
	tempcomplex <- complex(real = tempvec[odd_index], imaginary=tempvec[even_index])  # real and imaginary part
	return(array(tempcomplex, dim=c(chnum, ipnum)) )
}

#-------- Read PolariS Power Spectrum (A file) and integrate the spectrum referring given scan pattern
readPolariS_Scan <- function(fname, scan){
	head_size <- 128						# 128-byte Header
	file_size <- file.info(fname)$size		# File size [byte]
	file_ptr <- file(fname, "rb")			# Access to the file
	header   <- readBin(file_ptr, what=integer(), size=4, n=32)	# Read header information
	chnum <- header[16]						# Number of channels in the header
	ipnum <- (file_size - head_size) / (4* chnum)	# 1 record = sizeof(float) x chnum
	pattern <- rep(0, ipnum); pattern[scan] <- 1
	spec <- numeric(chnum)
	for(recIndex in 1:ipnum){
		spec <- spec + pattern[recIndex]* readBin(file_ptr, what=numeric(), n=chnum, size=4)
	}
	close(file_ptr)
	return(spec / sum(pattern))
}

#-------- Read PolariS Cross Power Spectrum (C file) and integrate the spectrum referring given scan pattern
readPolariS_XScan <- function(fname, scan){
	head_size <- 128						# 128-byte Header
	file_size <- file.info(fname)$size		# File size [byte]
	file_ptr <- file(fname, "rb")			# Access to the file
	header   <- readBin(file_ptr, what=integer(), size=4, n=32)	# Read header information
	chnum <- header[16]						# Number of channels in the header
	ipnum <- (file_size - head_size) / (8* chnum)	# 1 record = sizeof(float) x chnum x (re, im)
	pattern <- rep(0, ipnum); pattern[scan] <- 1
	even_index <- (0:(chnum-1))*2; odd_index <- even_index + 1; even_index <- odd_index + 1
	spec <- complex(chnum)
	for(recIndex in 1:ipnum){
		tempvec <- readBin(file_ptr, what=numeric(), n=2*chnum, size=4)
		if(pattern[recIndex] == 1){
			spec <- spec + complex(real = tempvec[odd_index], imaginary=tempvec[even_index])
		}
	}
	close(file_ptr)
	return(spec / sum(pattern))
}


#-------- Function to pick prefix that covers specified MJD
findPrefix <- function(mjdSec, prefix){
	mjdPrefix <- as.numeric(lapply(prefix, FUN=prefix2MJDsec))
	index <- which(mjdPrefix <= mjdSec)
	if( length(index) == 0){	return(-1)}
	return(max(index))
}



