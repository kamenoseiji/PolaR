library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/delaysearch.R", ssl.verifypeer = FALSE)))

#-------- Read PolariS Power Spectrum (A file)
readPolariS <- function(fname){
	head_size <- 128						# 128-byte Header
	file_size <- file.info(fname)$size		# File size [byte]
	file_ptr <- file(fname, "rb")			# Access to the file
	header <- readBin(file_ptr, what=integer(), size=4, n=32)	# Read header information
	chnum <- header[16]						# Number of channels in the header
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

#-------- Bunching given vector
bunch_vec <- function(vector, bunchNum){	# Function to bunch 1-d vector
	pudding <- (bunchNum - length(vector) %% bunchNum) %% bunchNum
	temp <- matrix( append(vector, rep(vector[length(vector)], pudding)), nrow=bunchNum)
	return( apply(temp, 2, mean))
}

wireGridPhaseCorr <- function(XY, ScanIndex, chRange=9:65544){
	# XY is the cross power spectrum
	# ScanIndex is the time index for WG observation (e.g. 1234:1421)
	# chRange is spectral channels range to be used (e.g. 9:65544)
	Xspec <- apply(XY[chRange, ScanIndex], 1, mean)	# CrossCorr
	delay <- delay_search(Xspec)
	return( list(delay=delay, phase=Arg(mean(delay_cal(Xspec, delay)))))
}

unpolDterm <- function( XX, YY, XY, scan, chRange, Tsys, delay, phase ){
	#-------- On- and Off-source power
	XX_on <- mean(XX[chRange, scan$on]); XX_off <- mean(XX[chRange, scan$off])			# Time and Spectral Averaging
	YY_on <- mean(YY[chRange, scan$on]); YY_off <- mean(YY[chRange, scan$off])			# Time and Spectral Averaging
	XY_on <- mean(delay_cal(apply(XY[chRange, scan$on], 1, mean), delay))
	XY_off <- mean(delay_cal(apply(XY[chRange, scan$off], 1, mean), delay))
	#-------- Gain = polaris unit / K
	GainXX <- XX_off / Tsys[1]		# Gx Gx* = <XX*> / <NxNx*>, unit:[polaris/K]
	GainYY <- YY_off / Tsys[2]		# Gy Gy* = <YY*> / <NyNy*>, unit:[polaris/K]
	Gx <- complex(modulus=sqrt(GainXX), argument=0)			# Complex Gain
	Gy <- complex(modulus=sqrt(GainYY), argument=phase)	# Complex Gain
	StokesI <- 0.5*((XX_on - XX_off)/GainXX + (YY_on - YY_off)/GainYY) 	# Antenna temperature of the unpolarized source [K]
	Dx_plus_Dy <- (XY_on - XY_off) / (Gx * Gy* StokesI)
	return(list(StokesI, Dx_plus_Dy))
}
