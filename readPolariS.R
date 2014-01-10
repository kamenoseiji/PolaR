readPolariS <- function(fname, chnum=131072){
	head_size <- 128						# 128-byte Header
	file_size <- file.info(fname)$size		# File size [byte]
	ipnum <- (file_size - head_size) / (4* chnum)	# 1 record = sizeof(float) x chnum
	file_ptr <- file(fname, "rb")
	header <- readBin(file_ptr, what=integer(), size=4, n=32)
	temparray <- array(readBin(file_ptr, what=numeric(), n=chnum* ipnum, size=4), dim=c(chnum, ipnum))
	close(file_ptr)
	return(temparray)
}

readBitDist <- function(fname, levelnum=16){	# Default = 4bit 16 levels
	head_size <- 128						# 128-byte Header
	file_size <- file.info(fname)$size		# File size [byte]
	ipnum <- (file_size - head_size) / (4* levelnum)	# 1 record = sizeof(int) x levelnum
	file_ptr <- file(fname, "rb")
	header <- readBin(file_ptr, what=integer(), size=4, n=32)
	temparray <- array(readBin(file_ptr, what=integer(), n=levelnum* ipnum, size=4), dim=c(levelnum, ipnum))
	close(file_ptr)
	return(temparray)
}

readPolariS_X <- function(fname, chnum=131072){
	head_size <- 128						# 128-byte Header
	file_size <- file.info(fname)$size		# File size [byte]
	ipnum <- (file_size - head_size) / (8* chnum)	# 1 record = sizeof(float) x chnum x (re, im)
	file_ptr <- file(fname, "rb")
	header <- readBin(file_ptr, what=integer(), size=4, n=head_size / 4)
	tempvec <- readBin(file_ptr, what=numeric(), n=2*chnum* ipnum, size=4)
	close(file_ptr)
	even_index <- (0:(chnum* ipnum-1))*2; odd_index <- even_index + 1; even_index <- odd_index + 1
	tempcomplex <- complex(real = tempvec[odd_index], imaginary=tempvec[even_index])  # real and imaginary part
	return(array(tempcomplex, dim=c(chnum, ipnum)) )
}

readPolariS_Scan <- function(fname, scan, chnum=131072){
	head_size <- 128						# 128-byte Header
	file_size <- file.info(fname)$size		# File size [byte]
	ipnum <- (file_size - head_size) / (4* chnum)	# 1 record = sizeof(float) x chnum
	pattern <- rep(0, ipnum); pattern[scan] <- 1
	file_ptr <- file(fname, "rb")
	header <- readBin(file_ptr, what=integer(), size=4, n=32)
	spec <- numeric(chnum)
	for(recIndex in 1:ipnum){
		spec <- spec + pattern[recIndex]* readBin(file_ptr, what=numeric(), n=chnum, size=4)
	}
	close(file_ptr)
	return(spec / sum(pattern))
}
		
readPolariS_XScan <- function(fname, scan, chnum=131072){
	head_size <- 128						# 128-byte Header
	file_size <- file.info(fname)$size		# File size [byte]
	ipnum <- (file_size - head_size) / (8* chnum)	# 1 record = sizeof(float) x chnum x (re, im)
	pattern <- rep(0, ipnum); pattern[scan] <- 1
	even_index <- (0:(chnum-1))*2; odd_index <- even_index + 1; even_index <- odd_index + 1
	file_ptr <- file(fname, "rb")
	header <- readBin(file_ptr, what=integer(), size=4, n=head_size / 4)
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

bunch_vec <- function(vector, bunchNum){	# Function to bunch 1-d vector
	pudding <- (bunchNum - length(vector) %% bunchNum) %% bunchNum
	temp <- matrix( append(vector, rep(vector[length(vector)], pudding)), nrow=bunchNum)
	return( apply(temp, 2, mean))
}

wireGridPhaseCorr <- function( prefix ){
	# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
	A00 <- apply(readPolariS(paste(prefix, '.A.00', sep=''))[9:65544,], 2, mean)	# AutoCorr, IF0, X-pol
	A01 <- apply(readPolariS(paste(prefix, '.A.01', sep=''))[9:65544,], 2, mean)	# AutoCorr, IF1, X-pol
	A02 <- apply(readPolariS(paste(prefix, '.A.02', sep=''))[9:65544,], 2, mean)	# AutoCorr, IF2, Y-pol
	A03 <- apply(readPolariS(paste(prefix, '.A.03', sep=''))[9:65544,], 2, mean)	# AutoCorr, IF3, Y-pol
	C00 <- apply(readPolariS_X(paste(prefix, '.C.00', sep=''))[9:65544,], 2, mean) / sqrt(A00* A02)	# CrossCorr
	C01 <- apply(readPolariS_X(paste(prefix, '.C.01', sep=''))[9:65544,], 2, mean) / sqrt(A01* A03)	# CrossCorr
	ScanIndex <- which(Mod(C00) > 0.2)										# Corr. Coeff. is significant
	repC00 <- as.complex(names(which.max(table(round(C00[ScanIndex], 2)))))	# Most dense points
	MScan <- which( Mod( C00 - repC00 ) < 0.1 )								# Set of the +45-deg measurements
	repC00 <- as.complex(names(which.max(table(round(C00[setdiff(ScanIndex, MScan) ], 2)))))	# 2nd dense points
	PScan <- which( Mod( C00 - repC00 ) < 0.1 )								# Set of the -45-deg measurements
	return( c(Arg(mean(c(C00[PScan], -C00[MScan]))),  Arg(mean(c(C01[PScan], -C01[MScan])))) )							# Phase difference [rad]
}

unpolDterm <- function( prefix, PhsDiff ){
	# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
	# PhsDiff is the phase difference between X and Y polarizations (determined by wireGridPhaseCorr()
	load(paste(prefix, '.scan', sep=''))		# Load scan pattern and Tsys
	#-------- Spectra of unpolarized source
	A00 <- readPolariS(paste(prefix, '.A.00', sep=''))[9:65544,]	# Autocorr CH00
	A01 <- readPolariS(paste(prefix, '.A.01', sep=''))[9:65544,]	# Autocorr CH00
	A02 <- readPolariS(paste(prefix, '.A.02', sep=''))[9:65544,]	# Autocorr CH02
	A03 <- readPolariS(paste(prefix, '.A.03', sep=''))[9:65544,]	# Autocorr CH02
	C00 <- readPolariS_X(paste(prefix, '.C.00', sep=''))[9:65544,]	# Cross Correlation
	C01 <- readPolariS_X(paste(prefix, '.C.01', sep=''))[9:65544,]	# Cross Correlation
	#-------- On- and Off-source power
	A00_On <- mean(A00[,onScan]); A00_Off <- mean(A00[,offScan])			# Time and Spectral Averaging
	A01_On <- mean(A01[,onScan]); A01_Off <- mean(A01[,offScan])			# Time and Spectral Averaging
	A02_On <- mean(A02[,onScan]); A02_Off <- mean(A02[,offScan])			# Time and Spectral Averaging
	A03_On <- mean(A03[,onScan]); A03_Off <- mean(A03[,offScan])			# Time and Spectral Averaging
	C00_On <- mean(C00[,onScan]); C00_Off <- mean(C00[,offScan])			# Time and Spectral Averaging
	C01_On <- mean(C01[,onScan]); C01_Off <- mean(C01[,offScan])			# Time and Spectral Averaging
	#-------- Gain = polaris unit / K
	GainXX <- c(A00_Off / mean(Tsys[[1]][offScan]), A01_Off / mean(Tsys[[2]][offScan]))		# Gx Gx* = <XX*> / <NxNx*>, unit:[polaris/K]
	GainYY <- c(A02_Off / mean(Tsys[[3]][offScan]), A03_Off / mean(Tsys[[4]][offScan]))		# Gy Gy* = <YY*> / <NyNy*>, unit:[polaris/K]
	Gx <- c(complex(modulus=sqrt(GainXX[1]), argument=0), 	complex(modulus=sqrt(GainXX[2]), argument=0))	# Complex Gain
	Gy <- c(complex(modulus=sqrt(GainYY[1]), argument=-PhsDiff[1]), complex(modulus=sqrt(GainYY[2]), argument=-PhsDiff[2]))	# Complex Gain
	StokesI <- c((A00_On - A00_Off)/GainXX[1] + (A02_On - A02_Off)/GainYY[1], 	# Antenna temperature of the unpolarized source [K]
				 (A01_On - A01_Off)/GainXX[2] + (A03_On - A03_Off)/GainYY[2])
	Dx_plus_Dy <- 2.0* c((C00_On - C00_Off) / (Gx[1] * Conj(Gy[1])* StokesI[1]), (C01_On - C01_Off) / (Gx[2] * Conj(Gy[2])* StokesI[2]))
	return(Dx_plus_Dy)
}
