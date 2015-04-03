# SpecInteg.R
# usage: Rscript SpecInteg.R [Scan.Rdata file name]
# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
#
library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readSAM45.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readPolariS.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/mjd.R", ssl.verifypeer = FALSE)))
setwd('.')
#-------- Function to produce scan pattern
scanSegment <- function( mjdSec ){
	scanStart <- c(min(mjdSec), mjdSec[which(diff(mjdSec) > 1) + 1])
	scanEnd <- c(mjdSec[which( diff(mjdSec) > 1)], max(mjdSec))
	return( data.frame( startMjd=scanStart, stopMjd=scanEnd ))
}

#-------- Function to integrate spectra referring scan pattern
integSegment <- function( prefix, postfix, IF_index, MJD ){
	prevFileIndex  <- 0
	for(scanIndex in 1:length(MJD[[1]])){
		startFileIndex <- findPrefix(MJD[[1]][scanIndex], prefix)
		endFileIndex   <- findPrefix(MJD[[2]][scanIndex], prefix)
		if( startFileIndex != prevFileIndex){
			cat(sprintf('New P=%d S=%d E=%d \n', prevFileIndex, startFileIndex, endFileIndex))
			if( postfix == 'C'){
				XP <- readPolariS_X(sprintf('%s.%s.%02d', prefix[startFileIndex], postfix, IF_index))
			} else {
				XP <- readPolariS(sprintf('%s.%s.%02d', prefix[startFileIndex], 'A', IF_index))
			}
		}
		if( endFileIndex > startFileIndex){
			cat(sprintf('Cont P=%d S=%d E=%d \n', prevFileIndex, startFileIndex, endFileIndex))
			if( postfix == 'C'){
				temp <- readPolariS_X(sprintf('%s.%s.%02d', prefix[endFileIndex], postfix, IF_index))
			} else {
				temp <- readPolariS(sprintf('%s.%s.%02d', prefix[endFileIndex], 'A', IF_index))
			}
			cat( dim(XP)); cat(' '); cat(dim(temp))
			XP <- cbind(XP, temp)
		}
		prevFileIndex <- startFileIndex
		startIndex <- MJD[[1]][scanIndex] - prefix2MJDsec(prefix[startFileIndex]) + 1
		stopIndex  <- MJD[[2]][scanIndex] - prefix2MJDsec(prefix[startFileIndex])
		# stopIndex  <- MJD[[2]][scanIndex] - prefix2MJDsec(prefix[startFileIndex]) + 1
		cat(sprintf("SCAN[%d]: MJD range=(%10.0f, %10.0f)  prefix-%s  scanRange=(%d, %d)\n", scanIndex, MJD[[1]][scanIndex], MJD[[2]][scanIndex], prefix[startFileIndex], startIndex, stopIndex))
		if(scanIndex == 1){
			spec <- apply(XP[,startIndex:stopIndex], 1, mean)
		} else {
			spec <- append(spec, apply(XP[,startIndex:stopIndex], 1, mean))
		}
	}
	return(matrix(spec, ncol=length(MJD[[1]])))
}

#-------- Procedures
args <- commandArgs()
prefix <- character(0)
scanFile <- args[6]
load(scanFile)

#-------- Scan Segments
R_index   <- which(Scan$scanType == 'R');  RMJD   <- scanSegment(Scan$mjdSec[R_index])
ON_index  <- which(Scan$scanType == 'ON'); onMJD  <- scanSegment(Scan$mjdSec[ON_index])
OFF_index <- which(Scan$scanType == 'OFF');offMJD <- scanSegment(Scan$mjdSec[OFF_index])

#-------- List prefix of PolariS data
Year <- mjd2doy(Scan$mjdSec[1])[[1]]
P00fileList <- system(  sprintf('ls %s*.P.00', Year), intern=T )
for(index in 1:length(P00fileList)){
	prefix[index] <- substr(P00fileList[index], 1, 13)
}
postFix <- ''
if(file.exists(sprintf('%s.C.%02dB', prefix[1], 0))){ postFix <- 'B' }
chnum <- GetChNum(sprintf('%s.C.%02d%s', prefix[1], 0, postFix))
#
#-------- Integrate Segments
on_C00 <- integSegment(prefix, 'C', 0, onMJD ); off_C00 <- integSegment(prefix, 'C', 0, offMJD ) # ; R_C00 <- integSegment(prefix, 'C', 0, RMJD )
on_C01 <- integSegment(prefix, 'C', 1, onMJD ); off_C01 <- integSegment(prefix, 'C', 1, offMJD ) # ; R_C01 <- integSegment(prefix, 'C', 1, RMJD )
on_A00 <- integSegment(prefix, 'A', 0, onMJD ); off_A00 <- integSegment(prefix, 'A', 0, offMJD ) # ; R_A00 <- integSegment(prefix, 'A', 0, RMJD )
on_A01 <- integSegment(prefix, 'A', 1, onMJD ); off_A01 <- integSegment(prefix, 'A', 1, offMJD ) # ; R_A01 <- integSegment(prefix, 'A', 1, RMJD )
on_A02 <- integSegment(prefix, 'A', 2, onMJD ); off_A02 <- integSegment(prefix, 'A', 2, offMJD ) # ; R_A02 <- integSegment(prefix, 'A', 2, RMJD )
on_A03 <- integSegment(prefix, 'A', 3, onMJD ); off_A03 <- integSegment(prefix, 'A', 3, offMJD ) # ; R_A03 <- integSegment(prefix, 'A', 3, RMJD )

#-------- Save into file
StartUTC <- mjd2doy(min(Scan$mjdSec[ON_index]))
fileName <- sprintf("%04d%03d%02d%02d%02d.SPEC.Rdata", StartUTC$year, StartUTC$doy, StartUTC$hour, StartUTC$min, StartUTC$sec)
#save(onMJD, on_C00, on_C01, on_A00, on_A01, on_A02, on_A03, offMJD, off_C00, off_C01, off_A00, off_A01, off_A02, off_A03, RMJD, R_C00, R_C01, R_A00, R_A01, R_A02, R_A03, file=fileName)
save(onMJD, on_C00, on_C01, on_A00, on_A01, on_A02, on_A03, offMJD, off_C00, off_C01, off_A00, off_A01, off_A02, off_A03, file=fileName)
cat('Segment-integrated spectra (uncal) are saved into '); cat(fileName); cat('\n')
