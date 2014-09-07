# SpecInteg.R
# usage: Rscript ScanPattern.R [Number of SAM45 file] [list of SAM45 file]
# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
#
library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readSAM45.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readPolariS.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/Qeff.R", ssl.verifypeer = FALSE)))
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
			cat(sprintf('New P=%d S=%d E=%d ', prevFileIndex, startFileIndex, endFileIndex))
			if( postfix == 'C'){
				XP <- readPolariS_X(sprintf('%s.%s.%02d', prefix[startFileIndex], postfix, IF_index))
			} else {
				XP <- readPolariS(sprintf('%s.%s.%02d', prefix[startFileIndex], 'A', IF_index))
			}
		}
		if( endFileIndex > startFileIndex){
			cat(sprintf('Cont P=%d S=%d E=%d ', prevFileIndex, startFileIndex, endFileIndex))
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
		stopIndex  <- MJD[[2]][scanIndex] - prefix2MJDsec(prefix[startFileIndex]) + 1
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
Year <- substr(strsplit(SAM45File[1], '\\.')[[1]][5], 1, 4)
P00fileList <- system(  sprintf('ls %s*.P.00', Year), intern=T )
for(index in 1:length(P00fileList)){
	prefix[index] <- substr(P00fileList[index], 1, 13)
}

#-------- List IF ID of PolariS data
P00fileList <- system(  sprintf('ls %s.P.*', prefix[1]), intern=T )
IF_ID <- integer(0)
for(index in 1:length(P00fileList)){
	IF_ID[index] <- as.integer(strsplit(P00fileList[index] , '\\.')[[1]][3])
}

#-------- Produce Scan data frame
for(fileIndex in 1:length(SAM45File)){
	tempScan <- scanPattern(SAM45File[fileIndex], prefix, IF_ID)
	if(fileIndex == 1){	Scan <- tempScan}
	else { Scan <- rbind(Scan, tempScan)}
}
#-------- Tsys
Scan <- scanTsys(Scan, 290.0)
save(Scan, file=sprintf("%s.Scan.Rdata", prefix[1]))

#-------- Plot
