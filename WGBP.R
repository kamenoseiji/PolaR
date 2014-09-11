# WGBP.R : Generate complex bandpass table using WG scan
# usage: Rscript WGBP.R [thresh] [prefix list]
# e.g. Rscript WGBP.R 15 2014107010610 2014107013610 2014107020610 2014107023610 2014107030610 2014107033610 2014107040610 2014107043610 2014107050610 2014107053610 2014107060610 2014107063610 2014107070610 2014107073610 2014107080610 2014107083610 2014107090610
 
#
library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readPolariS.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/date.R", ssl.verifypeer = FALSE)))
setwd('.')

#-------- Function to filter WG scans
scanGap <- function( XP, thresh, gapThresh ){
	WG_index <- which(Mod(XP$C00) > thresh)	# filter WG scans
	scanGap <- which( diff(WG_index) > gapThresh)
	scanEnd <- c(WG_index[scanGap], max(WG_index))
	scanStart <- c(min(WG_index), WG_index[scanGap+1])
	scanNum <- length(scanStart); XPamp <- numeric(0)
	for(scan_index in 1:scanNum){
		scanRange <- scanStart[scan_index]:scanEnd[scan_index]
		XPamp[scan_index] <- Mod( mean(XP$C00[scanRange]) + mean(XP$C00[scanRange]) )
	}
	return(data.frame( startMJD=XP$mjdSec[scanStart], endMJD=XP$mjdSec[scanEnd], XPamp = XPamp))
}

#-------- Function to pick prefix that covers specified MJD
#findPrefix <- function(mjdSec, prefix){
	#mjdPrefix <- as.numeric(lapply(prefix, FUN=prefix2MJDsec))
	#index <- which(mjdPrefix <= mjdSec)
	#if( length(index) == 0){	return(-1)}
	#return(max(index))
#}

#-------- Function to Generate BP table
makeBP <- function(scanXP, prefix){
	BPscanIndex <- which.max(scanXP$XPamp)
	file_index <- findPrefix(scanXP$startMJD[BPscanIndex], prefix)
	endPoint   <- min( c(scanXP$endMJD[BPscanIndex] - prefix2MJDsec(prefix[file_index]) + 1, 1800) )
	startPoint <- scanXP$startMJD[BPscanIndex] - prefix2MJDsec(prefix[file_index]) + 1
	C00 <- apply(readPolariS_X(sprintf('%s.C.%02d', prefix[file_index], 0))[,startPoint:endPoint], 1, mean)
	C01 <- apply(readPolariS_X(sprintf('%s.C.%02d', prefix[file_index], 1))[,startPoint:endPoint], 1, mean)
	BP <- data.frame( BP00 = BPtable(C00, length(C00)/6), BP01 =  BPtable(C01, length(C01)/6) )
	return(BP)
}

#-------- Procedures
args <- commandArgs()
thresh <- as.numeric(args[6])
gapThresh <- as.numeric(args[7])
prefix <- args[8:length(args)]
XPfname <- sprintf('%s.XP.Rdata', prefix[1])

#-------- Load XP and identify WG scans
load(XPfname)
scanXP <- scanGap(XP, thresh, gapThresh )

#-------- Generate Complex BP table
BP <- makeBP(scanXP, prefix)
save( BP, scanXP, file=sprintf("%s.BP.Rdata", prefix[1]))

