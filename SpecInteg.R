# SpecInteg.R
# usage: Rscript SpecInteg.R [Scan.Rdata file name]
# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
#
RPATH <- '~/Programs/PolaR'
FuncList <- c('readSAM45', 'readPolariS', 'date')
source(sprintf('%s/loadModule.R', RPATH))
library(RCurl)

funcNum <- length(FuncList)
for( index in 1:funcNum){
    URL <- sprintf("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/%s.R", FuncList[index])
    Err <- try( eval(parse(text = getURL(URL, ssl.verifypeer = FALSE))), silent=FALSE)
}   
if(class(Err) == "try-error"){ loadLocal( RPATH, FuncList ) }

setwd('.')
IntegCommand <- "/usr/custom/bin/SpecInteg"

#-------- Function to Get PolariS file information
GetChnumRecnum <- function(fname, postfix){
    #-------- Read file header
    file_ptr <- file(fname, "rb")
    header <- readBin(file_ptr, what=integer(), size=4, n=32)
    close(file_ptr)
    #
    #-------- PolariS Header Dictionary
    head_size <- 128
    numCH_index <- 16
    numBit_index <- 11
    chnum <- header[numCH_index]
    levelnum <- 2^header[numBit_index]
    byteperrec <- 0
    if( postfix == 'A' ){   byteperrec <- 4* chnum }
    if( postfix == 'C' ){   byteperrec <- 8* chnum }
    if( postfix == 'P' ){   byteperrec <- 4*levelnum }
    if( byteperrec == 0 ){    return(list(chnum=0, ipnum=0))}
    file_size <- file.info(fname)$size - head_size
    return( list(chnum=chnum, ipnum=file_size / byteperrec) )
}
#-------- Function to produce scan pattern
scanSegment <- function( mjdSec ){
	scanStart <- c(min(mjdSec), mjdSec[which(diff(mjdSec) > 1) + 1])
	scanEnd <- c(mjdSec[which( diff(mjdSec) > 1)], max(mjdSec))
	return( data.frame( startMjd=scanStart, stopMjd=scanEnd ))
}

#-------- Function to integrate spectra referring scan pattern
integSegment <- function( prefix, chnum, ipnum, postfix, IF_index, MJD ){
    #-------- Loop for Scan
	for(scanIndex in 1:length(MJD[[1]])){
		startFileIndex <-findPrefix(MJD[[1]][scanIndex], prefix); endFileIndex <- findPrefix(MJD[[2]][scanIndex], prefix)
        integSec <- MJD[[2]][scanIndex] - MJD[[1]][scanIndex] + 1
        fileNum <- endFileIndex - startFileIndex + 1        # Number of files in the scan
        # cat(sprintf('Scan%d File=%d-%d integ=%d sec\n', 1, startFileIndex, endFileIndex, integSec))
        #-------- Loop for file
        fileCounter <- 0
        remainingIntegSec <- integSec
        spec <- numeric(0)
        while( remainingIntegSec > 0 ){
            startIndex <- max(0, MJD[[1]][scanIndex] - prefix2MJDsec(prefix[startFileIndex + fileCounter]) )
            stopIndex  <- min( startIndex + remainingIntegSec - 1, ipnum[startFileIndex] - 1)
            fileName <- sprintf('%s.%s.%02d', prefix[startFileIndex + fileCounter], postfix, IF_index)
            command_text <- sprintf('%s %s %d %d', IntegCommand, fileName, startIndex, stopIndex)
		    cat(sprintf("SCAN[%d]: MJD range=(%10.0f, %10.0f) %s scanRange=(%d, %d)\n", scanIndex, MJD[[1]][scanIndex], MJD[[2]][scanIndex], prefix[startFileIndex], startIndex, stopIndex))
            # cat(command_text); cat('\n')
            #-------- Throw time-integration process
            system(command_text, wait=T)
            remainingIntegSec <- remainingIntegSec - (stopIndex - startIndex + 1)
            #-------- Read Time-integrated spectrum
            file_size <- file.info('tmp.spec')$size
            tmpSpec <- readBin('tmp.spec', what=numeric(), n=file_size/4, size=4)
            if(postfix == 'C'){
                even_index <- (0:(chnum-1))*2; odd_index <- even_index + 1; even_index <- odd_index + 1
                tmpSpec <- complex(real = tmpSpec[odd_index], imaginary=tmpSpec[even_index])
            }
            if( fileCounter == 0 ){
                singleScanSpec <- tmpSpec
            } else {
                singleScanSpec <- singleScanSpec + tmpSpec
            }
            fileCounter <- fileCounter + 1
        }
		if(scanIndex == 1){ multiScanSpec <- singleScanSpec / integSec }
		if(scanIndex >  1){
            multiScanSpec <- append(multiScanSpec, (singleScanSpec / integSec))
            #cat(sprintf('ScanIndex=%d : Current length of SPEC = %d\n', scanIndex, length(multiScanSpec)))
        }
    }
    #cat( length(MJD[[1]]) )
	return(matrix(multiScanSpec, ncol=length(MJD[[1]])))
}

#-------- Procedures
args <- commandArgs(trailingOnly = T)
prefix <- character(0)
scanFile <- args[1]
load(scanFile)

#-------- Scan Segments
R_index   <- which(Scan$scanType == 'R');  RMJD   <- scanSegment(Scan$mjdSec[R_index])
ON_index  <- which(Scan$scanType == 'ON'); onMJD  <- scanSegment(Scan$mjdSec[ON_index])
OFF_index <- which(Scan$scanType == 'OFF');offMJD <- scanSegment(Scan$mjdSec[OFF_index])

#-------- List prefix of PolariS data
Year <- mjd2doy(Scan$mjdSec[1])[[1]]
P00fileList <- system(  sprintf('ls %s*.P.00', Year), intern=T )
ipnum <- integer(length(P00fileList))
for(index in 1:length(P00fileList)){
	prefix[index] <- substr(P00fileList[index], 1, 13)
    AfileSize <- GetChnumRecnum(sprintf('%s.A.00', prefix[index]), 'A')
    ipnum[index] <- AfileSize$ipnum
}
chNum <- AfileSize$chnum
#
#-------- Integrate Segments
on_C00  <- integSegment(prefix, chNum, ipnum, 'C', 0, onMJD )
on_C01  <- integSegment(prefix, chNum, ipnum, 'C', 1, onMJD )
on_A00  <- integSegment(prefix, chNum, ipnum, 'A', 0, onMJD )
on_A01  <- integSegment(prefix, chNum, ipnum, 'A', 1, onMJD )
on_A02  <- integSegment(prefix, chNum, ipnum, 'A', 2, onMJD )
on_A03  <- integSegment(prefix, chNum, ipnum, 'A', 3, onMJD )
off_C00 <- integSegment(prefix, chNum, ipnum, 'C', 0, offMJD )
off_C01 <- integSegment(prefix, chNum, ipnum, 'C', 1, offMJD )
off_A00 <- integSegment(prefix, chNum, ipnum, 'A', 0, offMJD )
off_A01 <- integSegment(prefix, chNum, ipnum, 'A', 1, offMJD )
off_A02 <- integSegment(prefix, chNum, ipnum, 'A', 2, offMJD )
off_A03 <- integSegment(prefix, chNum, ipnum, 'A', 3, offMJD )
#-------- Save into file
StartUTC <- mjd2doy(min(Scan$mjdSec[ON_index]))
fileName <- sprintf("%04d%03d%02d%02d%02d.SPEC.Rdata", StartUTC$year, StartUTC$doy, StartUTC$hour, StartUTC$min, StartUTC$sec)
save(onMJD, on_C00, on_C01, on_A00, on_A01, on_A02, on_A03, offMJD, off_C00, off_C01, off_A00, off_A01, off_A02, off_A03, file=fileName)
cat('Segment-integrated spectra (uncal) are saved into '); cat(fileName); cat('\n')
