# ScanPattern.R
# usage: Rscript ScanPattern.R [Number of SAM45 file] [list of SAM45 file]
# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
#
RPATH <- '~/Programs/PolaR'
FuncList <- c('readSAM45', 'readPolariS', 'PolariCalib', 'Qeff', 'date', 'plotTool')
source(sprintf('%s/loadModule.R', RPATH))
library(RCurl)

funcNum <- length(FuncList)
for( index in 1:funcNum){
    URL <- sprintf("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/%s.R", FuncList[index])
    Err <- try( eval(parse(text = getURL(URL, ssl.verifypeer = FALSE))), silent=FALSE)
}   
if(class(Err) == "try-error"){ loadLocal( RPATH, FuncList ) }
#-------- Parse command-line arguments
parseArg <- function( args ){
    argNum <- length(args)
    smoothWidth <- 128      # Set Default Value
    for( index in 1:argNum ){
        if(substr(args[index], 1,2) == "-S"){ SAM45File <- substring(args[index], 3)}
        if(substr(args[index], 1,2) == "-D"){ DtermFile <- substring(args[index], 3)}
        if(substr(args[index], 1,2) == "-B"){ BPFile    <- substring(args[index], 3)}
        if(substr(args[index], 1,2) == "-W"){ WGFile    <- substring(args[index], 3)}
        if(substr(args[index], 1,2) == "-T"){ ThreshFile <- substring(args[index], 3)}
    }
    return( list(SAM45File = SAM45File, DtermFile = DtermFile, BPFile = BPFile, WGFile = WGFile, ThreshFile = ThreshFile) )
}

#-------- Procedures
args <- parseArg(commandArgs(trailingOnly = T))
setwd('.')
load(args$DtermFile)
load(args$WGFile)
load(args$BPFile)
#-------- Smoothed Delay and Phase
delay00Fit <- smooth.spline(WG$mjdSec, WG$delay00, spar=0.25)
delay01Fit <- smooth.spline(WG$mjdSec, WG$delay01, spar=0.25)
Re00Fit <- smooth.spline(WG$mjdSec, Re(WG$Vis00), spar=0.25)
Im00Fit <- smooth.spline(WG$mjdSec, Im(WG$Vis00), spar=0.25)
Re01Fit <- smooth.spline(WG$mjdSec, Re(WG$Vis01), spar=0.25)
Im01Fit <- smooth.spline(WG$mjdSec, Im(WG$Vis01), spar=0.25)
#
#-------- List prefix of PolariS data
prefix <- character(0)
Year <- substr(strsplit(args$SAM45File, '\\.')[[1]][5], 1, 4)
A00fileList <- system(  sprintf('ls %s*.A.00', Year), intern=T )
ipnum <- integer(length(A00fileList))
for(index in 1:length(A00fileList)){
	prefix[index] <- substr(A00fileList[index], 1, 13)
	AfileSize <- GetChnumRecnum(sprintf('%s.A.00', prefix[index]), 'A')
	ipnum[index] <- AfileSize$ipnum
}
chNum <- AfileSize$chnum
freq <- (0:(chNum-1))/chNum* 4.0    # MHz
chSep <- 4.0 / chNum
mitigCH <- c(23359, 32769)
flagCH <- unique(c(1, 2, 4, mitigCH))
weight <- rep(1, chNum); weight[flagCH] <- 0.0
#-------- List IF ID of PolariS data
A00fileList <- system(  sprintf('ls %s.A.*', prefix[1]), intern=T )
IF_ID <- integer(0)
for(index in 1:length(A00fileList)){
	IF_ID[index] <- as.integer(strsplit(A00fileList[index] , '\\.')[[1]][3])
}
#-------- Produce Scan data frame
tmpScan <- scanPattern(args$SAM45File, prefix, IF_ID, args$ThreshFile); Scan <- scanTsys(tmpScan$scanDF, 280.0)
OnIndex <- which( Scan$scanType == 'ON');  onMJD  <- scanSegment(Scan$mjdSec[OnIndex])
OfIndex <- which( Scan$scanType == 'OFF'); offMJD <- scanSegment(Scan$mjdSec[OfIndex])
D_index <- which( D$Gxy02 + D$Gxy13 > 0.9* median( D$Gxy02 + D$Gxy13 )) # Flag pointing-error out
Rxy02 <- mean(D$Rxy02[D_index])
Rxy13 <- mean(D$Rxy13[D_index])
Dxy02 <- mean(D$XY02[D_index])
Dxy13 <- mean(D$XY13[D_index])
#--------
#AfileIndex <- seq(findPrefix(min(SAM45df$mjd_st), prefix), findPrefix(max(SAM45df$mjd_st), prefix))
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
#-------- Tsys at each scan
Tsys00 <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$Tsys00[OnIndex], spar=0.5), scanTime(onMJD))$y
Tsys01 <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$Tsys01[OnIndex], spar=0.5), scanTime(onMJD))$y
Tsys02 <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$Tsys02[OnIndex], spar=0.5), scanTime(onMJD))$y
Tsys03 <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$Tsys03[OnIndex], spar=0.5), scanTime(onMJD))$y
#-------- Parallactic angle
AZ <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$AZ[OnIndex], spar=0.5), scanTime(onMJD))$y
EL <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$EL[OnIndex], spar=0.5), scanTime(onMJD))$y
Pang <- -azel2pa(AZ, EL) + EL*pi/180 - pi/2
cs <- cos(Pang)
sn <- sin(Pang)
#-------- Amplitude calibration of Autocorr
Ta00 <- TaCalSpec(on_A00, off_A00, scanTime(onMJD), scanTime(offMJD), Tsys00, weight, mitigCH) / sqrt(Rxy02)
Ta01 <- TaCalSpec(on_A01, off_A01, scanTime(onMJD), scanTime(offMJD), Tsys01, weight, mitigCH) / sqrt(Rxy13)
Ta02 <- TaCalSpec(on_A02, off_A02, scanTime(onMJD), scanTime(offMJD), Tsys02, weight, mitigCH) * sqrt(Rxy02)
Ta03 <- TaCalSpec(on_A03, off_A03, scanTime(onMJD), scanTime(offMJD), Tsys03, weight, mitigCH) * sqrt(Rxy13)
#-------- Delay, phase, bandpass, and amplitude calibration for CrossCorr
Tx02 <- TxCalSpec(
    DelayPhaseCal( BPphsCal(on_C00, BP$BP00), scanTime(onMJD), delay00Fit, Re00Fit, Im00Fit),
    DelayPhaseCal( BPphsCal(off_C00, BP$BP00), scanTime(offMJD), delay00Fit, Re00Fit, Im00Fit),
    off_A00, off_A02, scanTime(onMJD), scanTime(offMJD), sqrt(Tsys00 * Tsys02), weight, mitigCH)
Tx13 <- TxCalSpec(
    DelayPhaseCal( BPphsCal(on_C01, BP$BP01), scanTime(onMJD), delay01Fit, Re01Fit, Im01Fit),
    DelayPhaseCal( BPphsCal(off_C01, BP$BP01), scanTime(offMJD), delay01Fit, Re01Fit, Im01Fit),
    off_A01, off_A03, scanTime(onMJD), scanTime(offMJD), sqrt(Tsys01 * Tsys03), weight, mitigCH)
StokesI02 <- 0.5*(Ta00 + Ta02)
StokesI13 <- 0.5*(Ta01 + Ta03)
StokesV02 <- Im(Tx02) - Im(Dxy02)* StokesI02
StokesV13 <- Im(Tx13) - Im(Dxy13)* StokesI13
RHCPspec <- (StokesI02 + StokesV02 + StokesI13 + StokesV13)/4
LHCPspec <- (StokesI02 - StokesV02 + StokesI13 - StokesV13)/4
#-------- Line Flux
lineRange <- c(1.58, 1.8)
chRange <- which( freq > lineRange[1] & freq < lineRange[2] )
RHCPflux <- apply( RHCPspec[chRange,], 2, sum)
LHCPflux <- apply( LHCPspec[chRange,], 2, sum)
#-------- Pointing
posAz <- posEl <- numeric(0)
for(scan_index in 1:length(onMJD$stopMjd)){
    posAz[scan_index] <- Scan$dAZ[Scan$mjdSec == floor(scanTime(onMJD)[scan_index])]
    posEl[scan_index] <- Scan$dEL[Scan$mjdSec == floor(scanTime(onMJD)[scan_index])]
}
pointingNum <- floor((length(offMJD$stopMjd) - 1) / 2)
for(scan_index in 1:pointingNum){
    posIndex <- (6*scan_index - 5):(6*scan_index)
    posOff <- Scan$dAZ[Scan$mjdSec == floor(scanTime(offMJD)[scan_index])]
    RDF <- data.frame(
        x = c(posAz[posIndex], -posOff, posOff, 0, 0),
        y = c(posEl[posIndex], 0, 0, -posOff, posOff),
        Z = c(RHCPflux[posIndex], 0, 0, 0, 0))
    LDF <- RDF; LDF$Z <- c(LHCPflux[posIndex], 0, 0, 0, 0)
    Rfit <- nls( formula = Z ~ a* exp( -0.5*( (x - x0)^2 + (y - y0)^2 )/ sigma^2), data = RDF, start=list(a=max(RDF$Z), x0=0.0, y0=0.0, sigma=max(posAz[posIndex])))
    Lfit <- nls( formula = Z ~ a* exp( -0.5*( (x - x0)^2 + (y - y0)^2 )/ sigma^2), data = LDF, start=list(a=max(RDF$Z), x0=0.0, y0=0.0, sigma=max(posAz[posIndex])))
    text_sd <- sprintf('Scan[%d] RHCP: AZ= %5.2f (%4.2f) / EL= %5.2f (%4.2f) Peak=%5.1f FWHM=%4.2f', scan_index, coef(Rfit)['x0'], sqrt(vcov(Rfit)['x0','x0']), coef(Rfit)['y0'], sqrt(vcov(Rfit)['y0','y0']), coef(Rfit)['a'], 2.0*sqrt(2.0* log(2.0))* coef(Rfit)['sigma'])
    cat(text_sd); cat('\n')
    text_sd <- sprintf('Scan[%d] LHCP: AZ= %5.2f (%4.2f) / EL= %5.2f (%4.2f) Peak=%5.1f FWHM=%4.2f', scan_index, coef(Lfit)['x0'], sqrt(vcov(Lfit)['x0','x0']), coef(Lfit)['y0'], sqrt(vcov(Lfit)['y0','y0']), coef(Lfit)['a'], 2.0*sqrt(2.0* log(2.0))* coef(Lfit)['sigma'])
    cat(text_sd); cat('\n')
}
