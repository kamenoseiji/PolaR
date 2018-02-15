# ScanPattern.R
# usage: Rscript ScanPattern.R [Number of SAM45 file] [list of SAM45 file]
# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
#
RPATH <- '~/Programs/PolaR'
FuncList <- c('readSAM45', 'readPolariS', 'PolariCalib', 'Qeff', 'date', 'plotTool')
source(sprintf('%s/loadModule.R', RPATH))
library(RCurl)
scanTime   <- function(MJD_df){ return((MJD_df[[1]] + MJD_df[[2]])/2 )}
funcNum <- length(FuncList)
for( index in 1:funcNum){
    URL <- sprintf("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/%s.R", FuncList[index])
    Err <- try( eval(parse(text = getURL(URL, ssl.verifypeer = FALSE))), silent=FALSE)
}   
if(class(Err) == "try-error"){ loadLocal( RPATH, FuncList ) }
#-------- Parse command-line arguments
parseArg <- function( args ){
    argNum <- length(args)
    lineFreq <- numeric(0)
    for( index in 1:argNum ){
        if(substr(args[index], 1,2) == "-S"){ SpecFile <- substring(args[index], 3)} # *SPEC.Rdata
        if(substr(args[index], 1,2) == "-s"){ ScanPattern <- substring(args[index], 3)} # *Scan.Rdata
        if(substr(args[index], 1,2) == "-D"){ DtermFile <- substring(args[index], 3)}   # *Dcomb.Rdata
        if(substr(args[index], 1,2) == "-B"){ BPFile    <- substring(args[index], 3)}   # *BP.Rdata
        if(substr(args[index], 1,2) == "-W"){ WGFile    <- substring(args[index], 3)}   # WG.Rdata
        if(substr(args[index], 1,2) == "-T"){ ThreshFile <- substring(args[index], 3)}  # *Thresh.Rdata
        if(substr(args[index], 1,2) == "-l"){ lineFreq[1] <- as.numeric(substring(args[index], 3))}
        if(substr(args[index], 1,2) == "-L"){ lineFreq[2] <- as.numeric(substring(args[index], 3))}
    }
    return( list(SpecFile = SpecFile, ScanFile = ScanPattern, DtermFile = DtermFile, BPFile = BPFile, WGFile = WGFile, ThreshFile = ThreshFile, lineFreq = lineFreq) )
}

#-------- Procedures
args <- parseArg(commandArgs(trailingOnly = T))
setwd('.')
load(args$SpecFile)
load(args$ScanFile)
load(args$DtermFile)
load(args$WGFile)
load(args$BPFile)
#-------- Frequency setup
chNum <- dim(on_A00)[1]
freq <- (0:(chNum-1))/chNum* 4.0    # MH
#-------- Smoothed Delay and Phase
numWGentry <- length(WG$mjdSec)
if(numWGentry < 4){
    WG <- data.frame(
        mjdSec  = c( min(WG$mjdSec)-7200, min(WG$mjdSec)-3600, WG$mjdSec, max(WG$mjdSec)+3600, max(WG$mjdSec)+7200),
        delay00 = c( WG$delay00[1], WG$delay00[1], WG$delay00, WG$delay00[numWGentry], WG$delay00[numWGentry]),
        delay01 = c( WG$delay01[1], WG$delay01[1], WG$delay01, WG$delay01[numWGentry], WG$delay01[numWGentry]),
        Vis00 = c( WG$Vis00[1], WG$Vis00[1], WG$Vis00, WG$Vis00[numWGentry], WG$Vis00[numWGentry]),
        Vis01 = c( WG$Vis01[1], WG$Vis01[1], WG$Vis01, WG$Vis01[numWGentry], WG$Vis01[numWGentry]))
}
delay00Fit <- smooth.spline(WG$mjdSec, WG$delay00, spar=0.25)
delay01Fit <- smooth.spline(WG$mjdSec, WG$delay01, spar=0.25)
Re00Fit <- smooth.spline(WG$mjdSec, Re(WG$Vis00), spar=0.25)
Im00Fit <- smooth.spline(WG$mjdSec, Im(WG$Vis00), spar=0.25)
Re01Fit <- smooth.spline(WG$mjdSec, Re(WG$Vis01), spar=0.25)
Im01Fit <- smooth.spline(WG$mjdSec, Im(WG$Vis01), spar=0.25)
#-------- Tsys at each scan
OnIndex <- which(Scan$scanType == 'ON')
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
mitigCH <- c(9:10)
weight <- rep(1, 65536)
Ta00 <- TaCalSpec(on_A00, off_A00, scanTime(onMJD), scanTime(offMJD), Tsys00, weight, mitigCH) / sqrt(mean(D$Rxy02))
Ta01 <- TaCalSpec(on_A01, off_A01, scanTime(onMJD), scanTime(offMJD), Tsys01, weight, mitigCH) / sqrt(mean(D$Rxy13))
Ta02 <- TaCalSpec(on_A02, off_A02, scanTime(onMJD), scanTime(offMJD), Tsys02, weight, mitigCH) * sqrt(mean(D$Rxy02))
Ta03 <- TaCalSpec(on_A03, off_A03, scanTime(onMJD), scanTime(offMJD), Tsys03, weight, mitigCH) * sqrt(mean(D$Rxy13))
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
StokesV02 <- Im(Tx02) - Im(mean(D$XY02))* StokesI02
StokesV13 <- Im(Tx13) - Im(mean(D$XY13))* StokesI13
RHCPspec <- (StokesI02 + StokesV02 + StokesI13 + StokesV13)/4
LHCPspec <- (StokesI02 - StokesV02 + StokesI13 - StokesV13)/4
#-------- Line Flux
chRange <- which( freq > args$lineFreq[1] & freq < args$lineFreq[2] )
RHCPflux <- apply( RHCPspec[chRange,], 2, sum)
LHCPflux <- apply( LHCPspec[chRange,], 2, sum)
#-------- Pointing
cAz <- cEl <- posAz <- posEl <- numeric(0)
for(scan_index in 1:length(onMJD$stopMjd)){
    posAz[scan_index] <- Scan$dAZ[Scan$mjdSec == floor(scanTime(onMJD)[scan_index])]
    posEl[scan_index] <- Scan$dEL[Scan$mjdSec == floor(scanTime(onMJD)[scan_index])]
    cAz[scan_index] <- Scan$AZ[Scan$mjdSec == floor(scanTime(onMJD)[scan_index])]
    cEl[scan_index] <- Scan$EL[Scan$mjdSec == floor(scanTime(onMJD)[scan_index])]
}

beamDF <- data.frame(mjdSrc=scanTime(onMJD), posAz=posAz, posEl=posEl, cAz=cAz, cEl=cEl, RHCPflux=RHCPflux, LHCPflux=LHCPflux)
save(beamDF, file=sprintf('%s.beamDF.Rdata', substr(args$ScanFile, 1, 13)))
#-------- Beam squint fitting
fit <- nls(data=beamDF, formula=(RHCPflux + LHCPflux) ~ a* exp(-0.5*( ( posAz - b)^2 + (posEl -c)^2 )/d), weights=(RHCPflux + LHCPflux)^2, start=list(a=max(beamDF$RHCP+beamDF$LHCP), b=0, c=0, d=400))
beamVar <- summary(fit)$coefficients['d','Estimate']
fit <- lm(data=beamDF, formula=(RHCPflux - LHCPflux)/(RHCPflux + LHCPflux) ~ posAz + posEl, weights=(RHCPflux + LHCPflux)^2)
beamSquint <- as.numeric(summary(fit)$coefficients[c('posAz', 'posEl'),'Estimate'])* beamVar
beamSquintErr <- as.numeric(summary(fit)$coefficients[c('posAz', 'posEl'),'Std. Error'])* beamVar
cat('  AZ     EL      PA   FWHM   Squint(Az) err  Sqruint(El) err  \n')
cat(sprintf('%6.1f  %5.1f  %6.1f  %5.1f     %5.2f  %5.2f     %5.2f  %5.2f \n', median(beamDF$cAz), median(beamDF$cEl), azel2pa(median(beamDF$cAz)/180*pi, median(beamDF$cEl)/180*pi)*180/pi, 2.0* sqrt(beamVar/2.0 / log(2)), beamSquint[1], beamSquintErr[1], beamSquint[2], beamSquintErr[2]))
