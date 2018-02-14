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
save(beamDF, file='beamDF.Rdata')

if(0){
pointingNum <- floor((length(offMJD$stopMjd) - 1) / 2)
cat('#  AZ   EL  |R dAZ   err     dEL   err    Peak  FWHM   |L  dAZ   err     dEL   err    Peak  FWHM  Squint(Az,El) \n')
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
    if( scan_index == 1){
        SquintDF <- data.frame(
            EL=cEl[scan_index],
            Raz=coef(Rfit)['x0'], Rel=coef(Rfit)['y0'], Rpeak=coef(Rfit)['a'], RFWHM=2.0*sqrt(2.0* log(2.0))* coef(Rfit)['sigma'],
            Raze=sqrt(vcov(Rfit)['x0','x0']), Rele=sqrt(vcov(Rfit)['y0','y0']),
            Laz=coef(Lfit)['x0'], Lel=coef(Lfit)['y0'], Lpeak=coef(Lfit)['a'], LFWHM=2.0*sqrt(2.0* log(2.0))* coef(Lfit)['sigma'],
            Laze=sqrt(vcov(Lfit)['x0','x0']), Lele=sqrt(vcov(Lfit)['y0','y0']))
    } else {
        SquintDF <- rbind(SquintDF,  data.frame(
            EL=cEl[scan_index],
            Raz=coef(Rfit)['x0'], Rel=coef(Rfit)['y0'], Rpeak=coef(Rfit)['a'], RFWHM=2.0*sqrt(2.0* log(2.0))* coef(Rfit)['sigma'],
            Raze=sqrt(vcov(Rfit)['x0','x0']), Rele=sqrt(vcov(Rfit)['y0','y0']),
            Laz=coef(Lfit)['x0'], Lel=coef(Lfit)['y0'], Lpeak=coef(Lfit)['a'], LFWHM=2.0*sqrt(2.0* log(2.0))* coef(Lfit)['sigma'],
            Laze=sqrt(vcov(Lfit)['x0','x0']), Lele=sqrt(vcov(Lfit)['y0','y0'])))
    }
    text_sd <- sprintf('%d %4.1f %4.1f %5.2f (%4.2f)  %5.2f (%4.2f) %5.1f %4.2f ', scan_index, cAz[scan_index], cEl[scan_index], coef(Rfit)['x0'], sqrt(vcov(Rfit)['x0','x0']), coef(Rfit)['y0'], sqrt(vcov(Rfit)['y0','y0']), coef(Rfit)['a'], 2.0*sqrt(2.0* log(2.0))* coef(Rfit)['sigma'])
    cat(text_sd)
    text_sd <- sprintf('   %5.2f (%4.2f)  %5.2f (%4.2f) %5.1f %4.2f ', coef(Lfit)['x0'], sqrt(vcov(Lfit)['x0','x0']), coef(Lfit)['y0'], sqrt(vcov(Lfit)['y0','y0']), coef(Lfit)['a'], 2.0*sqrt(2.0* log(2.0))* coef(Lfit)['sigma'])
    cat(text_sd)
    text_sd <- sprintf('(%5.2f %5.2f)', coef(Rfit)['x0'] - coef(Lfit)['x0'], coef(Rfit)['y0'] - coef(Lfit)['y0']) 
    cat(text_sd); cat('\n')
}
save(SquintDF, file=sprintf('%s.Squint.Rdata', args$SAM45File))
pdf(sprintf("%s.Squint.pdf", args$SAM45File))
plot( bunch_vec(freq, 32), bunch_vec(RHCPspec[,posIndex[5]], 32), type='s', col='red', ylim=c(-1, max(max(RHCPspec[chRange,]))), xlab='Frequency [MHz]', ylab='Ta* [K]', main=args$SAM45File)
lines(bunch_vec(freq, 32), bunch_vec(LHCPspec[,posIndex[5]], 32), type='s', col='blue')
abline(v=args$lineFreq[1]); abline(v=args$lineFreq[2])
dev.off()
}
