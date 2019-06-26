# ContStokes.R
# usage: Rscript ContStokes.R [Scan.Rdata file name] [SPEC.Rdata file name] [WG.Rdata file name] [BP file name]
#
RPATH <- '~/Programs/PolaR'
FuncList <- c('date', 'PolariCalib')
source(sprintf('%s/loadModule.R', RPATH))
library(RCurl)

funcNum <- length(FuncList)
for( index in 1:funcNum){
    URL <- sprintf("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/%s.R", FuncList[index])
    Err <- try( eval(parse(text = getURL(URL, ssl.verifypeer = FALSE))), silent=FALSE)
}   
if(class(Err) == "try-error"){ loadLocal( RPATH, FuncList ) }

setwd('.')
options(digits = 4)
#-------- ScanTime
scanTime   <- function(MJD_df){ return((MJD_df[[1]] + MJD_df[[2]])/2 ) }
scanIntegT <- function(MJD_df){ return( MJD_df[[2]] - MJD_df[[2]] + 1) }

#-------- BP Cal
BPphsCal <- function( SPEC, BP ){
	BPphs <- complex( modulus=1, argument=-Arg(BP))
	return( SPEC* BPphs )
}

#-------- Function to Calibrate Delay and Phase
DelayPhaseCal <- function( scanSpec, mjdSec, delayFit, ReFit, ImFit){
	temp <- scanSpec
	phase <- atan2(predict(ImFit, mjdSec)$y, predict(ReFit, mjdSec)$y)
	delay <- predict(delayFit, mjdSec)$y
	for(timeIndex in 1:ncol(scanSpec)){
		scanSpec[,timeIndex] <- delayPhase_cal(temp[,timeIndex], delay[timeIndex], -phase[timeIndex])
	}
	return(scanSpec)
}

#-------- On - Off subtraction for total power
TaCal <- function(OnPower, OffPower, Tsys, OnTime, OffTime, TsysTime ){
	if( length(OffTime) > 3){
		basePower <- predict(smooth.spline( OffTime, OffPower, spar=0.5 ), OnTime)$y
	} else {
		basePower <- rep(mean(OffPower),  length(OnTime))
	}
	baseTsys  <- predict(smooth.spline( TsysTime, Tsys, spar=0.5 ), OnTime)$y
	# cat(TsysTime); cat('\n')
	# cat(Tsys); cat('\n')
	cat(baseTsys); cat('\n')
	# cat(OnTime); cat('\n')
	# cat(OnPower); cat('\n')
	# cat(OffPower); cat('\n')
	plot( c(OnTime, OffTime), c(OnPower, OffPower), xlab='MJD [sec]', ylab='Relative Power', main='Channel-Averaged Power', type='n')
	points( OnTime, OnPower, pch=20, col='blue'); # lines( OnTime, OnPower, lwd=0.5, col='blue')
	points( OffTime, OffPower, pch=20, col='cyan')
	points( OnTime, basePower, pch=20, cex=0.5, col='gray' )
	arrows( OnTime, basePower, OnTime, OnPower, length=0.0, col='blue')
	return(baseTsys* (OnPower - basePower) / basePower)
}

#-------- On - Off subtraction for complex power
TxCal <- function(OnXpower, OffXpower, OffPower, Tsys, OnTime, OffTime, TsysTime){
	if( length(OffTime) > 3){
		baseXpower <- complex(
			real = predict(smooth.spline( OffTime, Re(OffXpower), spar=0.5 ), OnTime)$y,
			imaginary = predict(smooth.spline( OffTime, Im(OffXpower), spar=0.5 ), OnTime)$y )
		basePower <- predict(smooth.spline( OffTime, OffPower, spar=0.5 ), OnTime)$y
	} else {
		baseXpower <- rep( mean(OffXpower), length(OnTime))
		basePower <- rep(mean(OffPower),  length(OnTime))
	}
	baseTsys  <- predict(smooth.spline( TsysTime, Tsys, spar=0.5 ), OnTime)$y
	return(baseTsys* (OnXpower - baseXpower) / basePower)
}

#-------- Stokes Parameters
corr2Stokes <- function( XX, YY, XY, Dxy, pang ){
	cs <- cos(2.0* pang);	sn <- sin(2.0* pang)	# pang: parallactic angle [rad]
	StokesI <- 0.5* (XX + YY)
	StokesQ <- 0.5* cs* (XX - YY) - sn* (Re(XY) - StokesI* Re(Dxy))
	StokesU <- 0.5* sn* (XX - YY) + cs* (Re(XY) - StokesI* Re(Dxy))
	StokesV <- Im(XY) - StokesI* Im(Dxy)
	return( data.frame(I=StokesI, Q=StokesQ, U=StokesU, V=StokesV, PA=pang) )
}


#-------- Load Spec and Scan data
# args <- commandArgs(trailingOnly = T)
setwd('.')
load(args[1])	 #Load Scan file
load(args[2])	 #Load SPEC file
load(args[3])	 #Load D-term file
load(args[4])	 #Load delay file
load(args[5])	 #Load BP file

#-------- Smoothed Delay and Phase
delayNum <- length(WG$mjdSec)
if(delayNum < 5){
    WG <- rbind( WG[1,], WG[1,], WG[1,], WG, WG[delayNum,], WG[delayNum,], WG[delayNum,])
    delayNum <- length(WG$mjdSec)
    WG$mjdSec[3] <- WG$mjdSec[1] - 1800
    WG$mjdSec[2] <- WG$mjdSec[1] - 3600
    WG$mjdSec[1] <- WG$mjdSec[1] - 5400
    WG$mjdSec[delayNum -2] <- WG$mjdSec[delayNum] + 1800
    WG$mjdSec[delayNum -1] <- WG$mjdSec[delayNum] + 3600
    WG$mjdSec[delayNum] <- WG$mjdSec[delayNum] + 5400
}
delay00Fit <- smooth.spline(WG$mjdSec, WG$delay00, spar=0.25); delay01Fit <- smooth.spline(WG$mjdSec, WG$delay01, spar=0.25)
Re00Fit <- smooth.spline(WG$mjdSec, Re(WG$Vis00), spar=0.25); Im00Fit <- smooth.spline(WG$mjdSec, Im(WG$Vis00), spar=0.25)
Re01Fit <- smooth.spline(WG$mjdSec, Re(WG$Vis01), spar=0.25); Im01Fit <- smooth.spline(WG$mjdSec, Im(WG$Vis01), spar=0.25)
#
#-------- Initial Parameters
chNum <- dim(on_A00)[1]
chRange <- floor(0.05*chNum):floor(0.95*chNum)

OnIndex <- which(Scan$scanType == 'ON')
OfIndex <- which(Scan$scanType == 'OFF')

D_index <- which( D$Gxy02 + D$Gxy13 > 0.9* median( D$Gxy02 + D$Gxy13 ))	# Flag pointing-error out
Rxy02 <- mean(D$Rxy02[D_index])
Rxy13 <- mean(D$Rxy13[D_index])
Dxy02 <- mean(D$XY02[D_index])
Dxy13 <- mean(D$XY13[D_index])
#-------- Amplitude calibration of Autocorr
fileName <- sprintf("%s.ContPower.pdf", strsplit(args[2], "\\.")[[1]][1])
pdf(fileName)
Ta00 <- TaCal(colSums(on_A00[chRange,]), colSums(off_A00[chRange,]), Scan$Tsys00[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]) / sqrt(Rxy02) #; cat(Scan$Tsys00[OfIndex]); cat('\n') # ; cat(Ta00); cat('\n')
Ta01 <- TaCal(colSums(on_A01[chRange,]), colSums(off_A01[chRange,]), Scan$Tsys01[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]) / sqrt(Rxy13) #; cat(Scan$Tsys01[OfIndex]); cat('\n') # ; cat(Ta01); cat('\n')
Ta02 <- TaCal(colSums(on_A02[chRange,]), colSums(off_A02[chRange,]), Scan$Tsys02[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]) * sqrt(Rxy02) #; cat(Scan$Tsys02[OfIndex]); cat('\n') # ; cat(Ta02); cat('\n')
Ta03 <- TaCal(colSums(on_A03[chRange,]), colSums(off_A03[chRange,]), Scan$Tsys03[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]) * sqrt(Rxy13) #; cat(Scan$Tsys03[OfIndex]); cat('\n') # ; cat(Ta03); cat('\n')

#-------- Delay, phase, bandpass, and amplitude calibration for CrossCorr
Tx02 <- TxCal( colSums( DelayPhaseCal( BPphsCal(on_C00, BP$BP00), scanTime(onMJD), delay00Fit, Re00Fit, Im00Fit)[chRange,]),
               colSums( DelayPhaseCal( BPphsCal(off_C00, BP$BP00), scanTime(offMJD), delay00Fit, Re00Fit, Im00Fit)[chRange,]),
			   sqrt( colSums(off_A00[chRange,])* colSums(off_A02[chRange,]) ),
			   sqrt(Scan$Tsys00[OfIndex]*  Scan$Tsys02[OfIndex]),
			   scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex])
Tx13 <- TxCal( colSums( DelayPhaseCal( BPphsCal(on_C01, BP$BP01), scanTime(onMJD), delay01Fit, Re01Fit, Im01Fit)[chRange,]),
               colSums( DelayPhaseCal( BPphsCal(off_C01, BP$BP01), scanTime(offMJD), delay01Fit, Re01Fit, Im01Fit)[chRange,]),
			   sqrt( colSums(off_A01[chRange,])* colSums(off_A03[chRange,]) ), 
			   sqrt(Scan$Tsys01[OfIndex]*  Scan$Tsys03[OfIndex]),
			   scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex])
#-------- Parallactic angle
AZ <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$AZ[OnIndex], spar=0.5), scanTime(onMJD))$y
EL <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$EL[OnIndex], spar=0.5), scanTime(onMJD))$y
Pang <- -azel2pa(AZ, EL) + EL*pi/180 - pi/2 
#-------- Correlation -> Stokes Parameters
if( length( scanTime(onMJD) ) > 3){
	Stokes00 <- corr2Stokes( predict(smooth.spline(scanTime(onMJD), Ta00, spar=0.8), scanTime(onMJD))$y, predict(smooth.spline(scanTime(onMJD), Ta02, spar=0.8), scanTime(onMJD))$y, Tx02, Dxy02, Pang)
	Stokes01 <- corr2Stokes( predict(smooth.spline(scanTime(onMJD), Ta01, spar=0.8), scanTime(onMJD))$y, predict(smooth.spline(scanTime(onMJD), Ta03, spar=0.8), scanTime(onMJD))$y, Tx13, Dxy13, Pang)
} else {
	Stokes00 <- corr2Stokes(Ta00, Ta02, Tx02, Dxy02, Pang)
	Stokes01 <- corr2Stokes(Ta01, Ta03, Tx13, Dxy13, Pang)
}

Stokes00$p <- sqrt(Stokes00$Q^2 + Stokes00$U^2); Stokes00$EVPA <- atan2(Stokes00$U, Stokes00$Q)*90/pi; Stokes00$mjdSec <- scanTime(onMJD); Stokes00$AZ <- AZ; Stokes00$EL <- EL
Stokes01$p <- sqrt(Stokes01$Q^2 + Stokes01$U^2); Stokes01$EVPA <- atan2(Stokes01$U, Stokes01$Q)*90/pi; Stokes01$mjdSec <- scanTime(onMJD); Stokes01$AZ <- AZ; Stokes01$EL <- EL

#cat(sprintf('AZ=%4.1f  EL=%4.1f  PA=%4.1f (deg)\n', AZ, EL, Pang*180/pi))
fileName <- sprintf("%s.STOKES.Rdata", strsplit(args[2], "\\.")[[1]][1])
save(Stokes00, Stokes01, file=fileName)
Stokes00
cat( mean(Stokes00$I), mean(Stokes00$Q), mean(Stokes00$U), mean(Stokes00$V), sqrt( mean(Stokes00$Q)^2 + mean(Stokes00$U)^2), atan2( mean(Stokes00$U), mean(Stokes00$Q))*90/pi); cat('\n')
cat( sd(Stokes00$I), sd(Stokes00$Q), sd(Stokes00$U), sd(Stokes00$V), sd(Stokes00$p), sd(Stokes00$EVPA)); cat('\n')
Stokes01
cat( mean(Stokes01$I), mean(Stokes01$Q), mean(Stokes01$U), mean(Stokes01$V), sqrt( mean(Stokes01$Q)^2 + mean(Stokes01$U)^2), atan2( mean(Stokes01$U), mean(Stokes01$Q))*90/pi); cat('\n')
cat( sd(Stokes01$I), sd(Stokes01$Q), sd(Stokes01$U), sd(Stokes01$V), sd(Stokes01$p), sd(Stokes01$EVPA)); cat('\n')
cat('Stokes parameters (D-term calibrated) was saved into '); cat(fileName); cat('\n')
dev.off()
