# ContStokes.R
# usage: Rscript ContStokes.R [Scan.Rdata file name] [SPEC.Rdata file name] [WG.Rdata file name] [BP file name]
#
library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/mjd.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/delaysearch.R", ssl.verifypeer = FALSE)))
setwd('.')

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
	basePower <- predict(smooth.spline( OffTime, OffPower, spar=0.8 ), OnTime)$y
	baseTsys  <- predict(smooth.spline( TsysTime, Tsys, spar=0.8 ), OnTime)$y
	# cat(TsysTime); cat('\n')
	# cat(Tsys); cat('\n')
	cat(baseTsys); cat('\n')
	# cat(OnTime); cat('\n')
	# cat(OnPower); cat('\n')
	# cat(OffPower); cat('\n')
	plot( c(OnTime, OffTime), c(OnPower, OffPower), xlab='MJD [sec]', ylab='Relative Power', type='n')
	points( OnTime, OnPower, pch=20, col='blue')
	points( OffTime, OffPower, pch=20, col='cyan')
	points( OnTime, basePower, pch=20, cex=0.5, col='gray' )
	return(baseTsys* (OnPower - basePower) / basePower)
}

#-------- On - Off subtraction for complex power
TxCal <- function(OnXpower, OffXpower, OffPower, Tsys, OnTime, OffTime, TsysTime ){
	baseXpower <- complex(
		real = predict(smooth.spline( OffTime, Re(OffXpower), spar=0.5 ), OnTime)$y,
		imaginary = predict(smooth.spline( OffTime, Im(OffXpower), spar=0.5 ), OnTime)$y )
	basePower <- predict(smooth.spline( OffTime, OffPower, spar=0.5 ), OnTime)$y
	baseTsys  <- predict(smooth.spline( TsysTime, Tsys, spar=0.5 ), OnTime)$y
	return(baseTsys* (OnXpower - baseXpower) / basePower)
}

#-------- Stokes Parameters
corr2Stokes <- function( XX, YY, XY, pang ){
	cs <- cos(2.0* pang);	sn <- sin(2.0* pang)	# pang: parallactic angle [rad]
	StokesI <- 0.5* (XX + YY)
	StokesQ <- 0.5* cs* (XX - YY) - sn* Re(XY)
	StokesU <- 0.5* sn* (XX - YY) + cs* Re(XY)
	StokesV <- Im(XY)
	return( data.frame(I=StokesI, Q=StokesQ, U=StokesU, V=StokesV) )
}


#-------- Load Spec and Scan data
args <- commandArgs(trailingOnly = T)
setwd('.')
#setwd('/Volumes/SSD/PolariS/20150317/')
#args <- c('2015076052841.Scan.Rdata', '2015076052912.SPEC.Rdata', '2015076035301.WG.Rdata', '2015076035301.BP.Rdata')
load(args[1])	 #Load Scan file
load(args[2])	 #Load SPEC file
load(args[3])	 #Load delay file
load(args[4])	 #Load BP file


#
#-------- Smoothed Delay and Phase
delay00Fit <- smooth.spline(WG$mjdSec, WG$delay00, spar=0.25); delay01Fit <- smooth.spline(WG$mjdSec, WG$delay01, spar=0.2)
Re00Fit <- smooth.spline(WG$mjdSec, Re(WG$Vis00), spar=0.25); Im00Fit <- smooth.spline(WG$mjdSec, Im(WG$Vis00), spar=0.2)
Re01Fit <- smooth.spline(WG$mjdSec, Re(WG$Vis01), spar=0.25); Im01Fit <- smooth.spline(WG$mjdSec, Im(WG$Vis01), spar=0.2)
#
#-------- Initial Parameters
chNum <- dim(on_A00)[1]
chRange <- floor(0.05*chNum):floor(0.95*chNum)

OnIndex <- which(Scan$scanType == 'ON')
OfIndex <- which(Scan$scanType == 'OFF')

#-------- Amplitude calibration of Autocorr
Ta00 <- TaCal( colSums(on_A00[chRange,]),  colSums(off_A00[chRange,]), Scan$Tsys00[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]) #; cat(Scan$Tsys00[OfIndex]); cat('\n') # ; cat(Ta00); cat('\n')
Ta01 <- TaCal( colSums(on_A01[chRange,]),  colSums(off_A01[chRange,]), Scan$Tsys01[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]) #; cat(Scan$Tsys01[OfIndex]); cat('\n') # ; cat(Ta01); cat('\n')
Ta02 <- TaCal( colSums(on_A02[chRange,]),  colSums(off_A02[chRange,]), Scan$Tsys02[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]) #; cat(Scan$Tsys02[OfIndex]); cat('\n') # ; cat(Ta02); cat('\n')
Ta03 <- TaCal( colSums(on_A03[chRange,]),  colSums(off_A03[chRange,]), Scan$Tsys03[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]) #; cat(Scan$Tsys03[OfIndex]); cat('\n') # ; cat(Ta03); cat('\n')

#-------- Delay, phase, bandpass, and amplitude calibration for CrossCorr
Tx00 <- TxCal( colSums( DelayPhaseCal( BPphsCal(on_C00, BP$BP00), scanTime(onMJD), delay00Fit, Re00Fit, Im00Fit)[chRange,]),
               colSums( DelayPhaseCal( BPphsCal(off_C00, BP$BP00), scanTime(offMJD), delay00Fit, Re00Fit, Im00Fit)[chRange,]),
			   sqrt( colSums(off_A00[chRange,])* colSums(off_A02[chRange,]) ),
			   sqrt(Scan$Tsys00[OfIndex]*  Scan$Tsys02[OfIndex]),
			   scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex])
Tx01 <- TxCal( colSums( DelayPhaseCal( BPphsCal(on_C01, BP$BP01), scanTime(onMJD), delay01Fit, Re01Fit, Im01Fit)[chRange,]),
               colSums( DelayPhaseCal( BPphsCal(off_C01, BP$BP01), scanTime(offMJD), delay01Fit, Re01Fit, Im01Fit)[chRange,]),
			   sqrt( colSums(off_A01[chRange,])* colSums(off_A03[chRange,]) ),
			   sqrt(Scan$Tsys01[OfIndex]*  Scan$Tsys03[OfIndex]), scanTime(onMJD),
			   scanTime(offMJD), Scan$mjdSec[OfIndex])
#-------- Parallactic angle
Pang <- -azel2pa(median(Scan$AZ[OnIndex]),  median(Scan$EL[OnIndex])) + median(Scan$EL[OnIndex])*pi/180 -pi/2	# Z45 receiver angle

#-------- Correlation -> Stokes Parameters
if( length( scanTime(onMJD) ) > 3){
	uncalStokes00 <- corr2Stokes( predict(smooth.spline(scanTime(onMJD), Ta00, spar=0.8), scanTime(onMJD))$y, predict(smooth.spline(scanTime(onMJD), Ta02, spar=0.8), scanTime(onMJD))$y, Tx00, Pang)
	uncalStokes01 <- corr2Stokes( predict(smooth.spline(scanTime(onMJD), Ta01, spar=0.8), scanTime(onMJD))$y, predict(smooth.spline(scanTime(onMJD), Ta03, spar=0.8), scanTime(onMJD))$y, Tx01, Pang)
} else {
	uncalStokes00 <- corr2Stokes(Ta00, Ta02, Tx00, Pang)
	uncalStokes01 <- corr2Stokes(Ta01, Ta03, Tx01, Pang)
}

uncalStokes00$p <- sqrt(uncalStokes00$Q^2 + uncalStokes00$U^2); uncalStokes00$EVPA <- atan2(uncalStokes00$U, uncalStokes00$Q)*90/pi; uncalStokes00$mjdSec <- scanTime(onMJD)
uncalStokes01$p <- sqrt(uncalStokes01$Q^2 + uncalStokes01$U^2); uncalStokes01$EVPA <- atan2(uncalStokes01$U, uncalStokes01$Q)*90/pi; uncalStokes01$mjdSec <- scanTime(onMJD)

cat(sprintf('AZ=%4.1f  EL=%4.1f  PA=%4.1f (deg)\n', median(Scan$AZ[OnIndex]), median(Scan$EL[OnIndex]), Pang*180/pi))
fileName <- sprintf("%s.uncalStokes.Rdata", onMJD[[1]][1])
save(uncalStokes00, uncalStokes01, file=fileName)
uncalStokes00
cat( mean(uncalStokes00$I), mean(uncalStokes00$Q), mean(uncalStokes00$U), mean(uncalStokes00$V), sqrt( mean(uncalStokes00$Q)^2 + mean(uncalStokes00$U)^2), atan2( mean(uncalStokes00$U), mean(uncalStokes00$Q))*90/pi); cat('\n')
cat( sd(uncalStokes00$I), sd(uncalStokes00$Q), sd(uncalStokes00$U), sd(uncalStokes00$V), sd(uncalStokes00$p), sd(uncalStokes00$EVPA)); cat('\n')
uncalStokes01
cat( mean(uncalStokes01$I), mean(uncalStokes01$Q), mean(uncalStokes01$U), mean(uncalStokes01$V), sqrt( mean(uncalStokes01$Q)^2 + mean(uncalStokes01$U)^2), atan2( mean(uncalStokes01$U), mean(uncalStokes01$Q))*90/pi); cat('\n')
cat( sd(uncalStokes01$I), sd(uncalStokes01$Q), sd(uncalStokes01$U), sd(uncalStokes01$V), sd(uncalStokes01$p), sd(uncalStokes01$EVPA)); cat('\n')
cat('Stokes parameters (D-term uncaled) was saved into '); cat(fileName); cat('\n')