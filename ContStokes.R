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
	plot( c(OnTime, OffTime), c(OnPower, OffPower), xlab='MJD [sec]', ylab='Relative Power', type='n')
	points( OnTime, OnPower, pch=20, col='blue')
	points( OffTime, OffPower, pch=20, col='cyan')
	points( OnTime, basePower, pch=20, cex=0.5, col='gray' )
	return(baseTsys* (OnPower - basePower) / basePower)
}
#TaCal <- function(OnPower, OffPower, Tsys, OnTime, OffTime, TsysTime ){
	#basePower <- predict(smooth.spline( OffTime, OffPower, spar=0.0 ), OnTime)$y
	#baseTsys  <- predict(smooth.spline( TsysTime, Tsys, spar=0.0 ), OnTime)$y
	#plot( c(OnTime, OffTime), c(OnPower, OffPower), xlab='MJD [sec]', ylab='Relative Power', type='n')
	#points( OnTime, OnPower, pch=20, col='blue')
	#points( OffTime, OffPower, pch=20, col='cyan')
	#lines( OnTime, basePower )
	#return(baseTsys* (OnPower - basePower) / basePower)
#}

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
	StokesI <- 0.5* (XX                             + YY)
	StokesQ <- 0.5* (cs*XX - sn*XY - sn*Conj(XY) - cs*YY)
	StokesU <- 0.5* (sn*XX + cs*XY + cs*Conj(XY) - sn*YY)
	StokesV <- 0.5* (    (0-1i)*XY + (0+1i)*Conj(XY)    )
	return( data.frame(I=StokesI, Q=Re(StokesQ), U=Re(StokesU), V=Re(StokesV)) )
}


#-------- Load Spec and Scan data
args <- commandArgs()
load(args[6])	 #Load Scan file
load(args[7])	 #Load SPEC file
load(args[8])	 #Load delay file
load(args[9])	 #Load BP file
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
R_Index <- which(Scan$scanType == 'R')

Ta00 <- TaCal( Scan$power00[OnIndex],  Scan$power00[OfIndex], Scan$Tsys00[OfIndex], Scan$mjdSec[OnIndex], Scan$mjdSec[OfIndex], Scan$mjdSec[OfIndex]); cat(Ta00); cat('\n')
Ta01 <- TaCal( Scan$power01[OnIndex],  Scan$power01[OfIndex], Scan$Tsys01[OfIndex], Scan$mjdSec[OnIndex], Scan$mjdSec[OfIndex], Scan$mjdSec[OfIndex]); cat(Ta01); cat('\n')
Ta02 <- TaCal( Scan$power02[OnIndex],  Scan$power02[OfIndex], Scan$Tsys02[OfIndex], Scan$mjdSec[OnIndex], Scan$mjdSec[OfIndex], Scan$mjdSec[OfIndex]); cat(Ta02); cat('\n')
Ta03 <- TaCal( Scan$power03[OnIndex],  Scan$power03[OfIndex], Scan$Tsys03[OfIndex], Scan$mjdSec[OnIndex], Scan$mjdSec[OfIndex], Scan$mjdSec[OfIndex]); cat(Ta03); cat('\n')

#Ta00 <- TaCal( colSums(on_A00[chRange,]),  colSums(off_A00[chRange,]), Scan$Tsys00[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]); cat(Ta00); cat('\n')
#Ta01 <- TaCal( colSums(on_A01[chRange,]),  colSums(off_A01[chRange,]), Scan$Tsys01[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]); cat(Ta01); cat('\n')
#Ta02 <- TaCal( colSums(on_A02[chRange,]),  colSums(off_A02[chRange,]), Scan$Tsys02[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]); cat(Ta02); cat('\n')
#Ta03 <- TaCal( colSums(on_A03[chRange,]),  colSums(off_A03[chRange,]), Scan$Tsys03[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]); cat(Ta03); cat('\n')

Tx00 <- TxCal( colSums(DelayPhaseCal(BPphsCal(on_C00, BP$BP00), scanTime(onMJD), delay00Fit, Re00Fit, Im00Fit)[chRange,]),  colSums(DelayPhaseCal(BPphsCal(off_C00, BP$BP00), scanTime(offMJD), delay00Fit, Re00Fit, Im00Fit)[chRange,]), colSums(off_A00[chRange,]), Scan$Tsys00[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]); cat(Tx00); cat('\n')
Tx01 <- TxCal( colSums(DelayPhaseCal(BPphsCal(on_C01, BP$BP01), scanTime(onMJD), delay01Fit, Re01Fit, Im01Fit)[chRange,]),  colSums(DelayPhaseCal(BPphsCal(off_C01, BP$BP01), scanTime(offMJD), delay01Fit, Re01Fit, Im01Fit)[chRange,]), colSums(off_A01[chRange,]), Scan$Tsys01[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]); cat(Tx01); cat('\n')

Pang <- -azel2pa(median(Scan$AZ[OnIndex]),  median(Scan$EL[OnIndex])) + median(Scan$EL[OnIndex])*pi/180 -pi/2	# Z45 receiver angle
uncalStokes00 <- corr2Stokes( Ta00, Ta02, Tx00, Pang); uncalStokes00$p <- sqrt(uncalStokes00$Q^2 + uncalStokes00$U^2); uncalStokes00$EVPA <- atan2(uncalStokes00$U, uncalStokes00$Q)*90/pi; uncalStokes00$mjdSec <- scanTime(onMJD)
uncalStokes01 <- corr2Stokes( Ta01, Ta03, Tx01, Pang); uncalStokes01$p <- sqrt(uncalStokes01$Q^2 + uncalStokes01$U^2); uncalStokes01$EVPA <- atan2(uncalStokes01$U, uncalStokes01$Q)*90/pi; uncalStokes01$mjdSec <- scanTime(onMJD)
cat(sprintf('AZ=%4.1f  EL=%4.1f  PA=%4.1f (deg)\n', median(Scan$AZ[OnIndex]), median(Scan$EL[OnIndex]), Pang*180/pi))
fileName <- sprintf("%s.uncalStokes.Rdata", onMJD[[1]][1])
save(uncalStokes00, uncalStokes01, file=fileName)
uncalStokes00	
uncalStokes01
cat('Stokes parameters (D-term uncaled) was saved into '); cat(fileName); cat('\n')