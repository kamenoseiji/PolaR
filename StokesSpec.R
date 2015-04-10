# StokeSpec
# usage: Rscript StokesSpec [Scan.Rdata file name] [SPEC.Rdata file name] [WG.Rdata file name] [BP file name]
#
library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/Qeff.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readSAM45.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/mjd.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/delaysearch.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readPolariS.R", ssl.verifypeer = FALSE)))
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

#-------- Spectral Mitigation
SPmitigation <- function( spec, SPCH ){
	for(index in 1:length(SPCH)){
		spec[SPCH[index]] <- mean(spec[c((SPCH[index]-8):(SPCH[index]-1), (SPCH[index]+1):(SPCH[index]+8))])
	}
	return(spec)
}
#-------- Smoothed Bandpass Calibration
SBCspec <- function( spec, knotNum, weight, mitigCH){
	return(predict(smooth.spline( SPmitigation(spec, mitigCH), w=weight, all.knots=F, nknots=knotNum), 1:length(spec))$y)
}

#-------- On - Off subtraction for total power
TaCal <- function(OnSpec, OffSpec, mjdOn, mjdOff, Tsys, weight, mitigCH){
	#-------- Parameters of Smoothed Bandpass Calibration
	chNum <- nrow( OnSpec )
	#flagCH <- c(1, 2, 4, mitigCH)
	#availCH <- which(weight == 1.0)
	smoothWidth <- 512; knotNum <- floor(chNum / smoothWidth)

	TA <- matrix( nrow=chNum, ncol=length(mjdOn) )
	for( timeIndex in 1:length(mjdOn) ){
		#-------- if scan starts on-source
		tryCatch(
			{ off_before <- max(which( mjdOff < mjdOn[timeIndex] ))},
			warning = function(e){ off_before = which.min(mjdOff)},
			silent = TRUE )
		#-------- if scan ends on-source
		tryCatch(
			{ off_after  <- min(which( mjdOff > mjdOn[timeIndex] )) },
			warning = function(e){ off_after = which.max(mjdOff)},
			silent = TRUE )
		#
		#---- Off-source power spectra
		SmoothOffSpec <- SBCspec( 0.5*(OffSpec[,off_before] + OffSpec[,off_after]), knotNum, weight, mitigCH)
		TA[,timeIndex] <- Tsys[timeIndex]* (SPmitigation(OnSpec[,timeIndex], mitigCH) - SmoothOffSpec)/SmoothOffSpec
	}
	return(TA)
}

#-------- On - Off subtraction for complex power
TxCal <- function(OnXspec, OffXspec, OffSpec0, OffSpec1, mjdOn, mjdOff, Tsys, weight, mitigCH ){
	# Cross-power spectra are assumed to be calibrated in terms of Delay and Phase
	chNum <- nrow( OnXspec )
	mitigCH <- c(7256, 16385, 32769, 47522)
	flagCH <- c(1, 2, 4, mitigCH)
	weight <- rep(1, chNum); weight[flagCH] <- 0.0
	availCH <- which(weight == 1.0)
	smoothWidth <- 512; knotNum <- floor(chNum / smoothWidth)
	TX <- matrix( nrow=chNum, ncol=length(mjdOn) )
	for( timeIndex in 1:length(mjdOn) ){
		#-------- if scan starts on-source
		tryCatch(
			{ off_before <- max(which( mjdOff < mjdOn[timeIndex] ))},
			warning = function(e){ off_before = which.min(mjdOff)},
			silent = TRUE )
		#-------- if scan ends on-source
		tryCatch(
			{ off_after  <- min(which( mjdOff > mjdOn[timeIndex] )) },
			warning = function(e){ off_after = which.max(mjdOff)},
			silent = TRUE )
		#
		#---- Off-source power spectra
		SmoothOffSpec <- sqrt( SBCspec( 0.5*(OffSpec0[,off_before] + OffSpec0[,off_after]), knotNum, weight, mitigCH)* SBCspec( 0.5*(OffSpec1[,off_before] + OffSpec1[,off_after]), knotNum, weight, mitigCH) )
		
		#---- Off-source cross-power spectra
		SmoothOffXspec <- complex(
			real = SBCspec(0.5*Re(OffXspec[,off_before] + OffXspec[,off_after]), knotNum, weight, mitigCH),
			imaginary = SBCspec(0.5*Im(OffXspec[,off_before] + OffXspec[,off_after]), knotNum, weight, mitigCH))
			
		#---- Scaling
		TX[, timeIndex] <- Tsys[timeIndex]* (OnXspec[,timeIndex] - SmoothOffXspec) / SmoothOffSpec
	}
	return(TX)
}
#-------- Doppler Tracking
deDopp <- function( spec, chShift ){
	timeNum <- ncol(spec)
	chNum <- nrow(spec)
	for( time_index in 1:timeNum ){
		if(chShift[time_index] > 0){
			spec[,time_index] <- c(spec[(chNum - chShift[time_index] + 1):chNum, time_index], spec[1:(chNum - chShift[time_index]), time_index])
		}
		if(chShift[time_index] < 0){
			spec[,time_index] <- c(spec[(chShift[time_index] + 1):chNum, time_index], spec[1:chShift[time_index], time_index])
		}
	}
	return(spec)
}

#-------- Load Spec and Scan data
args <- commandArgs(trailingOnly = T)
setwd('.')
#setwd('/Volumes/SSD/PolariS/20140417/')
#args <- c('2014107013853.Scan.Rdata', '2014107013932.SPEC.Rdata', '2014107040319.Dcomb.Rdata', '2014107010610.WG.Rdata', '2014107010610.BP.Rdata')
#args <- c('2015076052841.Scan.Rdata', '2015076052912.SPEC.Rdata', '2015076035301.WG.Rdata', '2015076035301.BP.Rdata')
#args <- c('2015076043056.Scan.Rdata', '2015076043132.SPEC.Rdata', '2015076035301.WG.Rdata', '2015076035301.BP.Rdata')
load(args[1])	 #Load Scan file
load(args[2])	 #Load SPEC file
load(args[3])	 #Load D-term file
load(args[4])	 #Load delay file
load(args[5])	 #Load BP file
#
#-------- Smoothed Delay and Phase
delay00Fit <- smooth.spline(WG$mjdSec, WG$delay00, spar=0.25)
delay01Fit <- smooth.spline(WG$mjdSec, WG$delay01, spar=0.25)
Re00Fit <- smooth.spline(WG$mjdSec, Re(WG$Vis00), spar=0.25)
Im00Fit <- smooth.spline(WG$mjdSec, Im(WG$Vis00), spar=0.25)
Re01Fit <- smooth.spline(WG$mjdSec, Re(WG$Vis01), spar=0.25)
Im01Fit <- smooth.spline(WG$mjdSec, Im(WG$Vis01), spar=0.25)
#
#-------- Initial Parameters
chNum <- dim(on_A00)[1]
chRange <- floor(0.05*chNum):floor(0.95*chNum)
freq <- (0:(chNum-1))/chNum* 4.0	# MHz
chSep <- 4.0 / chNum
#mitigCH <- c(7256, 16385, 32769, 47522)
mitigCH <- c(16385)
flagCH <- c(1, 2, 4, mitigCH)
weight <- rep(1, chNum); weight[flagCH] <- 0.0
# smoothWidth <- 512; knotNum <- floor(chNum / smoothWidth)
smoothWidth <- 384; knotNum <- floor(chNum / smoothWidth)
#-------- Scan Pattern
OnIndex <- which(Scan$scanType == 'ON')
OfIndex <- which(Scan$scanType == 'OFF')
D_index <- which( D$Gxy02 + D$Gxy13 > 0.9* median( D$Gxy02 + D$Gxy13 ))	# Flag pointing-error out
Rxy02 <- mean(D$Rxy02[D_index])
Rxy13 <- mean(D$Rxy13[D_index])
Dxy02 <- mean(D$XY02[D_index])
Dxy13 <- mean(D$XY13[D_index])
#-------- Tsys at each scan
Tsys00 <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$Tsys00[OnIndex], spar=0.5), scanTime(onMJD))$y
Tsys01 <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$Tsys01[OnIndex], spar=0.5), scanTime(onMJD))$y
Tsys02 <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$Tsys02[OnIndex], spar=0.5), scanTime(onMJD))$y
Tsys03 <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$Tsys03[OnIndex], spar=0.5), scanTime(onMJD))$y

#-------- Radial Velocity
#DeltaFreq <- 317.00 - 206.00	# [MHz], LO(CH1) - LO(CH2)	// 2014.4.17
DeltaFreq <- 317.50 - 206.20	# [MHz], LO(CH1) - LO(CH2)  // 2015.3.15 - 17
Vrad <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$Vrad[OnIndex], spar=0.5), scanTime(onMJD))$y
chShift <- round(DeltaFreq* Vrad / 299792458 / chSep)		# Number of channels to shift

#-------- Parallactic angle
AZ <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$AZ[OnIndex], spar=0.5), scanTime(onMJD))$y
EL <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$EL[OnIndex], spar=0.5), scanTime(onMJD))$y
Pang <- -azel2pa(AZ, EL) + EL*pi/180 - pi/2 
cs <- cos(Pang)
sn <- sin(Pang)
#-------- Amplitude calibration of Autocorr
Ta00 <- TaCal(on_A00, off_A00, scanTime(onMJD), scanTime(offMJD), Tsys00, weight, mitigCH) / sqrt(Rxy02)
Ta01 <- TaCal(on_A01, off_A01, scanTime(onMJD), scanTime(offMJD), Tsys01, weight, mitigCH) / sqrt(Rxy13)
Ta02 <- TaCal(on_A02, off_A02, scanTime(onMJD), scanTime(offMJD), Tsys02, weight, mitigCH) * sqrt(Rxy02)
Ta03 <- TaCal(on_A03, off_A03, scanTime(onMJD), scanTime(offMJD), Tsys03, weight, mitigCH) * sqrt(Rxy13)

#-------- Delay, phase, bandpass, and amplitude calibration for CrossCorr
Tx02 <- TxCal( 
	DelayPhaseCal( BPphsCal(on_C00, BP$BP00), scanTime(onMJD), delay00Fit, Re00Fit, Im00Fit),
	DelayPhaseCal( BPphsCal(off_C00, BP$BP00), scanTime(offMJD), delay00Fit, Re00Fit, Im00Fit),
	off_A00, off_A02, scanTime(onMJD), scanTime(offMJD), sqrt(Tsys00 * Tsys02), weight, mitigCH)

Tx13 <- TxCal( 
	DelayPhaseCal( BPphsCal(on_C01, BP$BP01), scanTime(onMJD), delay01Fit, Re01Fit, Im01Fit),
	DelayPhaseCal( BPphsCal(off_C01, BP$BP01), scanTime(offMJD), delay01Fit, Re01Fit, Im01Fit),
	off_A01, off_A03, scanTime(onMJD), scanTime(offMJD), sqrt(Tsys01 * Tsys03), weight, mitigCH)
#-------- Doppler Tracking
Ta00 <- deDopp(Ta00, chShift)
Ta02 <- deDopp(Ta02, chShift)
Tx02 <- deDopp(Tx02, chShift)
#-------- Time Integration to determine Stokes spectrum
weight0 <- 1.0 / Tsys00^2
weight1 <- 1.0 / Tsys01^2
weight2 <- 1.0 / Tsys02^2
weight3 <- 1.0 / Tsys03^2
weight02 <- 1.0 / (Tsys00 * Tsys02)
weight13 <- 1.0 / (Tsys01 * Tsys03)
StokesI02 <- 0.5*(rowSums(Ta00* weight0)/sum(weight0) + rowSums(Ta02* weight2)/sum(weight2))
StokesI13 <- 0.5*(rowSums(Ta01* weight1)/sum(weight1) + rowSums(Ta03* weight3)/sum(weight3))

StokesQ02 <- 0.5* (rowSums(Ta00* cs* weight0)/sum(weight0) - rowSums(Ta02* cs* weight2)/sum(weight2)) - rowSums((Re(Tx02) - Re(Dxy02)* StokesI02)* sn* weight02) / sum(weight02)
StokesQ13 <- 0.5* (rowSums(Ta01* cs* weight1)/sum(weight1) - rowSums(Ta03* cs* weight3)/sum(weight3)) - rowSums((Re(Tx13) - Re(Dxy13)* StokesI13)* sn* weight13) / sum(weight13)

StokesU02 <- 0.5* (rowSums(Ta00* sn* weight0)/sum(weight0) - rowSums(Ta02* sn* weight2)/sum(weight2)) + rowSums((Re(Tx02) - Re(Dxy02)* StokesI02)* cs* weight02) / sum(weight02)
StokesU13 <- 0.5* (rowSums(Ta01* sn* weight1)/sum(weight1) - rowSums(Ta03* sn* weight3)/sum(weight3)) + rowSums((Re(Tx13) - Re(Dxy13)* StokesI13)* cs* weight13) / sum(weight13)

StokesV02 <- rowSums(Im(Tx02)* weight02)/sum(weight02) - Im(Dxy02)* StokesI02
StokesV13 <- rowSums(Im(Tx13)* weight13)/sum(weight13) - Im(Dxy13)* StokesI13

#-------- Save into file
fileName <- sprintf("%s.STOKES.Rdata", strsplit(args[2], "\\.")[[1]][1])
save(StokesI02, StokesI13, StokesQ02, StokesQ13, StokesU02, StokesU13, StokesV02, StokesV13, file=fileName)
cat('Stokes spectra are saved into '); cat(fileName); cat('\n')
