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
TaCal <- function(OnPower, OffPower, Tsys, OnTime, OffTime, TsysTime ){
	basePower <- predict(smooth.spline( OffTime, OffPower, spar=0.8 ), OnTime)$y
	baseTsys  <- predict(smooth.spline( TsysTime, Tsys, spar=0.8 ), OnTime)$y
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
#args <- commandArgs()
#load(args[6])	 #Load Scan file
#load(args[7])	 #Load SPEC file
#load(args[8])	 #Load delay file
#load(args[9])	 #Load BP file
setwd('/Volumes/SSD/PolariS/20150317/')
load('2015076094523.Scan.Rdata')
load('2015076094555.SPEC.Rdata')
load('2015076035301.WG.Rdata')
load('2015076035301.BP.Rdata')
#load('2015076043056.Scan.Rdata')
#load('2015076043132.SPEC.Rdata')
#load('2015076035301.WG.Rdata')
#load('2015076035301.BP.Rdata')
#setwd('/Volumes/SSD/PolariS/20150316/')
#load('2015075034443.Scan.Rdata')
#load('2015075034519.SPEC.Rdata')
#load('2015075031846.XP.Rdata')
#load('2015075031846.WG.Rdata')
#load('2015075031846.BP.Rdata')
#setwd('/Volumes/SSD/PolariS/20150315/')
#load('2015074030852.BP.Rdata')
#load('2015074030852.WG.Rdata')
#load('2015074030852.XP.Rdata')
#load('2015074040644.Scan.Rdata')
#load('2015074040721.SPEC.Rdata')
#
#-------- Smoothed Delay and Phase
delay00Fit <- smooth.spline(WG$mjdSec, WG$delay00, spar=0.25); delay01Fit <- smooth.spline(WG$mjdSec, WG$delay01, spar=0.2)
Re00Fit <- smooth.spline(WG$mjdSec, Re(WG$Vis00), spar=0.25); Im00Fit <- smooth.spline(WG$mjdSec, Im(WG$Vis00), spar=0.2)
Re01Fit <- smooth.spline(WG$mjdSec, Re(WG$Vis01), spar=0.25); Im01Fit <- smooth.spline(WG$mjdSec, Im(WG$Vis01), spar=0.2)
#
#-------- Scan Pattern
mjdOn  <- floor(0.5*(onMJD$startMjd + onMJD$stopMjd))
mjdOff <- floor(0.5*(offMJD$startMjd + offMJD$stopMjd))

#-------- Initial Parameters
chNum <- dim(on_A00)[1]
chRange <- floor(0.05*chNum):floor(0.95*chNum)


#-------- Bandpass Calibration
BPphs00 <- complex( modulus=rep(1, length(BP$BP00)), argument=-Arg(BP$BP00))	# Phase-only bandpass table
BPphs01 <- complex( modulus=rep(1, length(BP$BP01)), argument=-Arg(BP$BP01))
on_C00 <- on_C00* BPphs00; off_C00 <- off_C00* BPphs00; R_C00 <- R_C00* BPphs00
on_C01 <- on_C01* BPphs01; off_C01 <- off_C01* BPphs01; R_C01 <- R_C01* BPphs01

#-------- Delay and Phase Calibration

delayCal_onC00 <- DelayPhaseCal( on_C00, mjdOn, delay00Fit, Re00Fit, Im00Fit )
delayCal_onC01 <- DelayPhaseCal( on_C01, mjdOn, delay01Fit, Re01Fit, Im01Fit )
delayCal_offC00 <- DelayPhaseCal( off_C00, mjdOff, delay00Fit, Re00Fit, Im00Fit )
delayCal_offC01 <- DelayPhaseCal( off_C01, mjdOff, delay01Fit, Re01Fit, Im01Fit )

if(0){
onTsys0 <- Scan$Tsys00[OnIndex]
onTsys1 <- Scan$Tsys01[OnIndex]
onTsys2 <- Scan$Tsys02[OnIndex]
onTsys3 <- Scan$Tsys03[OnIndex]
onTsys00 <- sqrt(onTsys0* onTsys2)
onTsys01 <- sqrt(onTsys1* onTsys3)
onweight0 <- 1.0/onTsys0^2
onweight1 <- 1.0/onTsys1^2
onweight2 <- 1.0/onTsys2^2
onweight3 <- 1.0/onTsys3^2
onweight00 <- sqrt(onweight0* onweight2)
onweight01 <- sqrt(onweight1* onweight3)

offTsys0 <- Scan$Tsys00[OfIndex]
offTsys1 <- Scan$Tsys01[OfIndex]
offTsys2 <- Scan$Tsys02[OfIndex]
offTsys3 <- Scan$Tsys03[OfIndex]
offTsys00 <- sqrt(offTsys0* offTsys2)
offTsys01 <- sqrt(offTsys1* offTsys3)
offweight0 <- 1.0/offTsys0^2
offweight1 <- 1.0/offTsys1^2
offweight2 <- 1.0/offTsys2^2
offweight3 <- 1.0/offTsys3^2
offweight00 <- sqrt(offweight0* offweight2)
offweight01 <- sqrt(offweight1* offweight3)

#-------- Smoothed Bandpass Calibration
mitigCH <- c(7256, 16385, 32769, 47522)
flagCH <- c(1, 2, 4, mitigCH)
weight <- rep(1, chNum); weight[flagCH] <- 0.0
availCH <- which(weight == 1.0)
smoothWidth <- 512; knotNum <- floor(chNum / smoothWidth)
freq <- freq <- (0:(chNum-1))/chNum* 4.0	# MHz
chSep <- 4.0 / chNum

TA00 <- on_A00; TA01 <- on_A01; TA02 <- on_A02; TA03 <- on_A03; TX00 <- delayCal_onC00; TX01 <- delayCal_onC01;
for(timeIndex in 1:length(mjdOn)){
	tryCatch(
		{ off_before <- max(which( mjdOff < mjdOn[timeIndex] ))},
		warning = function(e){ off_before = which.min(mjdOff)},
		silent = TRUE )
	tryCatch(
		{ off_after  <- min(which( mjdOff > mjdOn[timeIndex] )) },
		warning = function(e){ off_after = which.max(mjdOff)},
		silent = TRUE )
	#
	#-------- Power Spectra
	smoothOffA00 <- SBCspec( 0.5*(off_A00[,off_before] + off_A00[,off_after]), knotNum, weight, mitigCH)
	TA00[,timeIndex] <- onTsys0[timeIndex]* (SPmitigation(on_A00[,timeIndex], mitigCH) - smoothOffA00)/smoothOffA00
	smoothOffA01 <- SBCspec( 0.5*(off_A01[,off_before] + off_A01[,off_after]), knotNum, weight, mitigCH)
	TA01[,timeIndex] <- onTsys1[timeIndex]* (SPmitigation(on_A01[,timeIndex], mitigCH) - smoothOffA01)/smoothOffA01
	smoothOffA02 <- SBCspec( 0.5*(off_A02[,off_before] + off_A02[,off_after]), knotNum, weight, mitigCH)
	TA02[,timeIndex] <- onTsys2[timeIndex]* (SPmitigation(on_A02[,timeIndex], mitigCH) - smoothOffA02)/smoothOffA02
	smoothOffA03 <- SBCspec( 0.5*(off_A03[,off_before] + off_A03[,off_after]), knotNum, weight, mitigCH)
	TA03[,timeIndex] <- onTsys3[timeIndex]* (SPmitigation(on_A03[,timeIndex], mitigCH) - smoothOffA03)/smoothOffA03
	#-------- Cross Power Spectra
	smoothOffC00 <- complex(
		real = SBCspec(0.5*Re(delayCal_offC00[,off_before] + delayCal_offC00[,off_after]), knotNum, weight, mitigCH),
		imaginary = SBCspec(0.5*Im(delayCal_offC00[,off_before] + delayCal_offC00[,off_after]), knotNum, weight, mitigCH))
	smoothOffC01 <- complex(
		real = SBCspec(0.5*Re(delayCal_offC01[,off_before] + delayCal_offC01[,off_after]), knotNum, weight, mitigCH),
		imaginary = SBCspec(0.5*Im(delayCal_offC01[,off_before] + delayCal_offC01[,off_after]), knotNum, weight, mitigCH))
	TX00[,timeIndex] <- onTsys00[timeIndex]* (delayCal_onC00[,timeIndex] - smoothOffC00)/ sqrt(smoothOffA00* smoothOffA02)
	TX01[,timeIndex] <- onTsys01[timeIndex]* (delayCal_onC01[,timeIndex] - smoothOffC01)/ sqrt(smoothOffA01* smoothOffA03)
}

#-------- Time Integration
StokesI02 <- 0.5*( apply( TA00, 1, mean ) + apply( TA02, 1, mean )) 
StokesI13 <- 0.5*( apply( TA01, 1, mean ) + apply( TA03, 1, mean )) 
integTX00 <- apply( TX00, 1, mean)
integTX01 <- apply( TX01, 1, mean)
obsStokesV02 <- Im(integTX00)
obsStokesV13 <- Im(integTX01)

#-------- Plot Stokes I for HC3N
#pdf('TMC-1_CCS_Stokes.pdf')
plotRange <- 25001:32768
fitStokesI <- smooth.spline(freq, StokesI02, w=weight, all.knots=F, nknots=4*knotNum)
predStokesV <- predict(fitStokesI, (freq + 5e-6))$y - predict(fitStokesI, (freq - 5e-6))$y
plotBunch <- 16
plot( bunch_vec(freq[plotRange],plotBunch)-0.5*plotBunch*chSep, bunch_vec(StokesI02[plotRange], plotBunch), type='s', xlab='Frequency [MHz]', ylab='Stokes I [K]', main='TMC-1 HC3N')
lines(freq[plotRange], predict(fitStokesI, freq[plotRange])$y, type='s', col='red')
#-------- Plot Stokes V for HC3N
lineRange <- 25001:32768
plotBunch <- 32
fit <- lm( formula=y~1+ x, data=data.frame(x=predStokesV[lineRange], y=obsStokesV02[lineRange]))
plot( bunch_vec(freq[plotRange], plotBunch) - 0.5*plotBunch*chSep, bunch_vec(obsStokesV02[plotRange], plotBunch) - fit$coefficients[1], type='s', ylim=c(-1e-1, 1e-1), xlab='Frequency [MHz]', ylab='Stokes V [K]', main='TMC-1 HC3N')
points( bunch_vec(freq[plotRange], plotBunch), bunch_vec(obsStokesV02[plotRange], plotBunch) - fit$coefficients[1], pch=20, cex=0.5)
lines(freq[plotRange]-0.5*chSep, predStokesV[plotRange]*fit$coefficients[2], type='s', col='red')
legend(min(freq[plotRange]), 0.1, legend=sprintf('Zeeman Shift = %5.1f ± %4.1f Hz', 10*fit[[1]][2], 10*summary(fit)[[4]][4]))


#-------- Plot Stokes I for CCS
lineRange <- 25001:32768
fitStokesI <- smooth.spline(freq, StokesI13, w=weight, all.knots=F, nknots=4*knotNum)
predStokesV <- predict(fitStokesI, (freq + 5e-6))$y - predict(fitStokesI, (freq - 5e-6))$y
plotBunch <- 16
plot( bunch_vec(freq[plotRange],plotBunch)-0.5*plotBunch*chSep, bunch_vec(StokesI13[plotRange], plotBunch), type='s', xlab='Frequency [MHz]', ylab='Stokes I [K]', main='TMC-1 CCS')
lines(freq[plotRange], predict(fitStokesI, freq[plotRange])$y, type='s', col='red')
#-------- Plot Stokes V for CCS
lineRange <- 31131:35226
plotBunch <- 32
fit <- lm( formula=y~1+ x, data=data.frame(x=predStokesV[lineRange], y=obsStokesV13[lineRange]))
plot( bunch_vec(freq[plotRange], plotBunch) - 0.5*plotBunch*chSep, bunch_vec(obsStokesV13[plotRange], plotBunch) - fit$coefficients[1], type='s', ylim=c(-1e-1, 1e-1), xlab='Frequency [MHz]', ylab='Stokes V [K]', main='TMC-1 CCS')
points( bunch_vec(freq[plotRange], plotBunch), bunch_vec(obsStokesV13[plotRange], plotBunch) - fit$coefficients[1], pch=20, cex=0.5)
lines(freq[plotRange]-0.5*chSep, predStokesV[plotRange]*fit$coefficients[2], type='s', col='red')
legend(min(freq[plotRange]), 0.1, legend=sprintf('Zeeman Shift = %5.1f ± %4.1f Hz', 10*fit[[1]][2], 10*summary(fit)[[4]][4]))

}
