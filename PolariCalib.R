#-------- Bunching given vector
bunch_vec <- function(vector, bunchNum){	# Function to bunch 1-d vector
	pudding <- (bunchNum - length(vector) %% bunchNum) %% bunchNum
	temp <- matrix( append(vector, rep(vector[length(vector)], pudding)), nrow=bunchNum)
	return( apply(temp, 2, mean))
}

#-------- Median Window Filter (mwf)
medianWindowFilter <- function(vector, width){
    if(width %% 2 == 0){ width <- width + 1 }
    halfWidth <- (width + 1)/2
    vecLen <- length(vector)
    beginWindow <- c( rep(1, halfWidth), seq(2, vecLen - width + 1), rep(vecLen - width + 1, halfWidth - 1))
    endWindow   <- beginWindow + width - 1
    mwf <- numeric(vecLen)
    for(index in 1:vecLen){
        mwf[index] <- median(vector[beginWindow[index]:endWindow[index]])
    }
    return(mwf)
}
#-------------------------------------- Function to calculate cross spectrum from correlation function
corr2spec <- function( corr ){
	nspec <- length(corr)/2
	spec <- fft(corr)[1:nspec]
	return(spec)
}
#-------------------------------------- Function to calculate correlation function from cross spectrum
spec2corr <- function(spec){
	nspec <- length(spec)
	tempcorr <- fft( c(spec, rep(0.0, nspec)), inverse=T)
	corr <- c(tail(tempcorr, nspec), head(tempcorr, nspec))
	return(corr)
	# tmpspec <- c(spec, 0, Conj(spec[nspec:2]))
	# corr <- Re(fft(tmpspec, inverse=TRUE) / nspec)
	# return(c(corr[(nspec+1):(2*nspec)], corr[1:nspec]))
}

#-------------------------------------- Delay Search
delay_search <- function( spec, weight=0 ){
	nspec <- length( spec )
    if(length(weight) < length(spec)){    weight <- rep(1.0, nspec)}
    spec <- spec * weight
	#-------- Search for delay
	delay <- which.max(Mod(spec2corr(spec))) - nspec - 1	# Coarse Delay
	trial_delay <- delay + seq(-1, 1, by=0.1); trial_amp <- numeric(0)
	for(i in 1:length(trial_delay)){ trial_amp[i] <- Mod(sum(delay_cal(spec, trial_delay[i]))) }
	fit <- lm( formula = trial_amp ~ trial_delay + I(trial_delay^2))
	return(-fit[[1]][[2]] / (2.0* fit[[1]][[3]]))
}

#-------------------------------------- Delay calibration
delay_cal <- function( spec, delay ){
	# spec : input spectrum (complex)
	# delay_cal() returns delay-calibrated spectrum
	#
	nspec <- length( spec )
	twiddle <- complex(modulus=rep(1, nspec), argument = delay* pi* seq((-nspec/2), (nspec/2 - 1), by=1) / nspec)
	return( spec * twiddle )
}

#-------------------------------------- Delay and Phase calibration
delayPhase_cal <- function( spec, delay, phase ){
	nspec <- length( spec )
	twiddle <- complex(modulus=rep(1, nspec), argument = phase + delay* pi* seq((-nspec/2), (nspec/2 - 1), by=1) / nspec)
	return( spec* twiddle )
}

#-------------------------------------- Complex Bandpass Table
BPtable <- function( spec, smoothCH ){
    mwfBP <- medianWindowFilter(Re(spec), smoothCH) + (0+1i)* medianWindowFilter(Im(spec), smoothCH)
    diffSpec <- Mod(spec - mwfBP)
    thresh <- 5.0* median(diffSpec)
    mask <- which(diffSpec < thresh)
    weight <- rep(0.0, length(spec)); weight[mask] <- 1.0
	delay <- delay_search(spec, weight)
	delCaledSpec <- delay_cal(spec, delay)
	return(smoothComplex(delCaledSpec/mean(delCaledSpec[mask]), smoothCH))
}

#-------------------------------------- Smoothing complex vector
smoothComplex <- function( spec, smoothCH ){
	nspec <- length(spec)
	reSpec <- predict(smooth.spline( Re(spec), all.knots=F, nknots=as.integer(nspec/smoothCH)), 1:length(spec))$y
	imSpec <- predict(smooth.spline( Im(spec), all.knots=F, nknots=as.integer(nspec/smoothCH)), 1:length(spec))$y
	return(complex(real=reSpec, imaginary=imSpec))
}
#-------- Delay and phase determination using Wire-Grid scans
wireGridPhaseCorr <- function(XY, ScanIndex, chRange=9:65544){
	# XY is the cross power spectrum
	# ScanIndex is the time index for WG observation (e.g. 1234:1421)
	# chRange is spectral channels range to be used (e.g. 9:65544)
	Xspec <- apply(XY[chRange, ScanIndex], 1, mean)	# CrossCorr
	delay <- delay_search(Xspec)
	return( list(delay=delay, phase=Arg(mean(delay_cal(Xspec, delay)))))
}

#-------- Determine Dx+Dy* using unpolarized source
unpolDterm <- function( XX, YY, XY, scan, chRange, Tsys, delay, phase ){
	#-------- On- and Off-source power
	XX_on <- mean(XX[chRange, scan$on]); XX_off <- mean(XX[chRange, scan$off])			# Time and Spectral Averaging
	YY_on <- mean(YY[chRange, scan$on]); YY_off <- mean(YY[chRange, scan$off])			# Time and Spectral Averaging
	XY_on <- mean(delay_cal(apply(XY[chRange, scan$on], 1, mean), delay))
	XY_off <- mean(delay_cal(apply(XY[chRange, scan$off], 1, mean), delay))
	#-------- Gain = polaris unit / K
	GainXX <- XX_off / Tsys[1]		# Gx Gx* = <XX*> / <NxNx*>, unit:[polaris/K]
	GainYY <- YY_off / Tsys[2]		# Gy Gy* = <YY*> / <NyNy*>, unit:[polaris/K]
	Gx <- complex(modulus=sqrt(GainXX), argument=0)			# Complex Gain
	Gy <- complex(modulus=sqrt(GainYY), argument=phase)	# Complex Gain
	StokesI <- 0.5*((XX_on - XX_off)/GainXX + (YY_on - YY_off)/GainYY) 	# Antenna temperature of the unpolarized source [K]
	Dx_plus_Dy <- (XY_on - XY_off) / (Gx * Gy* StokesI)
	return(list(StokesI, Dx_plus_Dy))
}

#-------- BP Cal
BPphsCal <- function( SPEC, BP ){
	BPphs <- complex( modulus=1, argument=-Arg(BP))
	return( SPEC* BPphs )
}

#-------- Function to Calibrate Delay and Phase
DelayPhaseCal <- function( scanSpec, mjdSec, delayFit, ReFit, ImFit){
	temp <- scanSpec
    if(summary(ImFit)[[1]][1] == 'lm()'){
	    phase <- atan2(as.numeric(predict(ImFit, data.frame(x=mjdSec))), as.numeric(predict(ReFit, data.frame(x=mjdSec))))
    } else {
	    phase <- atan2(predict(ImFit, x=mjdSec)$y, predict(ReFit, x=mjdSec)$y)
    }
    if(summary(delayFit)[[1]][1] == 'lm()'){
	    delay <- as.numeric(predict(delayFit, data.frame(x=mjdSec)))
    } else {
	    delay <- predict(delayFit, x=mjdSec)$y
    }
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
	return(baseTsys* (OnPower - basePower) / basePower)
}

#-------- On - Off subtraction for Spectral Data
TaCalSpec <- function(OnSpec, OffSpec, mjdOn, mjdOff, Tsys, weight, mitigCH){
    #-------- Parameters of Smoothed Bandpass Calibration
    chNum <- nrow( OnSpec )
    knotNum <- 120
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
TxCal <- function(OnXpower, OffXpower, OffPower, Tsys, OnTime, OffTime, TsysTime ){
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
#-------- On - Off subtraction for complex power
TxCalSpec <- function(OnXspec, OffXspec, OffSpec0, OffSpec1, mjdOn, mjdOff, Tsys, weight, mitigCH ){
    # Cross-power spectra are assumed to be calibrated in terms of Delay and Phase
    chNum <- nrow( OnXspec )
    # mitigCH <- c(7256, 16385, 32769, 47522)
    #flagCH <- c(1, 2, 4, mitigCH)
    # weight <- rep(1, chNum); weight[flagCH] <- 0.0
    # availCH <- which(weight == 1.0)
    knotNum <- 120
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
#-------- BP Cal
##BPphsCal <- function( SPEC, BP ){
#    BPphs <- complex( modulus=1, argument=-Arg(BP))
#    return( SPEC* BPphs )
#}
#-------- Function to Calibrate Delay and Phase
##DelayPhaseCal <- function( scanSpec, mjdSec, delayFit, ReFit, ImFit){
#    temp <- scanSpec
#    phase <- atan2(predict(ImFit, mjdSec)$y, predict(ReFit, mjdSec)$y)
#    delay <- predict(delayFit, mjdSec)$y
#    for(timeIndex in 1:ncol(scanSpec)){
#        scanSpec[,timeIndex] <- delayPhase_cal(temp[,timeIndex], delay[timeIndex], -phase[timeIndex])
#    }
#    return(scanSpec)
#}
#-------- Spectral Mitigation
SPmitigation <- function( spec, SPCH ){
    for(index in 1:length(SPCH)){
        spec[SPCH[index]] <- mean(spec[c((SPCH[index]-8):(SPCH[index]-1), (SPCH[index]+1):(SPCH[index]+8))])
    }
    return(spec)
}
#-------- Smoothed Bandpass Calibration
SBCspec <- function( spec, knotNum, weight, mitigCH){
    if( length(mitigCH) == 0){
        return(predict(smooth.spline( spec, w=weight, all.knots=F, nknots=knotNum), 1:length(spec))$y)
    } else {
        return(predict(smooth.spline( SPmitigation(spec, mitigCH), w=weight, all.knots=F, nknots=knotNum), 1:length(spec))$y)
    }
}
#-------- ScanTime
scanTime   <- function(MJD_df){ return((MJD_df[[1]] + MJD_df[[2]])/2 ) }
scanIntegT <- function(MJD_df){ return( MJD_df[[2]] - MJD_df[[2]] + 1) }
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

