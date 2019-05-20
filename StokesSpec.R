# StokeSpec
# usage: Rscript StokesSpec [Scan.Rdata file name] [SPEC.Rdata file name] [WG.Rdata file name] [BP file name]
# Options:
#  -b and -B : Beam squint in AZ and EL (at feed horn), i.e. dAZ1 and dEL1 
#  -o and -O : Beam offset in AZ and EL (at dish), i.e. dAZ0 and dEL0
#              Thus, the beam separation between RHCP and LHCP on the sky will be given as 
#                 dAZ = dAZ0 + dAZ1* cos EL + dEL1* sin(EL)
#                 dEL = dEL0 - dAZ1* sin EL + dEL1* cos(EL)
#
RPATH <- '~/Programs/PolaR'
FuncList <- c('readPolariS', 'readSAM45', 'date', 'Qeff', 'PolariCalib')
source(sprintf('%s/loadModule.R', RPATH))
library(RCurl)

funcNum <- length(FuncList)
for( index in 1:funcNum){
    URL <- sprintf("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/%s.R", FuncList[index])
    Err <- try( eval(parse(text = getURL(URL, ssl.verifypeer = FALSE))), silent=FALSE)
}   
if(class(Err) == "try-error"){ loadLocal( RPATH, FuncList ) }
setwd('.')

#-------- Parse command-line arguments
parseArg <- function( args ){
    argNum <- length(args)
    fileNum <- argNum
    smoothWidth <- 128      # Set Default Value
    DtermFile   <- FALSE    # Set Default Value
    for( index in 1:argNum ){
        if(substr(args[index], 1,2) == "-S"){ smoothWidth <- as.numeric(substring(args[index], 3));  fileNum <- fileNum - 1}
        if(substr(args[index], 1,2) == "-D"){ DtermFile <- substring(args[index], 3);  fileNum <- fileNum - 1}
        if(substr(args[index], 1,2) == "-b"){ SqX <- as.numeric(substring(args[index], 3));  fileNum <- fileNum - 1}        # Beam Squint (RHCP-LHCP) in AZ [arcsec]
        if(substr(args[index], 1,2) == "-B"){ SqY <- as.numeric(substring(args[index], 3));  fileNum <- fileNum - 1}        # Beam Squint (RHCP-LHCP) in EL [arcsec]
        if(substr(args[index], 1,2) == "-o"){ OfX <- as.numeric(substring(args[index], 3));  fileNum <- fileNum - 1}        # Beam offset (RHCP-LHCP) in AZ [arcsec]
        if(substr(args[index], 1,2) == "-O"){ OfY <- as.numeric(substring(args[index], 3));  fileNum <- fileNum - 1}        # Beam offset (RHCP-LHCP) in EL [arcsec]
        if(substr(args[index], 1,2) == "-v"){ VgRA  <- as.numeric(substring(args[index], 3));  fileNum <- fileNum - 1}      # Velocity Gradient in RA [Hz/arcsec]
        if(substr(args[index], 1,2) == "-V"){ VgDEC <- as.numeric(substring(args[index], 3));  fileNum <- fileNum - 1}      # Velocity Gradient in DEC [Hz/arcsec]
        if(substr(args[index], 1,2) == "-p"){ minPA <- as.numeric(substring(args[index], 3));  fileNum <- fileNum - 1}
        if(substr(args[index], 1,2) == "-P"){ maxPA <- as.numeric(substring(args[index], 3));  fileNum <- fileNum - 1}
    }
    fileName <- args[(argNum - fileNum + 1):argNum]
    return( list(smoothWidth = smoothWidth, DtermFile = DtermFile, SqX = SqX, SqY = SqY, OfX = OfX, OfY = OfY, VgRA = VgRA, VgDEC = VgDEC, minPA = pi*minPA/180, maxPA = pi*maxPA/180, fileName = fileName) )
}
#-------- Load Spec and Scan data
args <- parseArg(commandArgs(trailingOnly = T))
#args <- parseArg(c('-S64', '-D2016023002950.Dcomb.Rdata', '-o2.04', '-O0.39', '-b0.41', '-B1.25', '-v0', '-V0', '-p-180', '-P180', '2016022232525.Scan.Rdata','2016022232602.SPEC.Rdata', '2016022225647.WG.Rdata', '2016022225647.BP.Rdata'))
setwd('.')
cat('Loading spectral data...\n')
load(args$fileName[1])	 #Load Scan file
load(args$fileName[2])	 #Load SPEC file
load(args$fileName[3])	 #Load delay file
load(args$fileName[4])	 #Load BP file
if(args$DtermFile != FALSE){ load(args$DtermFile)} else {D <- data.frame(XY02 = c(0+0i), XY13 = c(0+0i), Rxy02 = c(1.0), Rxy13 = c(1.0), Gxy02 = c(1.0), Gxy13 = c(1.0))}
#
#-------- Smoothed Delay and Phase
cat('--- Delay and Phase Cal\n')
delayNum <- length(WG$mjdSec)
if(delayNum < 5){
    WG <- rbind( WG[1,], WG[1,], WG, WG[delayNum,], WG[delayNum,])
    WG$mjdSec[1] <- WG$mjdSec[1] - 3600
    WG$mjdSec[2] <- WG$mjdSec[2] - 1800
    WG$mjdSec[delayNum + 1] <- WG$mjdSec[delayNum + 1] + 1800
    WG$mjdSec[delayNum + 2] <- WG$mjdSec[delayNum + 2] + 3600
}
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
#mitigCH <- c(16385)
#mitigCH <- c(23359, 32769)
mitigCH <- c(60424)
flagCH <- unique(c(1, 2, 4, mitigCH))
weight <- rep(1, chNum); weight[flagCH] <- 0.0
# smoothWidth <- 512; knotNum <- floor(chNum / smoothWidth)
smoothWidth <- args$smoothWidth; knotNum <- floor(chNum / smoothWidth)
#-------- Scan Pattern
OnIndex <- which(Scan$scanType == 'ON')
OfIndex <- which(Scan$scanType == 'OFF')
D_index <- which( D$Gxy02 + D$Gxy13 > 0.9* median( D$Gxy02 + D$Gxy13 ))	# Flag pointing-error out
Rxy02 <- mean(D$Rxy02[D_index])
Rxy13 <- mean(D$Rxy13[D_index])
Dxy02 <- mean(D$XY02[D_index])
Dxy13 <- mean(D$XY13[D_index])
#-------- Tsys at each scan
cat('--- Tsys smoothing\n')
Tsys00 <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$Tsys00[OnIndex], spar=0.5), scanTime(onMJD))$y
Tsys01 <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$Tsys01[OnIndex], spar=0.5), scanTime(onMJD))$y
Tsys02 <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$Tsys02[OnIndex], spar=0.5), scanTime(onMJD))$y
Tsys03 <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$Tsys03[OnIndex], spar=0.5), scanTime(onMJD))$y

#-------- Radial Velocity
#DeltaFreq <- 317.00 - 206.00	# [MHz], LO(CH1) - LO(CH2)	// 2014.4.17
#DeltaFreq <- 317.50 - 206.20	# [MHz], LO(CH1) - LO(CH2)  // 2015.3.15 - 17
DeltaFreq <- 0.0
Vrad <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$Vrad[OnIndex], spar=0.5), scanTime(onMJD))$y
chShift <- round(DeltaFreq* Vrad / 299792458 / chSep)		# Number of channels to shift

#-------- Parallactic angle
AZ <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$AZ[OnIndex], spar=0.5), scanTime(onMJD))$y
EL <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$EL[OnIndex], spar=0.5), scanTime(onMJD))$y
Pang <- -azel2pa(AZ, EL) + EL*pi/180 - pi/2 
cs <- cos(Pang)
sn <- sin(Pang)
PAindex <- which( (azel2pa(AZ, EL) > args$minPA) & (azel2pa(AZ, EL) < args$maxPA) )
#-------- Beam Squint and Velocity Gradient
#BeamSquintAzEl <- c(args$SqX,  args$SqY)    # arcsec
VelocGradRADEC <- c(args$VgRA, args$VgDEC)  # Hz / arcsec
#BeamSquintAzEl <- c(-1.67, 1.18)    # arcsec
#VelocGradRADEC <- c(152, 115)       # Hz / arcsec
#-------- Amplitude calibration of Autocorr
cat('--- Amplitude and Phase Calibration\n')
Ta00 <- TaCalSpec(on_A00, off_A00, scanTime(onMJD), scanTime(offMJD), Tsys00, weight, mitigCH) * sqrt(Rxy02)
Ta01 <- TaCalSpec(on_A01, off_A01, scanTime(onMJD), scanTime(offMJD), Tsys01, weight, mitigCH) * sqrt(Rxy13)
Ta02 <- TaCalSpec(on_A02, off_A02, scanTime(onMJD), scanTime(offMJD), Tsys02, weight, mitigCH) / sqrt(Rxy02)
Ta03 <- TaCalSpec(on_A03, off_A03, scanTime(onMJD), scanTime(offMJD), Tsys03, weight, mitigCH) / sqrt(Rxy13)
#
#-------- Delay, phase, bandpass, and amplitude calibration for CrossCorr
Tx02 <- TxCalSpec( 
	DelayPhaseCal( BPphsCal(on_C00, BP$BP00), scanTime(onMJD), delay00Fit, Re00Fit, Im00Fit),
	DelayPhaseCal( BPphsCal(off_C00, BP$BP00), scanTime(offMJD), delay00Fit, Re00Fit, Im00Fit),
	off_A00, off_A02, scanTime(onMJD), scanTime(offMJD), sqrt(Tsys00 * Tsys02), weight, mitigCH)

Tx13 <- TxCalSpec( 
	DelayPhaseCal( BPphsCal(on_C01, BP$BP01), scanTime(onMJD), delay01Fit, Re01Fit, Im01Fit),
	DelayPhaseCal( BPphsCal(off_C01, BP$BP01), scanTime(offMJD), delay01Fit, Re01Fit, Im01Fit),
	off_A01, off_A03, scanTime(onMJD), scanTime(offMJD), sqrt(Tsys01 * Tsys03), weight, mitigCH)
#-------- Doppler Tracking
#Ta00 <- deDopp(Ta00, chShift)
#Ta02 <- deDopp(Ta02, chShift)
#Tx02 <- deDopp(Tx02, chShift)
#-------- Time Integration to determine Stokes spectrum
cat('--- Full Stokes with Beam Squint Correction\n')
weight0 <- rep(0.0, length(Tsys00)); weight0[PAindex] <- 1.0 / Tsys00[PAindex]^2
weight1 <- rep(0.0, length(Tsys01)); weight1[PAindex] <- 1.0 / Tsys01[PAindex]^2
weight2 <- rep(0.0, length(Tsys02)); weight2[PAindex] <- 1.0 / Tsys02[PAindex]^2
weight3 <- rep(0.0, length(Tsys03)); weight3[PAindex] <- 1.0 / Tsys03[PAindex]^2
weight02 <- rep(0.0, length(Tsys00)); weight02[PAindex] <- 1.0 / (Tsys00[PAindex] * Tsys02[PAindex])
weight13 <- rep(0.0, length(Tsys01)); weight13[PAindex] <- 1.0 / (Tsys01[PAindex] * Tsys03[PAindex])
StokesI02 <- 0.5*(rowSums(Ta00* weight0)/sum(weight0) + rowSums(Ta02* weight2)/sum(weight2))
StokesI13 <- 0.5*(rowSums(Ta01* weight1)/sum(weight1) + rowSums(Ta03* weight3)/sum(weight3))

StokesQ02 <- 0.5* (rowSums(Ta00* cs* weight0)/sum(weight0) - rowSums(Ta02* cs* weight2)/sum(weight2)) - rowSums((Re(Tx02) - Re(Dxy02)* StokesI02)* sn* weight02) / sum(weight02)
StokesQ13 <- 0.5* (rowSums(Ta01* cs* weight1)/sum(weight1) - rowSums(Ta03* cs* weight3)/sum(weight3)) - rowSums((Re(Tx13) - Re(Dxy13)* StokesI13)* sn* weight13) / sum(weight13)

StokesU02 <- 0.5* (rowSums(Ta00* sn* weight0)/sum(weight0) - rowSums(Ta02* sn* weight2)/sum(weight2)) + rowSums((Re(Tx02) - Re(Dxy02)* StokesI02)* cs* weight02) / sum(weight02)
StokesU13 <- 0.5* (rowSums(Ta01* sn* weight1)/sum(weight1) - rowSums(Ta03* sn* weight3)/sum(weight3)) + rowSums((Re(Tx13) - Re(Dxy13)* StokesI13)* cs* weight13) / sum(weight13)

#lineRange <- c(0.1, 2.6, 2.8, 3.9)
#chRange <- c( which(freq > lineRange[1] & freq < lineRange[2]), which(freq > lineRange[3] & freq < lineRange[4]))

fshift <- function(AZEL, sqX, sqY, ofX, ofY, vg){
    PA <- azel2pa(AZEL$AZ, AZEL$EL)
    PAel <- PA - pi* AZEL$EL / 180.0
    csPA <- cos(PA)
    snPA <- sin(PA)
    csPAel <- cos(PAel)
    snPAel <- sin(PAel)
    dRA <- -ofX*csPA + ofY* snPA - sqX*csPAel + sqY* snPAel
    dEL <-  ofX*snPA + ofY* csPA + sqX*snPAel + sqY* csPAel
    Fshift <- vg[1]* dRA + vg[2]* dEL
    return(Fshift)
}

weight <- rep(0.0, length(freq)); weight[chRange] <- 1.0
fitStokesI02 <- smooth.spline(freq, StokesI02, w=weight, all.knots=F, nknots=knotNum)
fitStokesI13 <- smooth.spline(freq, StokesI13, w=weight, all.knots=F, nknots=knotNum)

AZEL <- data.frame(AZ = AZ, EL = EL)
FShift <- fshift(AZEL, args$OfX, args$OfY, args$SqX, args$SqY, VelocGradRADEC)

for(scan_index in 1:length(Tsys00)){
    #PA <- azel2pa(AZ[scan_index], EL[scan_index]); csPA <- cos(PA); snPA <- sin(PA)
    #PAel <- PA + pi* EL[scan_index] / 180.0; csPAel <- cos(PAel); snPAel <- sin(PAel)
    #squintRADEC <- c( -args$OfX*csPA + args$OfY* snPA - args$SqX*csPAel + args$SqY* snPAel, args$OfX*snPA + args$OfY* csPA + args$SqX*snPAel + args$SqY* csPAel )
    #Fshift <- VelocGradRADEC %*% squintRADEC
    #cat(sprintf("Scan%d AZ=%5.1f EL=%5.1f PA=%5.1f Fshift=%5.1f Hz\n", scan_index, AZ[scan_index], EL[scan_index], PA, Fshift))
    cat(sprintf("Scan%d AZ=%5.1f EL=%5.1f PA=%5.1f Fshift=%5.1f Hz\n", scan_index, AZ[scan_index], EL[scan_index], azel2pa(AZ[scan_index], EL[scan_index]), FShift[scan_index]))
    Fshift <- FShift[scan_index]* 0.5e-6                                                                                    # Frequency Shift in MHz
    fakeStokesV02 <- predict(fitStokesI02, (freq + Fshift))$y - predict(fitStokesI02, (freq - Fshift))$y            # Velocity gradient correction
    fakeStokesV13 <- predict(fitStokesI13, (freq + Fshift))$y - predict(fitStokesI13, (freq - Fshift))$y            # Velocity gradient correction
    Tx02[chRange,scan_index] = Tx02[chRange,scan_index] - 0.5i * fakeStokesV02[chRange]
    Tx13[chRange,scan_index] = Tx13[chRange,scan_index] - 0.5i * fakeStokesV13[chRange]
}
StokesV02 <- rowSums(Im(Tx02)* weight02)/sum(weight02) - Im(Dxy02)* StokesI02
StokesV13 <- rowSums(Im(Tx13)* weight13)/sum(weight13) - Im(Dxy13)* StokesI13
#-------- Save into file
fileName <- sprintf("%s.STOKES.Rdata", strsplit(args$fileName[1], "\\.")[[1]][1])
save(StokesI02, StokesI13, StokesQ02, StokesQ13, StokesU02, StokesU13, StokesV02, StokesV13, file=fileName)
cat('Stokes spectra are saved into '); cat(fileName); cat('\n')
