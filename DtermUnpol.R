# ContStokes.R
# usage: Rscript ContStokes.R [Scan.Rdata file name] [SPEC.Rdata file name] [WG.Rdata file name] [BP file name]
#
RPATH <- '~/Programs/PolaR'
FuncList <- c('date', 'PolariCalib')
source(sprintf('%s/loadModule.R', RPATH))

Err <- try(loadGitHub( FuncList ), silent=FALSE)
if(class(Err) == "try-error"){ loadLocal( RPATH, FuncList ) }

setwd('.')
options(digits = 4)
#-------- ScanTime
scanTime   <- function(MJD_df){ return((MJD_df[[1]] + MJD_df[[2]])/2 ) }
scanIntegT <- function(MJD_df){ return( MJD_df[[2]] - MJD_df[[2]] + 1) }


#-------- Load Spec and Scan data
args <- commandArgs(trailingOnly = T)
setwd('.')
load(args[1])	 #Load Scan file
load(args[2])	 #Load SPEC file
load(args[3])	 #Load delay file
load(args[4])	 #Load BP file

#
#-------- Smoothed Delay and Phase
if( length(WG$mjdSec) > 3 ){
    delay00Fit <- smooth.spline(WG$mjdSec, WG$delay00, spar=0.25); delay01Fit <- smooth.spline(WG$mjdSec, WG$delay01, spar=0.25)
    Re00Fit <- smooth.spline(WG$mjdSec, Re(WG$Vis00), spar=0.25); Im00Fit <- smooth.spline(WG$mjdSec, Im(WG$Vis00), spar=0.25)
    Re01Fit <- smooth.spline(WG$mjdSec, Re(WG$Vis01), spar=0.25); Im01Fit <- smooth.spline(WG$mjdSec, Im(WG$Vis01), spar=0.25)
} else {
    delay00Fit <- lm(formula = y ~ x, data.frame(x=WG$mjdSec, y=WG$delay00));
    delay01Fit <- lm(formula = y ~ x, data.frame(x=WG$mjdSec, y=WG$delay01));
    Re00Fit <- lm(formula = y ~ x, data.frame(x=WG$mjdSec, y=Re(WG$Vis00)));
    Im00Fit <- lm(formula = y ~ x, data.frame(x=WG$mjdSec, y=Im(WG$Vis00)));
    Re01Fit <- lm(formula = y ~ x, data.frame(x=WG$mjdSec, y=Re(WG$Vis01)));
    Im01Fit <- lm(formula = y ~ x, data.frame(x=WG$mjdSec, y=Im(WG$Vis01)));
}
#
#-------- Initial Parameters
chNum <- dim(on_A00)[1]
chRange <- floor(0.05*chNum):floor(0.95*chNum)

OnIndex <- which(Scan$scanType == 'ON')
OfIndex <- which(Scan$scanType == 'OFF')

#-------- Amplitude calibration of Autocorr
fileName <- sprintf("%s.ContPower.pdf", strsplit(args[2], "\\.")[[1]][1])
pdf(fileName)
Ta00 <- TaCal( colSums(on_A00[chRange,]),  colSums(off_A00[chRange,]), Scan$Tsys00[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]) #; cat(Scan$Tsys00[OfIndex]); cat('\n') # ; cat(Ta00); cat('\n')
Ta01 <- TaCal( colSums(on_A01[chRange,]),  colSums(off_A01[chRange,]), Scan$Tsys01[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]) #; cat(Scan$Tsys01[OfIndex]); cat('\n') # ; cat(Ta01); cat('\n')
Ta02 <- TaCal( colSums(on_A02[chRange,]),  colSums(off_A02[chRange,]), Scan$Tsys02[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]) #; cat(Scan$Tsys02[OfIndex]); cat('\n') # ; cat(Ta02); cat('\n')
Ta03 <- TaCal( colSums(on_A03[chRange,]),  colSums(off_A03[chRange,]), Scan$Tsys03[OfIndex], scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex]) #; cat(Scan$Tsys03[OfIndex]); cat('\n') # ; cat(Ta03); cat('\n')

#-------- Delay, phase, bandpass, and amplitude calibration for CrossCorr
Tx02 <- TxCal( colSums( DelayPhaseCal( BPphsCal(on_C00, BP$BP00), scanTime(onMJD), delay00Fit, Re00Fit, Im00Fit)[chRange,]),
               colSums( DelayPhaseCal( BPphsCal(off_C00, BP$BP00), scanTime(offMJD), delay00Fit, Re00Fit, Im00Fit)[chRange,]),
			   sqrt( colSums(off_A00[chRange,])* colSums(off_A02[chRange,]) ),
			   sqrt(Scan$Tsys00[OfIndex]*  Scan$Tsys02[OfIndex]),
			   scanTime(onMJD), scanTime(offMJD), Scan$mjdSec[OfIndex])
Tx13 <- TxCal( colSums( DelayPhaseCal( BPphsCal(on_C01, BP$BP01), scanTime(onMJD), delay01Fit, Re01Fit, Im01Fit)[chRange,]),
               colSums( DelayPhaseCal( BPphsCal(off_C01, BP$BP01), scanTime(offMJD), delay01Fit, Re01Fit, Im01Fit)[chRange,]),
			   sqrt( colSums(off_A01[chRange,])* colSums(off_A03[chRange,]) ),
			   sqrt(Scan$Tsys01[OfIndex]*  Scan$Tsys03[OfIndex]), scanTime(onMJD),
			   scanTime(offMJD), Scan$mjdSec[OfIndex])
#-------- Parallactic angle
AZ <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$AZ[OnIndex], spar=0.5), scanTime(onMJD))$y
EL <- predict(smooth.spline(Scan$mjdSec[OnIndex], Scan$EL[OnIndex], spar=0.5), scanTime(onMJD))$y
Pang <- -azel2pa(AZ, EL) + EL*pi/180 - pi/2 
#-------- Dx + Dy*, assuming Q = U = V = 0
D <- data.frame( XY02 = Tx02/sqrt(Ta00*Ta02), XY13 = Tx13/sqrt(Ta01*Ta03), Rxy02 = Ta00/Ta02, Rxy13 = Ta01/Ta03, Gxy02 = sqrt(Ta00*Ta02), Gxy13 = sqrt(Ta01*Ta03), AZ = AZ, EL = EL, mjdSec = scanTime(onMJD) )
D
fileName <- sprintf("%s.Dterm.Rdata", strsplit(args[2], "\\.")[[1]][1])
save(D, file=fileName)
cat(sprintf('D-term is saved into %s\n', fileName))


