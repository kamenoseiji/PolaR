# StokeSpec
# usage: Rscript StokesSpec [Scan.Rdata file name] [SPEC.Rdata file name] [WG.Rdata file name] [BP file name]
#
library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readPolariS.R", ssl.verifypeer = FALSE)))
setwd('.')
#-------- Load Spec and Scan data
args <- commandArgs(trailingOnly = T)
setwd('.')
#setwd('/Volumes/SSD/PolariS/20150317/')
#args <- c('2015076052841.Scan.Rdata', '2015076052912.SPEC.Rdata', '2015076035301.WG.Rdata', '2015076035301.BP.Rdata')
#args <- c('2015076043056.Scan.Rdata', '2015076043132.SPEC.Rdata', '2015076035301.WG.Rdata', '2015076035301.BP.Rdata')
load(args[1])	 #Load Stokes file
# freq_track <- 45379033000 # [Hz]
# freq_center <- c(45490316000, 45379033000, 45490316000, 45379033000)
# freq_IF <- 5208	# [MHz]
# freq_1stLO <- freq_track - freq_IF*1e6	#[Hz]
# freq_BBC <- c(5317.0, 5206.0, 5317.0, 5206.0)
# freq_rest   <- freq_center
# velSrc      <- 5900 # m/s
chNum <- length(StokesI02)
SDrange <- 8193:16384
chRange <- floor(0.05*chNum):floor(0.95*chNum)
freq <- (0:(chNum-1))/chNum* 4.0	# MHz
# veloc0 <- (2.0 - freq)*1.0e6 /  freq_center[1] * 299792458 + 
chSep <- 4.0 / chNum
weight <- rep(0, chNum); weight[chRange] = 1.0
smoothWidth <- 256; knotNum <- floor(chNum / smoothWidth)
#-------- Save into file
PDFfilename <- sprintf("%s.Zeeman.pdf", strsplit(args[1], "\\.")[[1]][1])
#-------- Plot Stokes I for HC3N
pdf(PDFfilename)
plotFreq <- c(1.4, 2.0) # range in MHz
lineFreq <- c(1.6, 1.75)
baseFreq <- c(1.4, 1.5, 1.85, 2.0)
plotRange <- which.min(abs(freq - plotFreq[1])):which.min(abs(freq - plotFreq[2]))
lineRange <- which.min(abs(freq - lineFreq[1])):which.min(abs(freq - lineFreq[2]))
baseRange <- c(which.min(abs(freq - baseFreq[1])):which.min(abs(freq - baseFreq[2])), which.min(abs(freq - baseFreq[3])):which.min(abs(freq - baseFreq[4])))
StokesI02 <- StokesI02 - mean(StokesI02[baseRange])
StokesI13 <- StokesI13 - mean(StokesI13[baseRange])
fitStokesI <- smooth.spline(freq, StokesI02, w=weight, all.knots=F, nknots=4*knotNum)
fitStokesQ <- smooth.spline(freq, StokesQ02, w=weight, all.knots=F, nknots=4*knotNum)
fitStokesU <- smooth.spline(freq, StokesU02, w=weight, all.knots=F, nknots=4*knotNum)
predStokesV <- predict(fitStokesI, (freq + 5e-6))$y - predict(fitStokesI, (freq - 5e-6))$y
predStokesI <- predict(fitStokesI, freq)$y
predStokesQ <- predict(fitStokesQ, freq)$y
predStokesU <- predict(fitStokesU, freq)$y
plotBunch <- 32
fitBunch  <- 8
cols   <- c('black', 'blue', 'green')
labels <- c('I', 'Q', 'U')
plot( bunch_vec(freq[plotRange],plotBunch)-0.5*plotBunch*chSep, bunch_vec(StokesI02[plotRange], plotBunch), type='s', xlab='Frequency [MHz]', ylab='Stokes I [K]', main='OMC-2 CH3OH')
lines( bunch_vec(freq[plotRange],plotBunch)-0.5*plotBunch*chSep, bunch_vec(StokesQ02[plotRange], plotBunch), type='s', col=cols[2])
lines( bunch_vec(freq[plotRange],plotBunch)-0.5*plotBunch*chSep, bunch_vec(StokesU02[plotRange], plotBunch), type='s', col=cols[3])
lines( bunch_vec(freq[plotRange], fitBunch), bunch_vec(predict(fitStokesI, freq[plotRange])$y, fitBunch), col='orange')
legend("topleft", legend=labels, col=cols, lty=rep(1,3))
abline(h=0)
#-------- Plot Stokes V for HC3N
plotBunch <- 64
err <- sd(StokesV02[SDrange]) / sqrt(plotBunch)
fit <- lm(formula=z~1+x+y, data=data.frame(x=predStokesV[lineRange], y=predStokesI[lineRange], z=StokesV02[lineRange]))
summary(fit)
plotX <- bunch_vec(freq[plotRange], plotBunch)
plotY <- bunch_vec(StokesV02[plotRange] - fit$coefficients[1] - fit$coefficients[3]* predStokesI[plotRange], plotBunch)
plot( plotX, plotY , pch=20, ylim=c(-2.5e-1, 2.5e-1), xlab='Frequency [MHz]', ylab='Stokes V [K]', main='OMC-2 CH3OH')
arrows( plotX, plotY - err, plotX, plotY + err, angle=90, length=0)
lines( bunch_vec(freq[plotRange], fitBunch), bunch_vec(predStokesV[plotRange]*fit$coefficients[2], fitBunch), col='red')
legend(min(freq[plotRange]), 0.2, legend=sprintf('Zeeman Shift = %5.1f ± %4.1f Hz', 10*fit[[1]][2], 10*summary(fit)[[4]][5]))

#-------- Plot Stokes I for CH3OH
plotFreq <- c(1.4, 2.0) # range in MHz
lineFreq <- c(1.6, 1.75)
plotRange <- which.min(abs(freq - plotFreq[1])):which.min(abs(freq - plotFreq[2]))
lineRange <- which.min(abs(freq - lineFreq[1])):which.min(abs(freq - lineFreq[2]))
fitStokesI <- smooth.spline(freq, StokesI13, w=weight, all.knots=F, nknots=4*knotNum)
fitStokesQ <- smooth.spline(freq, StokesQ13, w=weight, all.knots=F, nknots=4*knotNum)
fitStokesU <- smooth.spline(freq, StokesU13, w=weight, all.knots=F, nknots=4*knotNum)
predStokesV <- predict(fitStokesI, (freq + 5e-6))$y - predict(fitStokesI, (freq - 5e-6))$y
predStokesI <- predict(fitStokesI, freq)$y
predStokesQ <- predict(fitStokesQ, freq)$y
predStokesU <- predict(fitStokesU, freq)$y
plotBunch <- 32
plot( bunch_vec(freq[plotRange],plotBunch)-0.5*plotBunch*chSep, bunch_vec(StokesI13[plotRange], plotBunch), type='s', xlab='Frequency [MHz]', ylab='Stokes I [K]', main='OMC-2 CH3OH')
lines( bunch_vec(freq[plotRange],plotBunch)-0.5*plotBunch*chSep, bunch_vec(StokesQ13[plotRange], plotBunch), type='s', col=cols[2])
lines( bunch_vec(freq[plotRange],plotBunch)-0.5*plotBunch*chSep, bunch_vec(StokesU13[plotRange], plotBunch), type='s', col=cols[3])
lines( bunch_vec(freq[plotRange], fitBunch), bunch_vec(predict(fitStokesI, freq[plotRange])$y, fitBunch), col='orange')
#-------- Plot Stokes V for CH3OH
plotBunch <- 64
err <- sd(StokesV13[SDrange]) / sqrt(plotBunch)
fit <- lm( formula=z~1+x+y, data=data.frame(x=predStokesV[lineRange], y=predStokesI[lineRange], z=StokesV13[lineRange]))
summary(fit)
plotX <- bunch_vec(freq[plotRange], plotBunch)
plotY <- bunch_vec(StokesV13[plotRange] - fit$coefficients[1] - fit$coefficients[3]* predStokesI[plotRange], plotBunch)
plot( plotX, plotY , pch=20, ylim=c(-2.5e-1, 2.5e-1), xlab='Frequency [MHz]', ylab='Stokes V [K]', main='OMC-2 CH3OH')
arrows( plotX, plotY - err, plotX, plotY + err, angle=90, length=0)
lines( bunch_vec(freq[plotRange], fitBunch), bunch_vec(predStokesV[plotRange]*fit$coefficients[2], fitBunch), col='red')
legend(min(freq[plotRange]), 0.2, legend=sprintf('Zeeman Shift = %5.1f ± %4.1f Hz', 10*fit[[1]][2], 10*summary(fit)[[4]][5]))
dev.off()

