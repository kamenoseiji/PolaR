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
smoothWidth <- 384; knotNum <- floor(chNum / smoothWidth)
#-------- Save into file
PDFfilename <- sprintf("%s.Zeeman.pdf", strsplit(args[1], "\\.")[[1]][1])
#-------- Plot Stokes I for HC3N
pdf(PDFfilename)
#plotFreq <- c(2.0, 2.5) # range in MHz
#lineFreq <- c(2.15, 2.4)
plotFreq <- c(1.5, 1.95) # range in MHz
lineFreq <- c(1.65, 1.87)
#plotFreq <- c(0.1, 0.5)
#lineFreq <- c(0.25, 0.35)
#plotFreq <- c(3.4, 3.7)
#lineFreq <- c(3.5, 3.62)
plotRange <- which.min(abs(freq - plotFreq[1])):which.min(abs(freq - plotFreq[2]))
lineRange <- which.min(abs(freq - lineFreq[1])):which.min(abs(freq - lineFreq[2]))
fitStokesI <- smooth.spline(freq, StokesI02, w=weight, all.knots=F, nknots=4*knotNum)
predStokesV <- predict(fitStokesI, (freq + 5e-6))$y - predict(fitStokesI, (freq - 5e-6))$y
predStokesI <- predict(fitStokesI, freq)$y
plotBunch <- 16
plot( bunch_vec(freq[plotRange],plotBunch)-0.5*plotBunch*chSep, bunch_vec(StokesI02[plotRange], plotBunch), type='s', xlab='Frequency [MHz]', ylab='Stokes I [K]', main='TMC-1 HC3N')
lines(freq[plotRange], predict(fitStokesI, freq[plotRange])$y, type='s', col='red')
#-------- Plot Stokes V for HC3N
plotBunch <- 32
err <- sd(StokesV02[SDrange]) / sqrt(plotBunch)
fit <- lm(formula=z~1+x+y, data=data.frame(x=predStokesV[lineRange], y=predStokesI[lineRange], z=StokesV02[lineRange]))
summary(fit)
plotX <- bunch_vec(freq[plotRange], plotBunch)
plotY <- bunch_vec(StokesV02[plotRange] - fit$coefficients[3]* predStokesI[plotRange] - fit$coefficients[1], plotBunch)
plot( plotX, plotY , pch=20, ylim=c(-1e-1, 1e-1), xlab='Frequency [MHz]', ylab='Stokes V [K]', main='TMC-1 HC3N')
arrows( plotX, plotY - err, plotX, plotY + err, angle=90, length=0)
lines(freq[plotRange]-0.5*chSep, predStokesV[plotRange]*fit$coefficients[2], type='s', col='red')
legend(min(freq[plotRange]), 0.1, legend=sprintf('Zeeman Shift = %5.1f ± %4.1f Hz', 10*fit[[1]][2], 10*summary(fit)[[4]][5]))

#-------- Plot Stokes I for CCS
#plotFreq <- c(1.8, 2.2) # range in MHz
#lineFreq <- c(1.92, 2.1)
plotFreq <- c(1.65, 1.95) # range in MHz
lineFreq <- c(1.74, 1.87)
plotRange <- which.min(abs(freq - plotFreq[1])):which.min(abs(freq - plotFreq[2]))
lineRange <- which.min(abs(freq - lineFreq[1])):which.min(abs(freq - lineFreq[2]))
fitStokesI <- smooth.spline(freq, StokesI13, w=weight, all.knots=F, nknots=4*knotNum)
predStokesV <- predict(fitStokesI, (freq + 5e-6))$y - predict(fitStokesI, (freq - 5e-6))$y
predStokesI <- predict(fitStokesI, freq)$y
plotBunch <- 16
plot( bunch_vec(freq[plotRange],plotBunch)-0.5*plotBunch*chSep, bunch_vec(StokesI13[plotRange], plotBunch), type='s', xlab='Frequency [MHz]', ylab='Stokes I [K]', main='TMC-1 CCS')
lines(freq[plotRange], predict(fitStokesI, freq[plotRange])$y, type='s', col='red')
#-------- Plot Stokes V for CCS
plotBunch <- 32
err <- sd(StokesV13[SDrange]) / sqrt(plotBunch)
fit <- lm( formula=z~1+x+y, data=data.frame(x=predStokesV[lineRange], y=predStokesI[lineRange], z=StokesV13[lineRange]))
summary(fit)
plotX <- bunch_vec(freq[plotRange], plotBunch)
plotY <- bunch_vec(StokesV13[plotRange] - fit$coefficients[3]* predStokesI[plotRange] - fit$coefficients[1], plotBunch)
plot( plotX, plotY , pch=20, ylim=c(-1e-1, 1e-1), xlab='Frequency [MHz]', ylab='Stokes V [K]', main='TMC-1 CCS')
arrows( plotX, plotY - err, plotX, plotY + err, angle=90, length=0)
lines(freq[plotRange]-0.5*chSep, predStokesV[plotRange]*fit$coefficients[2], type='s', col='red')
legend(min(freq[plotRange]), 0.1, legend=sprintf('Zeeman Shift = %5.1f ± %4.1f Hz', 10*fit[[1]][2], 10*summary(fit)[[4]][5]))
dev.off()

