library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readPolariS.R", ssl.verifypeer = FALSE)))
#
SDrange <- 8193:16384
shift02 <- 8192		# Shift 8192 ch
shift13 <- 3277		# Shift 3276.8 ch
load('/Volumes/SSD/PolariS/20140417/2014107013932.STOKES.Rdata')
Stokes_2014107013932 <- data.frame(
	I02 = c(StokesI02[(shift02+1):chNum], StokesI02[1:shift02]),
	I13 = c(StokesI13[(shift13+1):chNum], StokesI13[1:shift13]),
	Q02 = c(StokesQ02[(shift02+1):chNum], StokesQ02[1:shift02]),
	Q13 = c(StokesQ13[(shift13+1):chNum], StokesQ13[1:shift13]),
	U02 = c(StokesU02[(shift02+1):chNum], StokesU02[1:shift02]),
	U13 = c(StokesU13[(shift13+1):chNum], StokesU13[1:shift13]),
	V02 = c(StokesV02[(shift02+1):chNum], StokesV02[1:shift02]),
	V13 = c(StokesV13[(shift13+1):chNum], StokesV13[1:shift13]))
SD_2014107013932 <- sd(StokesV13[SDrange])
load('/Volumes/SSD/PolariS/20150315/2015074040721.STOKES.Rdata')
Stokes_2015074040721 <- data.frame(I02 = StokesI02, I13=StokesI13, Q02=StokesQ02, Q13=StokesQ13, U02=StokesU02, U13=StokesU13, V02=StokesV02, V13=StokesV13)
SD_2015074040721 <- sd(StokesV13[SDrange])
load('/Volumes/SSD/PolariS/20150316/2015075034519.STOKES.Rdata')
Stokes_2015075034519 <- data.frame(I02 = StokesI02, I13=StokesI13, Q02=StokesQ02, Q13=StokesQ13, U02=StokesU02, U13=StokesU13, V02=StokesV02, V13=StokesV13)
SD_2015075034519 <- sd(StokesV13[SDrange])
load('/Volumes/SSD/PolariS/20150317/2015076043132.STOKES.Rdata')
Stokes_2015076043132 <- data.frame(I02 = StokesI02, I13=StokesI13, Q02=StokesQ02, Q13=StokesQ13, U02=StokesU02, U13=StokesU13, V02=StokesV02, V13=StokesV13)
SD_2015076043132 <- sd(StokesV13[SDrange])
SD <- c(SD_2014107013932, SD_2015074040721, SD_2015075034519, SD_2015076043132)
WT <- 1/SD^2
#-------- Average
StokesI02 <- (WT[1]* Stokes_2014107013932$I02 + WT[2]* Stokes_2015074040721$I02 + WT[3]* Stokes_2015075034519$I02 + WT[4]* Stokes_2015076043132$I02) / sum(WT)
StokesI13 <- (WT[1]* Stokes_2014107013932$I13 + WT[2]* Stokes_2015074040721$I13 + WT[3]* Stokes_2015075034519$I13 + WT[4]* Stokes_2015076043132$I13) / sum(WT)
StokesV02 <- (WT[1]* Stokes_2014107013932$V02 + WT[2]* Stokes_2015074040721$V02 + WT[3]* Stokes_2015075034519$V02 + WT[4]* Stokes_2015076043132$V02) / sum(WT)
StokesV13 <- (WT[1]* Stokes_2014107013932$V13 + WT[2]* Stokes_2015074040721$V13 + WT[3]* Stokes_2015075034519$V13 + WT[4]* Stokes_2015076043132$V13) / sum(WT)

chNum <- length(StokesI02)
chRange <- floor(0.05*chNum):floor(0.95*chNum)
freq <- (0:(chNum-1))/chNum* 4.0	# MHz
chSep <- 4.0 / chNum
weight <- rep(1, chNum)
smoothWidth <- 384; knotNum <- floor(chNum / smoothWidth)
#-------- Save into file
PDFfilename <- "2015.Zeeman.pdf"
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
predStokesI <- predict(fitStokesI, freq)$y
predStokesV <- predict(fitStokesI, (freq + 5e-6))$y - predict(fitStokesI, (freq - 5e-6))$y
plotBunch <- 16
plot( bunch_vec(freq[plotRange],plotBunch)-0.5*plotBunch*chSep, bunch_vec(StokesI02[plotRange], plotBunch), type='s', xlab='Frequency [MHz]', ylab='Stokes I [K]', main='TMC-1 HC3N')
lines(freq[plotRange], predict(fitStokesI, freq[plotRange])$y, type='s', col='red')
#-------- Plot Stokes V for HC3N
plotBunch <- 32
fit <- lm(formula=z~1+x+y, data=data.frame(x=predStokesV[lineRange], y=predStokesI[lineRange], z=StokesV02[lineRange]))
plot( bunch_vec(freq[plotRange], plotBunch) - 0.5*plotBunch*chSep, bunch_vec(StokesV02[plotRange] - fit$coefficients[3]* predStokesI[plotRange] - fit$coefficients[1], plotBunch), type='s', ylim=c(-1e-1, 1e-1), xlab='Frequency [MHz]', ylab='Stokes V [K]', main='TMC-1 HC3N')
points( bunch_vec(freq[plotRange], plotBunch), bunch_vec(StokesV02[plotRange] - fit$coefficients[3]* predStokesI[plotRange] - fit$coefficients[1], plotBunch), pch=20, cex=0.5)
lines(freq[plotRange]-0.5*chSep, predStokesV[plotRange]*fit$coefficients[2], type='s', col='red')
legend(min(freq[plotRange]), 0.1, legend=sprintf('Zeeman Shift = %5.1f ± %4.1f Hz', 10*fit[[1]][2], 10*summary(fit)[[4]][5]))

#-------- Plot Stokes I for CCS
#plotFreq <- c(1.8, 2.2) # range in MHz
#lineFreq <- c(1.92, 2.1)
plotFreq <- c(1.65, 1.96) # range in MHz
lineFreq <- c(1.74, 1.87)
plotRange <- which.min(abs(freq - plotFreq[1])):which.min(abs(freq - plotFreq[2]))
lineRange <- which.min(abs(freq - lineFreq[1])):which.min(abs(freq - lineFreq[2]))
fitStokesI <- smooth.spline(freq, StokesI13, w=weight, all.knots=F, nknots=4*knotNum)
predStokesI <- predict(fitStokesI, freq)$y
predStokesV <- predict(fitStokesI, (freq + 5e-6))$y - predict(fitStokesI, (freq - 5e-6))$y
plotBunch <- 16
plot( bunch_vec(freq[plotRange],plotBunch)-0.5*plotBunch*chSep, bunch_vec(StokesI13[plotRange], plotBunch), type='s', xlab='Frequency [MHz]', ylab='Stokes I [K]', main='TMC-1 CCS')
lines(freq[plotRange], predict(fitStokesI, freq[plotRange])$y, type='s', col='red')
#-------- Plot Stokes V for CCS
plotBunch <- 220
err <- sd(StokesV13[SDrange]) / sqrt(plotBunch)
fit <- lm( formula=z~1+x+y, data=data.frame(x=predStokesV[lineRange], y=predStokesI[lineRange], z=StokesV13[lineRange]))
#plot( bunch_vec(freq[plotRange], plotBunch) - 0.5*plotBunch*chSep, bunch_vec(StokesV13[plotRange] - fit$coefficients[3]* predStokesI[plotRange] - fit$coefficients[1], plotBunch), type='s', ylim=c(-3e-2, 3e-2), xlab='Frequency [MHz]', ylab='Stokes V [K]', main='TMC-1 CCS')
plotX <- bunch_vec(freq[plotRange], plotBunch)
plotY <- bunch_vec(StokesV13[plotRange] - fit$coefficients[3]* predStokesI[plotRange] - fit$coefficients[1], plotBunch)
plot( plotX, plotY , pch=20, ylim=c(-2e-2, 2e-2), xlab='Frequency [MHz]', ylab='Stokes V [K]', main='TMC-1 CCS')	# plot points
arrows( plotX, plotY - err, plotX, plotY + err, angle=90, length=0)
#points( bunch_vec(freq[plotRange], plotBunch), bunch_vec(StokesV13[plotRange] - fit$coefficients[3]* predStokesI[plotRange] - fit$coefficients[1], plotBunch) , pch=20, cex=0.5)
lines(freq[plotRange]-0.5*chSep, predStokesV[plotRange]*fit$coefficients[2], type='s', col='red')
legend(min(freq[plotRange]), 0.03, legend=sprintf('Zeeman Shift = %5.1f ± %4.1f Hz', 10*fit[[1]][2], 10*summary(fit)[[4]][5]))
dev.off()
