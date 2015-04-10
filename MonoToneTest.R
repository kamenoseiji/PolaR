library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readPolariS.R", ssl.verifypeer = FALSE)))
setwd('/Volumes/SSD/PolariS/20150317/')
C00 <- readPolariS_X('2015076095301.C.00')
C01 <- readPolariS_X('2015076095301.C.01')

chNum <- nrow(C00)
freq <- (0:(chNum-1))*4e6/chNum
C00_PeakCH <- 29214
C01_PeakCH <- 28952

C00_timeIndex <- which( abs(C00[C00_PeakCH,]) > 8)
C01_timeIndex <- 220:250

Xspec00 <- apply(C00[,C00_timeIndex], 1, mean)
Xspec01 <- apply(C01[,C01_timeIndex], 1, mean)

C00_plotRange <- (C00_PeakCH-64):(C00_PeakCH+64)
C01_plotRange <- (C01_PeakCH-64):(C01_PeakCH+64)

pdf('MonoTone.pdf')
plot( freq[C00_plotRange], Mod(Xspec00[C00_plotRange]), xlab='Frequency [Hz]', ylab='Cross Power Spectrum', main='SG@45.490316 GHz', type='n')
lines( freq[C00_plotRange], abs(Re(Xspec00[C00_plotRange])), type='s', col='blue')
lines( freq[C00_plotRange], abs(Im(Xspec00[C00_plotRange])), type='s', col='red')

plot( freq[C01_plotRange], Mod(Xspec01[C01_plotRange]), xlab='Frequency [Hz]', ylab='Cross Power Spectrum', main='SG@45.379000 GHz', type='n')
lines( freq[C01_plotRange], abs(Re(Xspec01[C01_plotRange])), type='s', col='blue')
lines( freq[C01_plotRange], abs(Im(Xspec01[C01_plotRange])), type='s', col='red')
dev.off()
