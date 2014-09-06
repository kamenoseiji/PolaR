# WGScan.R
# usage: Rscript WGScan.R [list of prefix]
# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
#
library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readPolariS.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/date.R", ssl.verifypeer = FALSE)))
setwd('.')
#-------- Function to find Scan Pattern from SAM45 Log
scanXP <- function(prefix){
	PolarisFileNum <- length(prefix)
	C00 <- numeric(0); mjdSec <- numeric(0)
	for(index in 1:PolarisFileNum){
		cat( sprintf("Loading %s ….\n", prefix[index]))
		temp <- readPolariS_X(sprintf('%s.C.%02dB', prefix[index], 0))
		C00 <- append(C00, apply(Mod(apply(temp, 2, fft)), 2, max))
		mjdSec <- append(mjdSec, (prefix2MJDsec(prefix[index]) + seq(0, dim(temp)[2]-1, by=1)))
	}
	return(data.frame(mjdSec, C00))
}

#-------- Function to plot Amplitude and Phase in a page
amphi_plot <- function(Time, Vis, label){
	par.old <- par(no.readonly=TRUE)
	ampMax <- max(Mod(Vis))
	par(mfrow=c(2,1), oma=rep(0.5, 4), mar=c(0,4,4,4), xaxt='n')
	plot(Time, Mod(Vis), pch=20, cex=0.5, xlab='', ylab='Corr. Amplitude', main=label, ylim=c(0,ampMax))
	par(mar=c(4,4,0,4), xaxt='s')
	plot(Time, Arg(Vis), pch=20, cex=0.5, xlab='UT [hour]', ylab='Phase [rad]', ylim=c(-pi,pi))
	par(par.old)
	return()
}

#-------- Procedures
args <- commandArgs()
prefix <- args[6:length(args)]
XP <- scanXP(prefix)
save(XP, file=sprintf("%s.XP.Rdata", prefix[1]))
pdf(sprintf('%s.WG.pdf', prefix[1]))
amphi_plot( (XP$mjdSec%%86400)/3600, XP$C00, 'Wire Grid C00')
dev.off()
