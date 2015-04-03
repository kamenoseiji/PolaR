#-------- Function to plot Amplitude and Phase in a page
time_amphi_plot <- function(Time, Vis, label){
	par.old <- par(no.readonly=TRUE)
	ampMax <- max(Mod(Vis))
	par(mfrow=c(2,1), oma=rep(0.5, 4), mar=c(0,4,4,4), xaxt='n')
	plot(Time, Mod(Vis), pch=20, cex=0.5, xlab='', ylab='Corr. Amplitude', main=label, ylim=c(0,ampMax))
	par(mar=c(4,4,0,4), xaxt='s')
	plot(Time, Arg(Vis), pch=20, cex=0.5, xlab='UT [hour]', ylab='Phase [rad]', ylim=c(-pi,pi))
	par(par.old)
}
#-------- Function to plot Two values in a page
time_two_plot <- function(Time, value1, value2, label){
	par.old <- par(no.readonly=TRUE)
	Max1 <- max(value1); Max2 <- max(value2); Min1 <- min(value1); Min2 <- min(value2)
	par(mfrow=c(2,1), oma=rep(0.5, 4), mar=c(0,4,4,4), xaxt='n')
	plot(Time, value1, pch=20, cex=0.5, xlab='', ylab=label$value1, main=label$title, ylim=c(Min1,Max1))
	par(mar=c(4,4,0,4), xaxt='s')
	plot(Time, value2, pch=20, cex=0.5, xlab=label$time, ylab=label$value2, ylim=c(Min2, Max2))
	par(par.old)
}

#-------- Function to plot Amplitude and Phase of cross power spectrum in a page
xspec_amphi_plot <- function(freq, xspec, label){
	par.old <- par(no.readonly=TRUE)
	ampMax <- max(Mod(xspec))
	par(mfrow=c(2,1), oma=rep(0.5, 4), mar=c(0,4,4,4), xaxt='n')
	plot(freq, Mod(xspec), type='l', xlab='', ylab=label$amp, main=label$title, ylim=c(0,ampMax))
	par(mar=c(4,4,0,4), xaxt='s')
	plot(freq, Arg(xspec), xlab=label$freq, ylab=label$phase, ylim=c(-pi,pi), type='n')
	abline(h=0, col='red', lwd=0.5)
	points(freq, Arg(xspec), pch=20, cex=0.5)
	par(par.old)
}

plotTsys <- function( Time, Tsys, OnIndex, OfIndex, label){
	plot( Time[c(OnIndex, OfIndex)], Tsys[c(OnIndex, OfIndex)], xlab=label$Time, ylab=label$Tsys, main=label$Title, type='n')
	points( Time[OnIndex], Tsys[OnIndex], pch=20, cex=0.3, col='green')
	points( Time[OfIndex], Tsys[OfIndex], pch=20, cex=0.3, col='blue')
}
