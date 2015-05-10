#-------- Procedures
args <- commandArgs(trailingOnly = T)
StokeFile <- args[1:length(args)]
setwd('.')

#-------- Load Stokes Files
for(fileIndex in 1:length(StokeFile)){
	load( StokesFile[fileIndex] )
	if( fileIndex == 1 ){
		temp00 <- uncalStokes00; temp01 <- uncalStokes01
	} else {
		temp00 <- rbind(temp00, uncalStokes00); temp01 <- rbind(temp01, uncalStokes01)
	}
}
uncalStokes00 <- temp00; uncalStokes01 <- temp01
#-------- Plot Stokes Parameters
pdf('Stokes.pdf')
plot(uncalStokes00$mjdSec, uncalStokes00$I, pch=20, ylim=c(-1, max(uncalStokes00$I)), xlab='MJD [sec]', ylab='Stokes Parameters [K]'); abline(h=mean(uncalStokes00$I))
points(uncalStokes00$mjdSec, uncalStokes00$Q, pch=20, col=2); abline(h=mean(uncalStokes00$Q), col=2)
points(uncalStokes00$mjdSec, uncalStokes00$U, pch=20, col=3); abline(h=mean(uncalStokes00$U), col=3)
points(uncalStokes00$mjdSec, uncalStokes00$V, pch=20, col=4); abline(h=mean(uncalStokes00$V), col=4)
#
plot(uncalStokes00$mjdSec, uncalStokes00$EVPA, pch=20, ylim=c(-90, 90), xlab='MJD [sec]', ylab='EVPA [deg]'); abline(h=152-180)
dev.off()