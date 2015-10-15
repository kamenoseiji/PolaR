fileList <- commandArgs(trailingOnly = T)
fileNum <- length(fileList)
load(fileList[1]); DF <- SquintDF
for(file_index in 2:fileNum){
	load(fileList[file_index])
	DF <- rbind(DF, SquintDF)
}

SquintDF <- data.frame(
    dAZ = cospi(DF$EL/180.0)*(DF$Raz - DF$Laz) - sinpi(DF$EL/180.0)*(DF$Rel - DF$Lel),
    dEL = sinpi(DF$EL/180.0)*(DF$Raz - DF$Laz) + cospi(DF$EL/180.0)*(DF$Rel - DF$Lel),
    eAZ = sqrt( DF$Raze^2 + DF$Laze^2),
    eEL = sqrt( DF$Rele^2 + DF$Lele^2))
Waz <- 1.0/SquintDF$eAZ^2
Wel <- 1.0/SquintDF$eEL^2
Maz <- SquintDF$dAZ %*% Waz / sum(Waz); Eaz <- sqrt(1.0/sum(Waz))
Mel <- SquintDF$dEL %*% Wel / sum(Wel); Eel <- sqrt(1.0/sum(Wel))
text_sd <- sprintf("dAZ = %5.2f ± %5.2f  dEL = %5.2f ± %5.2f",  Maz, Eaz, Mel, Eel )
cat(text_sd); cat('\n')
pdf("Squint.pdf")
#Xlim <- c( min(SquintDF$dAZ - SquintDF$eAZ), max(SquintDF$dAZ + SquintDF$eAZ) )
#Ylim <- c( min(SquintDF$dEL - SquintDF$eEL), max(SquintDF$dEL + SquintDF$eEL) )
Xlim <- c(-5,5); Ylim <- c(-5,5)
plot(SquintDF$dAZ, SquintDF$dEL, asp=1, xlim=Xlim, ylim=Ylim, type='n', xlab='AZ [arcsec]', ylab='EL [arcsec]')
abline(v=0, lwd=0.25); abline(h=0, lwd=0.25)
arrows( SquintDF$dAZ, SquintDF$dEL - SquintDF$eEL, SquintDF$dAZ, SquintDF$dEL + SquintDF$eEL, length=0, col=gray(1.0 - Waz/max(Waz)))
arrows( SquintDF$dAZ - SquintDF$eAZ, SquintDF$dEL, SquintDF$dAZ + SquintDF$eAZ, SquintDF$dE, length=0, col=gray(1.0 - Waz/max(Waz)))
points(SquintDF$dAZ, SquintDF$dEL, pch=20, cex=Waz/max(Waz))
text( Xlim[1] - 1, Ylim[2] - 1, pos=4, text_sd)
dev.off()
