#-------- Parse command-line arguments
parseArg <- function( args ){
    argNum <- length(args)
    fileNum <- argNum
    for( index in 1:argNum ){
        if(substr(args[index], 1,2) == "-S"){ srcName <- substring(args[index], 3); fileNum <- fileNum - 1}
        if(substr(args[index], 1,2) == "-E"){ EVPA <- as.numeric(substring(args[index], 3)); fileNum <- fileNum - 1}
    }
    fileName <- args[(argNum - fileNum + 1):argNum]
    return( list(srcName = srcName, EVPA = EVPA, fileName = fileName) )
}
#
PlotStokes <- function(mjdSec, Value, label, ColIndex){
    points(mjdSec, Value, pch=20, cex=0.5, col=ColIndex)
    abline(h=mean(Value) , col=ColIndex)
    sd_text <- sprintf('%s = %5.2f ± %5.3f K', label, mean(Value), sd(Value) / sqrt(length(Value) - 1))
    text( min(mjdSec), mean(Value)+0.1, sd_text, pos=4, offset=1.5, cex=0.7, col=ColIndex)
    return()
}

#-------- Procedures
args <- parseArg(commandArgs(trailingOnly = T))
StokesFile <- args$fileName
cat(args$srcName); cat('\n')
cat(args$EVPA); cat('\n')
cat(StokesFile); cat('\n')
setwd('.')

#-------- Load Stokes Files
for(fileIndex in 1:length(StokesFile)){
	load( StokesFile[fileIndex] )
	if( fileIndex == 1 ){
		temp00 <- Stokes00; temp01 <- Stokes01
	} else {
		temp00 <- rbind(temp00, Stokes00); temp01 <- rbind(temp01, Stokes01)
	}
}
Stokes00 <- temp00; Stokes01 <- temp01
#-------- Plot Stokes Parameters
pdf(sprintf('%s.StokesPlot.pdf', strsplit(args$fileName, "\\.")[[1]][1]))
TaRange <- range( c(Stokes00[[1]], Stokes00[[2]], Stokes00[[3]], Stokes00[[4]]) )
plot(Stokes00$mjdSec, Stokes00[[1]], pch=20, ylim=TaRange, xlab='MJD [sec]', ylab='Stokes Parameters [K]', main=args$srcName, type='n')
Label <- c('I', 'Q', 'U', 'V')
#-------- For each Stokes Parameters
for(StokesIndex in 1:4){
    PlotStokes( c(Stokes00$mjdSec, Stokes01$mjdSec), c(Stokes00[[StokesIndex]], Stokes01[[StokesIndex]]), Label[StokesIndex], StokesIndex)
}
#
mjdSec <- c(Stokes00$mjdSec, Stokes00$mjdSec); Value <- c(Stokes00$EVPA, Stokes01$EVPA)
plot(mjdSec, Value, pch=20, ylim=c(-90, 90), xlab='MJD [sec]', ylab='EVPA [deg]', main=args$srcName); abline(h=args$EVPA-180)
sd_text <- sprintf('EVPA = %5.1f ± %5.2f deg', mean(Value), sd(Value) / sqrt(length(Value) - 1))
text( min(mjdSec), args$EVPA-180+2, sd_text, pos=4, offset=1.5, cex=0.7)
dev.off()
