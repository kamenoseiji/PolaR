# WGScan.R
# usage: Rscript WGScan.R [list of prefix]
# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
#
library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readSAM45.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readPolariS.R", ssl.verifypeer = FALSE)))
setwd('.')
#-------- Function to calculate Tsys from Scan Pattern
scanTsys <- function(Scan, Tamb){
	Range <- range( Scan$mjdSec[c(Scan$on, Scan$off, Scan$R)]); mjdRange <- Range[1]:Range[2]
	RPower <- predict(smooth.spline(Scan$mjdSec[Scan$R], Scan$power[Scan$R], spar=1.0), mjdRange)$y
	offPower <- predict(smooth.spline(Scan$mjdSec[Scan$off], Scan$power[Scan$off], spar=0.8), mjdRange)$y
	onPower <- predict(smooth.spline(Scan$mjdSec[Scan$on], Scan$power[Scan$on], spar=0.25), mjdRange)$y
	on_Tsys  <- Tamb / (RPower / onPower - 1.0)
	off_Tsys <- Tamb / (RPower / offPower - 1.0)
	return( data.frame(mjdSec = mjdRange, TsysOn = on_Tsys, TsysOff = off_Tsys))
}

#-------- Procedures
args <- commandArgs()
#SAM45File <- args[6:length(args)]
#SAM45File   <- c('SAM45.TMC1.np32802.proj5.20140417103636', 'SAM45.TMC1.np32802.proj5.20140417115456', 'SAM45.TMC1.np32802.proj5.20140417130456', 'SAM45.TMC1.np32802.proj5.20140417140954', 'SAM45.TMC1.np32802.proj5.20140417150105', 'SAM45.TMC1.np32802.proj5.20140417160609')

#prefix <- list(c('2014107013610', '2014107020610', '2014107023610'), c('2014107023610', '2014107030610', '2014107033610'), c('2014107033610', '2014107040610', '2014107043610'), c('2014107050610', '2014107053610'),  c('2014107053610', '2014107060610', '2014107063610'), c('2014107070610', '2014107073610', '2014107080610'))

SAM45File <- 'SAM45.TMC1.np32802.proj5.20140417103636'
prefix <- c('2014107013610', '2014107020610', '2014107023610')

for(IF_index in 0:3){
	for(fileIndex in 1:length(SAM45File)){
		Scan <- scanPattern(SAM45File[fileIndex], prefix[[fileIndex]], IF_index)
		temp <- scanTsys(Scan)
		if( fileIndex == 1 ){
			TsysDF <- temp
		} else {
			TsysDF <- rbind(TsysDF, temp)
		}
	}
	assign( paste("TsysDF", IF_index, sep=""), TsysDF )
}
save(TsysDF0, TsysDF1, TsysDF2, TsysDF3, file='2014107010610.Tsys.Rdata')
