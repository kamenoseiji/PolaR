# ScanPattern.R
# usage: Rscript ScanPattern.R [Number of SAM45 file] [list of SAM45 file]
# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
#
library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readSAM45.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readPolariS.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/Qeff.R", ssl.verifypeer = FALSE)))
setwd('.')
#-------- Function to calculate Tsys from Scan Pattern
scanTsys <- function(Scan, Tamb){
	R_index <- which(Scan$scanType == 'R')
	index <- which(Scan$scanType == 'OFF' | Scan$scanType == 'ON' | Scan$scanType == 'SKY')
	power_ptr <- grep('power', names(Scan)); IFnum <- length(power_ptr)
	for(IF_index in 1:IFnum){
		Tsys   <- rep(NA, length(Scan$mjdSec))
		RPower <- predict(smooth.spline(Scan$mjdSec[R_index], Scan[[power_ptr[IF_index]]][R_index], spar=1.0), Scan$mjdSec)$y
		Tsys[index]  <- Tamb / (RPower / Scan[[power_ptr[IF_index]]][index] - 1.0)
		Scan <- cbind(Scan, Tsys)
	}
	return(Scan)
}
		
	
	#RPower <- predict(smooth.spline(Scan$mjdSec[R_index], Scan$power[R_index], spar=1.0), Scan$mjdSec)$y
	
	
	#Range <- range( Scan$mjdSec[c(Scan$on, Scan$off, Scan$R)]); mjdRange <- Range[1]:Range[2]
	#RPower <- predict(smooth.spline(Scan$mjdSec[Scan$R], Scan$power[Scan$R], spar=1.0), mjdRange)$y
	#offPower <- predict(smooth.spline(Scan$mjdSec[Scan$off], Scan$power[Scan$off], spar=0.8), mjdRange)$y
	#onPower <- predict(smooth.spline(Scan$mjdSec[Scan$on], Scan$power[Scan$on], spar=0.25), mjdRange)$y
	#on_Tsys  <- Tamb / (RPower / onPower - 1.0)
	#off_Tsys <- Tamb / (RPower / offPower - 1.0)
	#return( data.frame(mjdSec = mjdRange, TsysOn = on_Tsys, TsysOff = off_Tsys))
#}

#-------- Procedures
#args <- commandArgs()
prefix <- character(0)
#SAM45File <- args[6:length(args)]
SAM45File <- "SAM45.TMC1.np32802.proj5.20140417103636"

#-------- List prefix of PolariS data
Year <- substr(strsplit(SAM45File[1], '\\.')[[1]][5], 1, 4)
P00fileList <- system(  sprintf('ls %s*.P.00', Year), intern=T )
for(index in 1:length(P00fileList)){
	prefix[index] <- substr(P00fileList[index], 1, 13)
}

#-------- List IF ID of PolariS data
P00fileList <- system(  sprintf('ls %s.P.*', prefix[1]), intern=T )
IF_ID <- integer(0)
for(index in 1:length(P00fileList)){
	IF_ID[index] <- as.integer(strsplit(P00fileList[index] , '\\.')[[1]][3])
}

#-------- Produce Scan data frame
for(fileIndex in 1:length(SAM45File)){
	tempScan <- scanPattern(SAM45File[fileIndex], prefix, IF_ID)
	if(fileIndex == 1){	Scan <- tempScan}
	else { Scan <- rbind(Scan, tempScan)}
}

#-------- Tsys
Scan <- scanTsys(Scan, 290.0)