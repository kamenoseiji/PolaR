# ScanPattern.R
# usage: Rscript ScanPattern.R [Number of SAM45 file] [list of SAM45 file]
# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
#
library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readSAM45.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readPolariS.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/Qeff.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/date.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/plotTool.R", ssl.verifypeer = FALSE)))
setwd('.')
#-------- Function to calculate Tsys from Scan Pattern
scanTsys <- function(Scan, Tamb){
	nameList <- names(Scan)
	R_index <- which(Scan$scanType == 'R')
	#-------- Median Window Filter
	medR <- median( Scan$power00[R_index])
	R_index <- which( abs( (Scan$power00 - medR)/medR ) < 0.02 )
	#-------- On and Off scans
	index <- which(Scan$scanType == 'OFF' | Scan$scanType == 'ON')
	OutOfR_index <- which( Scan$mjdSec > max(Scan$mjdSec[R_index]))
	power_ptr <- grep('power', nameList); IFnum <- length(power_ptr)
	for(IF_index in 1:IFnum){
		 IF_ID <- as.integer(strsplit(nameList[power_ptr[IF_index]], "power")[[1]][2])
		 Tsys   <- rep(NA, length(Scan$mjdSec))
		 RPower <- predict(smooth.spline(Scan$mjdSec[R_index], Scan[[power_ptr[IF_index]]][R_index], spar=1.0), Scan$mjdSec)$y
		 RPower[OutOfR_index] <- RPower[max(R_index)]
		 Tsys[index]  <- Tamb / (RPower[index] / Scan[[power_ptr[IF_index]]][index] - 1.0)
		 Scan <- cbind(Scan, Tsys)
		 nameList <- append(nameList, sprintf('Tsys%02d', IF_ID))
	}
	names(Scan) <- nameList
	return(Scan)
}

#-------- Procedures
args <- commandArgs(trailingOnly = T)
#args <- c('2015114041130.Thresh.Rdata', 'SAM45.S235.as708fn.proj2.20150424110741', 'SAM45.S235.as708fn.proj2.20150424112534')
#setwd('/Users/kameno/SVDATA/20150424')
prefix <- character(0)
threshFile <- args[1]
SAM45File <- args[2:length(args)]

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
	cat(sprintf('Processing %s ...\n', SAM45File[fileIndex]))
	tempScan <- scanPattern(SAM45File[fileIndex], prefix, IF_ID, threshFile)
	if(fileIndex == 1){	Scan <- tempScan}
	else { Scan <- rbind(Scan, tempScan)}
}
#-------- Tsys
Scan <- scanTsys(Scan, 280.0)
#-------- Save
StartUTC <- mjd2doy(Scan$mjdSec[1])
filePrefix <- sprintf("%04d%03d%02d%02d%02d", StartUTC$year, StartUTC$doy, StartUTC$hour, StartUTC$min, StartUTC$sec)
fileName <- sprintf("%s.Scan.Rdata", filePrefix)
save(Scan, file=fileName)
cat('Scan and Tsys records are saved into '); cat(fileName); cat('\n')
#-------- Plot
OnIndex <- which( Scan$scanType == 'ON')
OfIndex <- which( Scan$scanType == 'OFF')
pdf(sprintf('%s.tsys.pdf', filePrefix))
plotTsys( (Scan$mjdSec%%86400)/3600, Scan$Tsys00, OnIndex, OfIndex, list(Time='UT [hour]', Tsys='Tsys [K]', Title=sprintf('%s %s IF=%d', SAM45File[1], prefix[1], 0)))
plotTsys( (Scan$mjdSec%%86400)/3600, Scan$Tsys01, OnIndex, OfIndex, list(Time='UT [hour]', Tsys='Tsys [K]', Title=sprintf('%s %s IF=%d', SAM45File[1], prefix[1], 1)))
plotTsys( (Scan$mjdSec%%86400)/3600, Scan$Tsys02, OnIndex, OfIndex, list(Time='UT [hour]', Tsys='Tsys [K]', Title=sprintf('%s %s IF=%d', SAM45File[1], prefix[1], 2)))
plotTsys( (Scan$mjdSec%%86400)/3600, Scan$Tsys03, OnIndex, OfIndex, list(Time='UT [hour]', Tsys='Tsys [K]', Title=sprintf('%s %s IF=%d', SAM45File[1], prefix[1], 3)))
dev.off()
