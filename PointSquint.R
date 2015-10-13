# ScanPattern.R
# usage: Rscript ScanPattern.R [Number of SAM45 file] [list of SAM45 file]
# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
#
RPATH <- '~/Programs/PolaR'
FuncList <- c('readSAM45', 'readPolariS', 'PolariCalib', 'Qeff', 'date', 'plotTool')
source(sprintf('%s/loadModule.R', RPATH))
library(RCurl)

funcNum <- length(FuncList)
for( index in 1:funcNum){
    URL <- sprintf("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/%s.R", FuncList[index])
    Err <- try( eval(parse(text = getURL(URL, ssl.verifypeer = FALSE))), silent=FALSE)
}   
if(class(Err) == "try-error"){ loadLocal( RPATH, FuncList ) }

setwd('.')

SAM45File <- '/Volumes/HDD/PolariS/SAM45/SAM45.OMC2p5.as708fn.proj2.20150424124554'
#-------- List prefix of PolariS data
prefix <- character(0)
Year <- substr(strsplit(SAM45File, '\\.')[[1]][5], 1, 4)
A00fileList <- system(  sprintf('ls %s*.A.00', Year), intern=T )
ipnum <- integer(length(A00fileList))
for(index in 1:length(A00fileList)){
	prefix[index] <- substr(A00fileList[index], 1, 13)
	AfileSize <- GetChnumRecnum(sprintf('%s.A.00', prefix[index]), 'A')
	ipnum[index] <- AfileSize$ipnum
}
chNum <- AfileSize$chnum
#-------- List IF ID of PolariS data
A00fileList <- system(  sprintf('ls %s.A.*', prefix[1]), intern=T )
IF_ID <- integer(0)
for(index in 1:length(A00fileList)){
	IF_ID[index] <- as.integer(strsplit(A00fileList[index] , '\\.')[[1]][3])
}
#-------- Produce Scan data frame
Scan <- scanPointing(SAM45File)
OnIndex <- which( Scan$scanType == 'ON');  onMJD  <- scanSegment(Scan$mjdSec[OnIndex])
OfIndex <- which( Scan$scanType == 'OFF'); offMJD <- scanSegment(Scan$mjdSec[OfIndex])

#AfileIndex <- seq(findPrefix(min(SAM45df$mjd_st), prefix), findPrefix(max(SAM45df$mjd_st), prefix))
on_C00  <- integSegment(prefix, chNum, ipnum, 'C', 0, onMJD )
on_C01  <- integSegment(prefix, chNum, ipnum, 'C', 1, onMJD )
on_A00  <- integSegment(prefix, chNum, ipnum, 'A', 0, onMJD )
on_A01  <- integSegment(prefix, chNum, ipnum, 'A', 1, onMJD )
on_A02  <- integSegment(prefix, chNum, ipnum, 'A', 2, onMJD )
on_A03  <- integSegment(prefix, chNum, ipnum, 'A', 3, onMJD )
off_C00 <- integSegment(prefix, chNum, ipnum, 'C', 0, offMJD )
off_C01 <- integSegment(prefix, chNum, ipnum, 'C', 1, offMJD )
off_A00 <- integSegment(prefix, chNum, ipnum, 'A', 0, offMJD )
off_A01 <- integSegment(prefix, chNum, ipnum, 'A', 1, offMJD )
off_A02 <- integSegment(prefix, chNum, ipnum, 'A', 2, offMJD )
off_A03 <- integSegment(prefix, chNum, ipnum, 'A', 3, offMJD )

#pointingNum <- floor((length(offMJD$stopMjd) - 1) / 2)
#for index in 1:pointingNum{
#}


if(0){
#-------- Procedures
args <- commandArgs(trailingOnly = T)
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
Scan <- data.frame()
for(fileIndex in 1:length(SAM45File)){
	cat(sprintf('Processing %s ...\n', SAM45File[fileIndex]))
	tempScan <- scanPattern(SAM45File[fileIndex], prefix, IF_ID, threshFile)
    Scan <- rbind(Scan, tempScan$scanDF)
	# if(fileIndex == 1){	Scan <- tempScan$scanDF}
	# else { Scan <- rbind(Scan, tempScan$scanDF)}
}
Head <- tempScan$head2
#-------- Tsys
Scan <- scanTsys(Scan, 280.0)
#-------- Save
StartUTC <- mjd2doy(Scan$mjdSec[1])
filePrefix <- sprintf("%04d%03d%02d%02d%02d", StartUTC$year, StartUTC$doy, StartUTC$hour, StartUTC$min, StartUTC$sec)
fileName <- sprintf("%s.Scan.Rdata", filePrefix)
save(Scan, Head, file=fileName)
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
}
