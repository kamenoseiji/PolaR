# ADCThresh.R
# usage: Rscript ADCThresh.R prefix
# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
#
RPATH <- '~/Programs/PolaR'
FuncList <- c('readSAM45', 'readPolariS', 'Qeff')
source(sprintf('%s/loadModule.R', RPATH))
library(RCurl)

funcNum <- length(FuncList)
for( index in 1:funcNum){
    URL <- sprintf("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/%s.R", FuncList[index])
    Err <- try( eval(parse(text = getURL(URL, ssl.verifypeer = FALSE))), silent=FALSE)
}   
if(class(Err) == "try-error"){ loadLocal( RPATH, FuncList ) }

setwd('.')

#-------- Procedures
args <- commandArgs(trailingOnly = T)
#args <- c('2015110025209')
#setwd('/Volumes/SSD/PolariS/20150420/')
prefix <- args[1]

#-------- List IF ID of BitPower
P00fileList <- system(  sprintf('ls %s.P.*', prefix), intern=T )
IF_ID <- integer(0)
for(index in 1:length(P00fileList)){
	IF_ID[index] <- as.integer(strsplit(P00fileList[index] , '\\.')[[1]][3])
}

LevelNum <- 256
ThreshNum <- LevelNum - 1
Thresh <- matrix( nrow=ThreshNum, ncol=length(IF_ID))
for(if_index in 1:length(IF_ID)){
    readFile <- sprintf("%s.P.%02d", prefix, IF_ID[if_index])
	levelHisto <- matrix(readBitDist(readFile), nrow=LevelNum)
	avail_index <- which( levelHisto[1,] > 0.9* max(levelHisto[1,]))
	threshHolds <- apply(levelHisto[,avail_index], 2, threshLevel)
	avail_index <- which((colSums(is.infinite(threshHolds)) + colSums(is.na(threshHolds))) == 0)
	if(length(avail_index) == 0){ cat(sprintf("Warning: No available data in %s.P.%02d\n", prefix, IF_ID[if_index]))}
	Thresh[,if_index] <- apply(threshHolds[,avail_index], 1, mean)
}
fileName <- sprintf("%s.Thresh.Rdata", prefix)
save(Thresh, file=fileName)
cat( sprintf("Thresholds are saved in %s.\n", fileName))
