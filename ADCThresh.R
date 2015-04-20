# ADCThresh.R
# usage: Rscript ADCThresh.R prefix
# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
#
library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readSAM45.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readPolariS.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/Qeff.R", ssl.verifypeer = FALSE)))
setwd('.')

#-------- Procedures
args <- commandArgs(trailingOnly = T)
#args <- c('2015110025219')
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
	bitDist <- readBitDist(sprintf("%s.P.%02d", prefix, IF_ID[if_index]))
	thresh <- apply(bitDist, 2, threshLevel)
	avail_index <- which(colSums(is.infinite(thresh)) + colSums(is.na(thresh)) == 0)
	if(length(avail_index) == 0){ cat(sprintf("Warning: No available data in %s.P.%02d\n", prefix, IF_ID[if_index]))}
	Thresh[,if_index] <- apply(thresh[,avail_index], 1, mean)
}
fileName <- sprintf("%s.Thresh.Rdata", prefix)
save(Thresh, file=fileName)
cat( sprintf("Thresholds are saved in %s.\n", fileName))
