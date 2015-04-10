# WGScan.R
# usage: Rscript WGScan.R [list of prefix]
# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
#
library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readPolariS.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/date.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/plotTool.R", ssl.verifypeer = FALSE)))
setwd('.')
#-------- Function to find Scan Pattern from SAM45 Log
scanXP <- function(prefix){
	PolarisFileNum <- length(prefix)
	C00 <- numeric(0); mjdSec <- numeric(0)
	for(index in 1:PolarisFileNum){
		cat( sprintf("Loading %s ….\n", prefix[index]))
		if( file.exists(sprintf('%s.C.%02dB', prefix[index], 0)) ){
			temp <- readPolariS_X(sprintf('%s.C.%02dB', prefix[index], 0))
		} else {
			temp <- readPolariS_X(sprintf('%s.C.%02d', prefix[index], 0))
		}
		C00 <- append(C00, apply(Mod(apply(temp, 2, fft)), 2, max))
		mjdSec <- append(mjdSec, (prefix2MJDsec(prefix[index]) + seq(0, dim(temp)[2]-1, by=1)))
	}
	return(data.frame(mjdSec, C00))
}


#-------- Procedures
args <- commandArgs(trailingOnly = T)
prefix <- args[1:length(args)]
XP <- scanXP(prefix)
save(XP, file=sprintf("%s.XP.Rdata", prefix[1]))
pdf(sprintf('%s.WG.pdf', prefix[1]))
time_amphi_plot( (XP$mjdSec%%86400)/3600, XP$C00, 'Wire Grid C00')
dev.off()
