# ScanPattern.R
# usage: Rscript ScanPattern.R [SAM45 directory] [prefix]
# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
#
library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readPolariS.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readSAM45.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/Qeff.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/findSAM45Log.R", ssl.verifypeer = FALSE)))
#
args <- commandArgs()
SAM45LOG_dir <- args[6]
prefix <- args[7]
setwd(SAM45LOG_dir)
SAM45logFile <- paste(SAM45LOG_dir, findSAM45Log(SAM45LOG_dir, prefix), sep="")
T_hot <- 290.0			# Hot load temperature [K]
#
#-------- Read SAM45 Logging
SAM45Log <- readSAM45(SAM45logFile)
head1 <- SAM45Log[[1]]; head2 <- SAM45Log[[2]]; SAM45spec <- SAM45Log[[3]]; SAM45df <- SAM45Log[[4]]
#-------- Read PolariS total power data
bitPowerIF <- list(); IFNum <- 4
for(IF_index in 1:IFNum){
	bitDist <- apply(readBitDist( sprintf('%s.P.%02d', prefix, IF_index-1), 256), 2, bunchVec16)
	gaussResults <- apply(bitDist, 2, gauss4bit)
	bitPowerIF[[IF_index]] <- 1/gaussResults[1,]^2
}
bitPower <- data.frame(IF0=bitPowerIF[[1]], IF1=bitPowerIF[[2]], IF2=bitPowerIF[[3]], IF3=bitPowerIF[[4]])
#
#-------- PolariS time tags
mjdSec <- prefix2MJDsec(prefix) + seq(0, length(gaussResults[1,])-1, by=1)
#-------- On/Off/R scan pattern
onScan <- c(); offScan <- c(); RScan <- c();
arrayIndex <- which(SAM45df$cary_name == SAM45df$cary_name[1])
for(scanIndex in arrayIndex){
	polaris_index <- which((mjdSec >= (SAM45df$mjd_st[scanIndex]-1)) & (mjdSec <= SAM45df$mjd_ed[scanIndex]))
	if( (SAM45df$dAZ[scanIndex]^2 + SAM45df$dEL[scanIndex]^2 < 1.0) & (SAM45df$cscan_type[scanIndex] == 'ON') )	onScan <- append(onScan, polaris_index)	# On scan
	if( SAM45df$cscan_type[scanIndex] == 'OFF')	offScan <- append(offScan, (min(polaris_index)-5):(max(polaris_index)+1))								# Off scan
	if( SAM45df$cscan_type[scanIndex] == 'R' )	RScan <- append(RScan, tail(polaris_index, length(polaris_index)-1))									# Hot load
}
#-------- Tsys
TsysIF <- list()
for(IF_index in 1:IFNum){ TsysIF[[IF_index]] <- T_hot / ( mean(bitPowerIF[[IF_index]][RScan]) / bitPowerIF[[IF_index]] - 1.0); TsysIF[[IF_index]][RScan] <- NA }
# Tsys00 <- T_hot / ( mean(bitPower00[RScan]) / bitPower00 - 1.0); Tsys00[RScan] <- NA
# Tsys01 <- T_hot / ( mean(bitPower01[RScan]) / bitPower01 - 1.0); Tsys01[RScan] <- NA
# Tsys02 <- T_hot / ( mean(bitPower02[RScan]) / bitPower02 - 1.0); Tsys02[RScan] <- NA
# Tsys03 <- T_hot / ( mean(bitPower03[RScan]) / bitPower03 - 1.0); Tsys03[RScan] <- NA
Tsys <- data.frame(IF0=TsysIF[[1]], IF1=TsysIF[[2]], IF2=TsysIF[[3]], IF3=TsysIF[[4]])
#-------- OutPuts
save(head1, head2, SAM45df, mjdSec, onScan, offScan, RScan, SAM45spec, bitPower, Tsys, file=paste(prefix, ".scan.Rdata", sep=""))
