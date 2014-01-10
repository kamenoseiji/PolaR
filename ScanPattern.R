# ScanPattern.R
# usage: Rscript ScanPattern.R [prefix]
# prefix is YYYYDOYHHMMSS in the PolariS file name (e.g. 2013362105803)
#
SAM45LOG_dir <- "/Volumes/SSD/PolariS/SAM45/"
Rscript_dir  <- "/Volumes/SSD/PolariS/R/"
source(paste(Rscript_dir, "readPolariS.R", sep=""))
source(paste(Rscript_dir, "readSAM45.R", sep=""))
source(paste(Rscript_dir, "Qeff4bit.R", sep=""))
source(paste(Rscript_dir, "findSAM45Log.R", sep=""))
#
args <- commandArgs()
prefix <- args[6]
SAM45logFile <- paste(SAM45LOG_dir, findSAM45Log(SAM45LOG_dir, prefix), sep="")
T_hot <- 290.0			# Hot load temperature [K]
#
#-------- Read SAM45 Logging
SAM45Log <- readSAM45(SAM45logFile)
head1 <- SAM45Log[[1]]; head2 <- SAM45Log[[2]]; SAM45spec <- SAM45Log[[3]]; SAM45df <- SAM45Log[[4]]
#-------- Read PolariS total power data
bitDist00 <- readBitDist(paste(prefix, '.P.00', sep=''))
bitDist01 <- readBitDist(paste(prefix, '.P.01', sep=''))
bitDist02 <- readBitDist(paste(prefix, '.P.02', sep=''))
bitDist03 <- readBitDist(paste(prefix, '.P.03', sep=''))
#
gaussResults <- apply(bitDist00, 2, gauss4bit); bitPower00 <- 1/gaussResults[1,]^2
gaussResults <- apply(bitDist01, 2, gauss4bit); bitPower01 <- 1/gaussResults[1,]^2
gaussResults <- apply(bitDist02, 2, gauss4bit); bitPower02 <- 1/gaussResults[1,]^2
gaussResults <- apply(bitDist03, 2, gauss4bit); bitPower03 <- 1/gaussResults[1,]^2
bitPower <- data.frame(IF0=bitPower00, IF1=bitPower01, IF2=bitPower02, IF3=bitPower03)
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
Tsys00 <- T_hot / ( mean(bitPower00[RScan]) / bitPower00 - 1.0); Tsys00[RScan] <- NA
Tsys01 <- T_hot / ( mean(bitPower01[RScan]) / bitPower01 - 1.0); Tsys01[RScan] <- NA
Tsys02 <- T_hot / ( mean(bitPower02[RScan]) / bitPower02 - 1.0); Tsys02[RScan] <- NA
Tsys03 <- T_hot / ( mean(bitPower03[RScan]) / bitPower03 - 1.0); Tsys03[RScan] <- NA
Tsys <- data.frame(IF0=Tsys00, IF1=Tsys01, IF2=Tsys02, IF3=Tsys03)
#-------- OutPuts
save(head1, head2, SAM45df, mjdSec, onScan, offScan, RScan, SAM45spec, bitPower, Tsys, file=paste(prefix, ".scan", sep=""))
