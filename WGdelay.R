# WGdelay.R : Generate delay and phase table using WG scans
# usage: Rscript WGdelay.R [prefix list]
# e.g. Rscript WGdelay.R 2014107010610 2014107013610 2014107020610 2014107023610 2014107030610 2014107033610 2014107040610 2014107043610 2014107050610 2014107053610 2014107060610 2014107063610 2014107070610 2014107073610 2014107080610 2014107083610 2014107090610
 
#
library(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readPolariS.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/date.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/mjd.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/plotTool.R", ssl.verifypeer = FALSE)))
setwd('.')

#-------- Procedures
args <- commandArgs(trailingOnly = T)
prefix <- args[1:length(args)]
XPfname <- sprintf('%s.XP.Rdata', prefix[1])
BPfname <- sprintf('%s.BP.Rdata', prefix[1])

#-------- Load XP, scanXP, and BP tables
load(XPfname)
load(BPfname)
postFix <- ''
if(file.exists(sprintf('%s.C.%02dB', prefix[1], 0))){ postFix <- 'B' }
chnum <- GetChNum(sprintf('%s.C.%02d%s', prefix[1], 0, postFix))
bunchNum <- length(BP$BP00) / chnum
BP <- data.frame(BP00 = bunch_vec(BP$BP00, bunchNum), BP01 = bunch_vec(BP$BP01, bunchNum))
chRange <- floor(chnum*0.05):floor(chnum*0.95)

#-------- Delay determination
pdf(sprintf('%s.delay.pdf', prefix[1]))
cat('UT        Delay  Amp     Phase    Delay  Amp     Phase\n')
delayC00 <- numeric(0); delayC01 <- numeric(0); C00Vis <- complex(0); C01Vis <- complex(0); mjdSec <- numeric(0)
for(index in 1:length(scanXP$startMJD)){
	startFileIndex <- findPrefix(scanXP$startMJD[index], prefix)
	endFileIndex   <- findPrefix(scanXP$endMJD[index], prefix)
	mjdSec[index] <- mean(c(scanXP$startMJD[index], scanXP$endMJD[index])); tempTime <- mjd2doy(mjdSec[index])
	for(file_index in startFileIndex:endFileIndex){
		endPoint   <- min( c(scanXP$endMJD[index] - prefix2MJDsec(prefix[file_index]) + 1, 1800) )
		if( file_index == startFileIndex){
			startPoint <- scanXP$startMJD[index] - prefix2MJDsec(prefix[file_index]) + 1
			C00 <- readPolariS_X(sprintf('%s.C.%02d%s', prefix[file_index], 0, postFix))[,startPoint:endPoint] / BP$BP00
			C01 <- readPolariS_X(sprintf('%s.C.%02d%s', prefix[file_index], 1, postFix))[,startPoint:endPoint] / BP$BP01
		} else {
			startPoint <- 1
			temp <- readPolariS_X(sprintf('%s.C.%02d%s', prefix[file_index], 0, postFix))[,startPoint:endPoint] / BP$BP00; C00 <- cbind(C00, temp)
			temp <- readPolariS_X(sprintf('%s.C.%02d%s', prefix[file_index], 1, postFix))[,startPoint:endPoint] / BP$BP01; C01 <- cbind(C01, temp)
		}
	}
	integRange <- which( Mod(apply(C00, 2, mean)) > median(Mod(apply(C00, 2, mean))) )
	delayC00[index] <- delay_search( apply(C00[chRange,integRange], 1, mean)) * chnum/length(chRange)
	delayC01[index] <- delay_search( apply(C01[chRange,integRange], 1, mean)) * chnum/length(chRange)
	temp00 <- delay_cal(apply(C00[chRange, integRange], 1, mean), delayC00[index] * length(chRange) / chnum )
	temp01 <- delay_cal(apply(C01[chRange, integRange], 1, mean), delayC01[index] * length(chRange) / chnum )
	C00Vis[index] <- mean( temp00 ); C01Vis[index] <- mean( temp01 )
	cat(sprintf('%02d:%02d:%02.0f  %6.4f %6.4f %6.3f    %6.4f %6.4f %6.3f\n', tempTime$hour, tempTime$min, tempTime$sec, delayC00[index], Mod(C00Vis[index]), Arg(C00Vis[index]), delayC01[index], Mod(C01Vis[index]), Arg(C01Vis[index])))
	xspec_amphi_plot( 1:length(temp00), temp00, list(freq='CH', amp='Amplitude', phase='phase [rad]', title=sprintf('%s IF=0 scan=%d %02d:%02d:%02.0f', prefix[1], index, tempTime$hour, tempTime$min, tempTime$sec)))
	xspec_amphi_plot( 1:length(temp00), temp01, list(freq='CH', amp='Amplitude', phase='phase [rad]', title=sprintf('%s IF=1 scan=%d %02d:%02d:%02.0f', prefix[1], index, tempTime$hour, tempTime$min, tempTime$sec)))
}
WG <- data.frame( mjdSec = mjdSec, delay00 = delayC00, delay01 = delayC01, Vis00 = C00Vis, Vis01 = C01Vis)
save(WG, file=sprintf("%s.WG.Rdata", prefix[1]))
#-------- Plot
time_two_plot( (WG$mjdSec%%86400)/3600, WG$delay00, WG$delay01, list(time='UT [hour]', value1='CH0 XY delay [sample]', value2='CH1 XY delay [sample]', title=prefix[1]) )
time_two_plot( (WG$mjdSec%%86400)/3600, Arg(WG$Vis00), Arg(WG$Vis01), list(time='UT [hour]', value1='CH0 XY phase [rad]', value2='CH1 XY phase [rad]', title=prefix[1]) )
dev.off()
