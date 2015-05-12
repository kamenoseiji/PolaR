# DelayInterP
# usage: Rscript ContStokes.R [Scan.Rdata file name] [SPEC.Rdata file name] [WG.Rdata file name] [BP file name]
#
#library(RCurl)
setwd('.')
options(digits = 4)
#-------- InterLeave
interLeave <- function( vector ){
    dataNum <- length(vector)
    results <- numeric(2* dataNum - 1)
    for(index in 1:(dataNum - 1)){
        results[2*index - 1] <- vector[index]
        results[2*index] <- 0.5* (vector[index] + vector[index + 1])
    }
    results[2* dataNum - 1] <- vector[dataNum]
    return(results)
}
#-------- Load Spec and Scan data
args <- commandArgs(trailingOnly = T)
setwd('.')
load(args[1])	 #Load delay file
#
#-------- Smoothed Delay and Phase
dataNum <- length(WG$mjdSec)
outDF <- data.frame(
    mjdSec = sort( c(WG$mjdSec, WG$mjdSec[1:(dataNum-1)] + diff(WG$mjdSec)/2)),
    delay00 = interLeave(WG$delay00),
    delay01 = interLeave(WG$delay01),
    Vis00 = complex( real=interLeave(Re(WG$Vis00)), imaginary=interLeave(Im(WG$Vis00))),
    Vis01 = complex( real=interLeave(Re(WG$Vis01)), imaginary=interLeave(Im(WG$Vis01))) )

WG <- outDF
fileName <- sprintf("%s.Delay.Rdata", strsplit(args[1], "\\.")[[1]][1])
save(WG, file=fileName)
cat(sprintf('InterLeaved Delay is saved into %s\n', fileName))
