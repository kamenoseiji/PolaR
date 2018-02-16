# PlotBP
# usage: Rscript PlotBP [BP file name]
#
setwd('.')
#-------- Load BP file
BPfile <- commandArgs(trailingOnly = T)
setwd('.')
load(BPfile)
IFnum = length(BP)
BPprefix <- strsplit(BPfile, 'Rdata')[[1]]
#-------- Plot in PDF
pdf(paste(BPprefix, 'pdf', sep=''))
par(oma = c(0, 0, 0, 0), mfcol = c(IFnum,2))
for(index in 1:IFnum){
    amp <- abs(BP[[index]])
    phs <- Arg(BP[[index]])* 180.0/pi
    plot(amp, ylim=c(0, max(amp)), type='l', ylab='Amp', xlab='', main=paste(BPprefix, sprintf('IF%d', index)))
    plot(phs, ylim=c(-180, 180), pch=20, cex=0.2, ylab='Phase [deg]', xlab='Channel')
}
dev.off()
