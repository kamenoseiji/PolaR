SDrange <- 1025:8192
load('/Volumes/SSD/PolariS/20150418/2015108024946.STOKES.Rdata')
tmpI02 <- StokesI02
tmpQ02 <- StokesQ02
tmpU02 <- StokesU02
tmpV02 <- StokesV02
tmpI13 <- StokesI13
tmpQ13 <- StokesQ13
tmpU13 <- StokesU13
tmpV13 <- StokesV13
SD02 <- sd(tmpV02[SDrange])
SDtmp13 <- sd(tmpV13[SDrange])
load('/Volumes/SSD/PolariS/20150419/2015109020627.STOKES.Rdata')
SD02 <- c(sd(tmpV02[SDrange]), sd(StokesV02[SDrange]))
SD13 <- c(sd(tmpV13[SDrange]), sd(StokesV13[SDrange]))
WT02 <- 1.0 / SD02^2
WT13 <- 1.0 / SD13^2

StokesI02 <- (WT02[1]* tmpI02 + WT02[2]* StokesI02) / sum(WT02)
StokesQ02 <- (WT02[1]* tmpQ02 + WT02[2]* StokesQ02) / sum(WT02)
StokesU02 <- (WT02[1]* tmpU02 + WT02[2]* StokesU02) / sum(WT02)
StokesV02 <- (WT02[1]* tmpV02 + WT02[2]* StokesV02) / sum(WT02)

StokesI13 <- (WT13[1]* tmpI13 + WT13[2]* StokesI13) / sum(WT13)
StokesQ13 <- (WT13[1]* tmpQ13 + WT13[2]* StokesQ13) / sum(WT13)
StokesU13 <- (WT13[1]* tmpU13 + WT13[2]* StokesU13) / sum(WT13)
StokesV13 <- (WT13[1]* tmpV13 + WT13[2]* StokesV13) / sum(WT13)

save(StokesI02, StokesI13, StokesQ02, StokesQ13, StokesU02, StokesU13, StokesV02, StokesV13, file='COMB.Stokes.Rdata')