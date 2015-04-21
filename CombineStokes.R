averageSpec <- function( spec, weight ){
    fileNum <- length(weight)
    return( rowSums(matrix( spec, ncol=fileNum ) %*% diag(weight)) / sum(weight) )
}

args <- commandArgs(trailingOnly = T)
setwd('.')
fileNum <- length(args)
SDrange <- 8193:16384
I02 <- numeric(0); I13 <- numeric(0)
Q02 <- numeric(0); Q13 <- numeric(0)
U02 <- numeric(0); U13 <- numeric(0)
V02 <- numeric(0); V13 <- numeric(0)
WT <- numeric(fileNum)
for( file_index in 1:fileNum){
    load(args[file_index])
    I02 <- append(I02, StokesI02); I13 <- append(I13, StokesI13)
    Q02 <- append(Q02, StokesQ02); Q13 <- append(Q13, StokesQ13)
    U02 <- append(U02, StokesU02); U13 <- append(U13, StokesU13)
    V02 <- append(V02, StokesV02); V13 <- append(V13, StokesV13)
    WT[file_index] <- 1.0/var(StokesV02[SDrange])
}
#
StokesI02 <- averageSpec(I02, WT); StokesI13 <- averageSpec(I13, WT)
StokesQ02 <- averageSpec(Q02, WT); StokesQ13 <- averageSpec(Q13, WT)
StokesU02 <- averageSpec(U02, WT); StokesU13 <- averageSpec(U13, WT)
StokesV02 <- averageSpec(V02, WT); StokesV13 <- averageSpec(V13, WT)

fileName <- sprintf("%s.SCOMB.Rdata", strsplit(args[1], "\\.")[[1]][1])
save(StokesI02, StokesI13, StokesQ02, StokesQ13, StokesU02, StokesU13, StokesV02, StokesV13, file=fileName)
cat('Combined Stokes spectra are saved into '); cat(fileName); cat('\n')
