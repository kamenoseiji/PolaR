#-------- Parse command-line arguments
parseArg <- function( args ){
    argNum <- length(args)
    chShift <- rep(0, 1024)
    IFselect <- rep(1, 1024)
    Mirror  <- rep(1, 1024)     # Mirror parity (X-Y cable inversion)
    fileName <- c()
    fileIndex <- 0
    for( index in 1:argNum ){
        switch( substr(args[index], 1,2),
            "-C" = chShift[fileIndex] <- round(as.numeric(substring(args[index], 3))),
            "-I" = IFselect[fileIndex] <- as.numeric(substring(args[index], 3)),
            "-M" = Mirror[fileIndex] <- -1,
            { fileIndex = fileIndex + 1; fileName[fileIndex] <- args[index]}
        )
    }
    fileNum <- length(fileName)
    return( list(fileName = fileName, chShift = chShift[1:fileNum], IFselect = IFselect[1:fileNum], Mirror = Mirror[1:fileNum]) )
}

ShiftCH <- function( spec, shift){
    if(shift == 0){ return(spec)}
    chNum <- length(spec)
    if(shift > 0){ return( c( spec[(chNum - shift + 1):chNum], spec[1:(chNum - shift)])) }
    return( c(spec[(-shift + 1):chNum], spec[1:(-shift)]))
}

averageSpec <- function( spec, weight ){
    fileNum <- length(weight)
    return( rowSums(matrix( spec, ncol=fileNum ) %*% diag(weight)) / sum(weight) )
}

argList <- parseArg(commandArgs(trailingOnly = T))
#cat(argList$fileName); cat('\n')
#cat(argList$chShift); cat('\n')
#cat(argList$IFselect); cat('\n')
#cat(argList$Mirror); cat('\n')
fileNum <- length(argList$fileName)
setwd('.')
SDrange <- 8193:16384
StokesI <- numeric(0)
StokesQ <- numeric(0)
StokesU <- numeric(0)
StokesV <- numeric(0)
WT <- numeric(fileNum)
for( file_index in 1:fileNum){
    load(argList$fileName[file_index])
    if( argList$IFselect[file_index] == 2 ){
        StokesI02 <- StokesI13
        StokesQ02 <- StokesQ13
        StokesU02 <- StokesU13
        StokesV02 <- StokesV13
    }
    StokesI <- append(StokesI, ShiftCH(StokesI02, argList$chShift[file_index]))
    StokesQ <- append(StokesQ, argList$Mirror[file_index]* ShiftCH(StokesQ02, argList$chShift[file_index]))
    StokesU <- append(StokesU, ShiftCH(StokesU02, argList$chShift[file_index]))
    StokesV <- append(StokesV, argList$Mirror[file_index]* ShiftCH(StokesV02, argList$chShift[file_index]))
    WT[file_index] <- 1.0/var(StokesV02[SDrange])
    cat(sprintf('%s : SD=%6.2e K\n', argList$fileName[file_index], sd(StokesV02[SDrange])))
}
#
StokesI02 <- averageSpec(StokesI, WT)
StokesQ02 <- averageSpec(StokesQ, WT)
StokesU02 <- averageSpec(StokesU, WT)
StokesV02 <- averageSpec(StokesV, WT)

fileName <- sprintf("%s.SCOMB.%d.Rdata", strsplit(argList$fileName[1], "\\.")[[1]][1], argList$IFselect[1])
save(StokesI02, StokesQ02, StokesU02, StokesV02, file=fileName)
cat('Combined Stokes spectra are saved into '); cat(fileName); cat('\n')
