args <- commandArgs(trailingOnly = T)
setwd('.')
fileNum <- length(args)

load(args[1])
D_comb <- D
for(file_index in 2:fileNum){
	load(args[2])
	D_comb <- rbind(D_comb, D)
}
D <- D_comb
fileName <- sprintf("%s.Dcomb.Rdata", strsplit(args[1], "\\.")[[1]][1])
save(D, file=fileName)
