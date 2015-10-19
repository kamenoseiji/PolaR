library(FITSio)
fileName <- commandArgs(trailingOnly = T)
setwd('.')
FITS <- readFITS(fileName)
ax1 <- axVec(1, FITS$axDat)
ax2 <- axVec(2, FITS$axDat)
ax3 <- axVec(3, FITS$axDat)
xlab <- FITS$axDat$ctype[1]
ylab <- FITS$axDat$ctype[2]
mapCenter <- c( 15*(04 + 41/60 + 42.47/3600), 25 + 41/60 + 27.1/3600)
velocRange <- c(5.2, 6.5)
velocIndex <- which( ax3 > velocRange[1] & ax3 < velocRange[2] )
MOM0 <- apply(FITS$imDat[,,velocIndex], c(1,2), sum)
thresh <- 5.0
mapIndex <- which( MOM0 > thresh ); MOM0[MOM0 < thresh] <- NA #flagIndex <- which( MOM0 < thresh ); MOM0[flagIndex] <- NA
MOM1 <- matrix(nrow=nrow(FITS$imDat), ncol=ncol(FITS$imDat))
for(index in mapIndex){
	row_index <- (index - 1) %% nrow(MOM0) + 1
	col_index <- ceiling( index / nrow(MOM0) )
	MOM1[index] <- FITS$imDat[ row_index, col_index, velocIndex] %*% ax3[velocIndex] / MOM0[index]
}
pdf(sprintf('%s.pdf', fileName))
image(-ax1+mapCenter[1], ax2-mapCenter[2], MOM0, xlab = xlab, ylab = ylab, main=fileName, col = topo.colors(32)); points(0,0, pch=3)
image(-ax1+mapCenter[1], ax2-mapCenter[2] ,MOM1, xlab = xlab, ylab = ylab, main=fileName, col = cm.colors(32)); points(0,0, pch=3)
DF <- data.frame( x=ax1[(mapIndex - 1) %% nrow(MOM0) + 1] - mapCenter[1], y=ax2[ceiling( mapIndex / nrow(MOM0) )] - mapCenter[2], z=MOM1[mapIndex], w=MOM0[mapIndex])
fit <- lm(formula = z ~ x + y, data=DF, weights = w)
text_sd <- sprintf('Verocity Gradient = (%5.2f %5.2f) km/s/deg', coef(fit)['x']/ cospi(mapCenter[2]/180), coef(fit)['y'])
text(-0.1, 0.15, pos=4, text_sd) 
cat(text_sd); cat('\n')
dev.off()
