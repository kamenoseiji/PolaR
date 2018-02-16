library(FITSio)
fileName <- commandArgs(trailingOnly = T)
setwd('.')
FITS <- readFITS(fileName)
ax1 <- axVec(1, FITS$axDat)
ax2 <- axVec(2, FITS$axDat)
ax3 <- axVec(3, FITS$axDat)
xlab <- FITS$axDat$ctype[1]
ylab <- FITS$axDat$ctype[2]
mapCenter <- c( 15*(04 + 41/60 + 43.87/3600), 25 + 41/60 + 17.7/3600)
velGradRadius <- 74.0/3600  # 74 arcsec (FWHM x 2)
velocRange <- c(5.2, 6.5)
velocIndex <- which( ax3 > velocRange[1] & ax3 < velocRange[2] )
MOM0 <- apply(FITS$imDat[,,velocIndex], c(1,2), sum)
PeakMap <- apply(FITS$imDat[,,velocIndex], c(1,2), max)
thresh <- 0.5; mapIndex <- which(PeakMap > thresh); PeakMap[which(PeakMap < thresh)] <- NA; MOM0[which(MOM0 < thresh)] <- NA
Vupper <- Vlower<- matrix(nrow=nrow(FITS$imDat), ncol=ncol(FITS$imDat))
for(index in mapIndex){
    row_index <- (index - 1) %% nrow(PeakMap) + 1
    col_index <- ceiling( index / nrow(PeakMap) )
    halfLineRange <- which( FITS$imDat[row_index,col_index,velocIndex] > 0.5* PeakMap[row_index,col_index] )
    Vupper[row_index, col_index] <- ax3[velocIndex[max(halfLineRange)]]
    Vlower[row_index, col_index] <- ax3[velocIndex[min(halfLineRange)]]
}
pdf(sprintf('%s.pdf', fileName))
#--------  Peak Intensity Map
image(-ax1+mapCenter[1], ax2-mapCenter[2], PeakMap, xlab = xlab, ylab = ylab, main=fileName, col = topo.colors(32)); points(0,0, pch=3)

#--------  Moment-0 (integrated intensity) Map
image(-ax1+mapCenter[1], ax2-mapCenter[2], MOM0, xlab = xlab, ylab = ylab, main=fileName, col = topo.colors(32)); points(0,0, pch=3)

#-------- Velocity Gradient for Vupper
image(-ax1+mapCenter[1], ax2-mapCenter[2] ,Vupper, xlab = xlab, ylab = ylab, main=fileName, col = cm.colors(32)); points(0,0, pch=3)
DF <- data.frame( x=ax1[(mapIndex - 1) %% nrow(Vupper) + 1] - mapCenter[1], y=ax2[ceiling( mapIndex / nrow(Vupper) )] - mapCenter[2], z=Vupper[mapIndex], w=PeakMap[mapIndex])
DF$w <- DF$w / (DF$x^2 + DF$y^2)
DF[(DF$x^2 + DF$y^2) > velGradRadius^2,]$w <- 0
fit <- lm(formula = z ~ x + y, data=DF, weights = w)
text_sd <- sprintf('Verocity Gradient = (%5.2f %5.2f) km/s/deg', coef(fit)['x']/ cospi(mapCenter[2]/180), coef(fit)['y'])
text(-0.1, 0.15, pos=4, text_sd) 
cat(text_sd); cat('\n')

#-------- Velocity Gradient for Vlower
image(-ax1+mapCenter[1], ax2-mapCenter[2] ,Vlower, xlab = xlab, ylab = ylab, main=fileName, col = cm.colors(32)); points(0,0, pch=3)
DF <- data.frame( x=ax1[(mapIndex - 1) %% nrow(Vlower) + 1] - mapCenter[1], y=ax2[ceiling( mapIndex / nrow(Vlower) )] - mapCenter[2], z=Vlower[mapIndex], w=PeakMap[mapIndex])
DF$w <- DF$w / (DF$x^2 + DF$y^2)
DF[(DF$x^2 + DF$y^2) > velGradRadius^2,]$w <- 0
fit <- lm(formula = z ~ x + y, data=DF, weights = w)
text_sd <- sprintf('Verocity Gradient = (%5.2f %5.2f) km/s/deg', coef(fit)['x']/ cospi(mapCenter[2]/180), coef(fit)['y'])
text(-0.1, 0.15, pos=4, text_sd) 
cat(text_sd); cat('\n')

#-------- Velocity Gradient for mean(Vupper, Vlower)
image(-ax1+mapCenter[1], ax2-mapCenter[2] ,0.5*(Vupper + Vlower), xlab = xlab, ylab = ylab, main=fileName, col = cm.colors(32)); points(0,0, pch=3)
DF <- data.frame( x=ax1[(mapIndex - 1) %% nrow(Vlower) + 1] - mapCenter[1], y=ax2[ceiling( mapIndex / nrow(Vlower) )] - mapCenter[2], z=0.5*(Vupper[mapIndex]+Vlower[mapIndex]), w=PeakMap[mapIndex])
DF$w <- DF$w / (DF$x^2 + DF$y^2)
DF[(DF$x^2 + DF$y^2) > velGradRadius^2,]$w <- 0
fit <- lm(formula = z ~ x + y, data=DF, weights = w)
text_sd <- sprintf('Verocity Gradient = (%5.2f %5.2f) km/s/deg', coef(fit)['x']/ cospi(mapCenter[2]/180), coef(fit)['y'])
text(-0.1, 0.15, pos=4, text_sd) 
cat(text_sd); cat('\n')
dev.off()


#-------- Velocity Gradient for Moment-1 map
if(0){
matrix(nrow=nrow(FITS$imDat), ncol=ncol(FITS$imDat))
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
}
