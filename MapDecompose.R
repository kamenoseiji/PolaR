library(FITSio)
fileName <- commandArgs(trailingOnly = T)
setwd('.')
FITS <- readFITS(fileName)
ax1 <- axVec(1, FITS$axDat)
ax2 <- axVec(2, FITS$axDat)
ax3 <- axVec(3, FITS$axDat) * 1.0e-3
xlab <- FITS$axDat$ctype[1]
ylab <- FITS$axDat$ctype[2]
mapCenter <- c( 277.4884, -1.945744)
velGradRadius <- 74.0/3600  # 74 arcsec (FWHM x 2)
velocRange <- c(6.2, 8.9)
velocIndex <- which( ax3 > velocRange[1] & ax3 < velocRange[2] )
veloc <- ax3[velocIndex]
MOM0 <- apply(FITS$imDat[,,velocIndex], c(1,2), sum)
PeakMap <- apply(FITS$imDat[,,velocIndex], c(1,2), max)
thresh <- 0.5; mapIndex <- which(PeakMap > thresh); PeakMap[which(PeakMap < thresh)] <- NA; MOM0[which(MOM0 < thresh)] <- NA
PeakIndex <- which.max(MOM0)
row_Peak <- (PeakIndex - 1) %% nrow(MOM0) + 1
col_Peak <- ceiling( PeakIndex / nrow(MOM0) )
pdf(sprintf('%s.pdf', fileName))
V1 <- V2 <- E1 <- E2 <- Pl <- Pm <- numeric(9)
for( row_index in 1:3 ){
    for( col_index in 1:3 ){
        index <-  3 * (row_index - 1)  + col_index
        DF <- data.frame(veloc=veloc, Ta=FITS$imDat[row_Peak + row_index -2, col_Peak + col_index -2, velocIndex])
        plot(DF, pch=20, xlab='velocity [km/s]', ylab='Ta* [K]', main=sprintf('%s [%.1f, %.1f]', fileName, 3600*(ax1[row_Peak + row_index-2] - ax1[row_Peak]), 3600*(ax2[col_Peak + col_index-2] - ax2[col_Peak])))
        fit <- nls(formula=Ta ~ a* exp(-0.5*((veloc - b)/c)^2) + d* exp(-0.5*((veloc - e)/f)^2), DF, start=list(a=0.72*max(DF$Ta), b=7.18, c=0.3, d=0.6*max(DF$Ta), e=7.8, f=0.33), control=nls.control(maxiter=1000))
        lines(DF$veloc, predict(fit))
        text_sd <- sprintf('V1 = %.3f (%.3f) km/s  V2 = %.3f (%.3f) km/s', coef(summary(fit))[2,1], coef(summary(fit))[2,2], coef(summary(fit))[5,1], coef(summary(fit))[5,2])
        text(6.5, max(DF$Ta), text_sd, pos=4)
        V1[index] <- coef(fit)['b']; V2[index] <- coef(fit)['e']
        E1[index] <- coef(summary(fit))[2,2]; E2[index] <- coef(summary(fit))[5,2]
        Pl[index] <- ax1[row_Peak + row_index-2] - ax1[row_Peak]; Pm[index] <- ax2[col_Peak + col_index-2] - ax2[col_Peak]
    }
}
dev.off()
#-------- Velocity Gradient for V1
DF <- data.frame( x=Pl, y=Pm, z=V1, w=1.0/E1^2)
fit <- lm(formula = z ~ x + y, data=DF, weights = w)
text_sd <- sprintf('Component 1: Velocity Gradient = %5.2f %5.2f (%.1f %.1f) km/s/deg', coef(fit)['x'], coef(fit)['y'], summary(fit)$coefficient[2,2], summary(fit)$coefficient[3,2])
cat(text_sd); cat('\n')

#-------- Velocity Gradient for V2
DF <- data.frame( x=Pl, y=Pm, z=V2, w=1.0/E2^2)
fit <- lm(formula = z ~ x + y, data=DF, weights = w)
text_sd <- sprintf('Component 2: Velocity Gradient = %5.2f %5.2f (%.1f %.1f) km/s/deg', coef(fit)['x'], coef(fit)['y'], summary(fit)$coefficient[2,2], summary(fit)$coefficient[3,2])
cat(text_sd); cat('\n')
