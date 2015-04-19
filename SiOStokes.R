load('/Volumes/SSD/PolariS/20150417SiO/2015107022058.STOKES.Rdata')
freq <- (0:65535)*4e6/65536
source('~/Programs/PolaR/readPolariS.R')
Q <- bunch_vec(StokesQ02, 16)
U <- bunch_vec(StokesU02, 16)
V <- bunch_vec(StokesV02, 16)
I <- bunch_vec(StokesI02, 16)
freqCenter <- 4.312209e10
veloc <- (34000.0 - (bunch_vec(freq, 16) - 2.0e6) / freqCenter * 299792458) * 1.0e-3
lineRange <- which( (veloc < 40) & (veloc > 27))
plotRange <- which( (veloc < 45) & (veloc > 10))
fit <- lm( formula=z~0+x, data=data.frame(x=U[lineRange],z=V[lineRange]) )
# rotation <- atan2( fit$coefficients[2], fit$coefficients[1] )
rotation <- atan2( mean(V[lineRange]), mean(U[lineRange]))
cs <- cos(rotation)
sn <- sin(rotation)
Uc <- cs* U - sn* V
Vc <- sn* U + cs* V
cols <- c('black', 'blue', 'green', 'red')
Stokes <- c('I', 'Q', 'U', 'V')
plot(veloc[plotRange], I[plotRange], type='l', ylim=c(-40, 250), xlab='LSR Velocity [km/s]', ylab='Ta [K]', main='NML Tau SiO Maser (v=1, J=1-0)', col=cols[1])
lines( veloc[plotRange], Q[plotRange], col=cols[2])
lines( veloc[plotRange], Uc[plotRange], col=cols[3])
lines( veloc[plotRange], Vc[plotRange], col=cols[4])
legend('topright', legend=Stokes, col=cols, lty=rep(1, 4))