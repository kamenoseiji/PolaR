#-------- Venus
#if(0){
load('/Volumes/SSD/PolariS/20150317/4933282665.uncalStokes.Rdata')	# Venus
Stokes00 <- uncalStokes00; Stokes01 <- uncalStokes01
load('/Volumes/SSD/PolariS/20150316/4933193431.uncalStokes.Rdata')
Stokes00 <- rbind(Stokes00, uncalStokes00); Stokes01 <- rbind(Stokes01, uncalStokes01)
load('/Volumes/SSD/PolariS/20150315/4933108401.uncalStokes.Rdata')
Stokes00 <- rbind(Stokes00, uncalStokes00); Stokes01 <- rbind(Stokes01, uncalStokes01)
load('/Volumes/SSD/PolariS/20150315/4933123243.uncalStokes.Rdata')
Stokes00 <- rbind(Stokes00, uncalStokes00); Stokes01 <- rbind(Stokes01, uncalStokes01)
#-------- Jupiter
load('/Volumes/SSD/PolariS/20150316/4933217002.uncalStokes.Rdata')
Stokes00 <- rbind(Stokes00, uncalStokes00); Stokes01 <- rbind(Stokes01, uncalStokes01)
load('/Volumes/SSD/PolariS/20150317/4933302995.uncalStokes.Rdata')
Stokes00 <- rbind(Stokes00, uncalStokes00); Stokes01 <- rbind(Stokes01, uncalStokes01)
#}
if(0){
#-------- Crab nebula
load('/Volumes/SSD/PolariS/20150317/4933302355.uncalStokes.Rdata')
Stokes00 <- uncalStokes00; Stokes01 <- uncalStokes01
load('/Volumes/SSD/PolariS/20150317/4933286952.uncalStokes.Rdata')
Stokes00 <- rbind(Stokes00, uncalStokes00); Stokes01 <- rbind(Stokes01, uncalStokes01)
load('/Volumes/SSD/PolariS/20150316/4933214077.uncalStokes.Rdata')
Stokes00 <- rbind(Stokes00, uncalStokes00); Stokes01 <- rbind(Stokes01, uncalStokes01)
load('/Volumes/SSD/PolariS/20150316/4933206170.uncalStokes.Rdata')
Stokes00 <- rbind(Stokes00, uncalStokes00); Stokes01 <- rbind(Stokes01, uncalStokes01)
load('/Volumes/SSD/PolariS/20150316/4933197813.uncalStokes.Rdata')
Stokes00 <- rbind(Stokes00, uncalStokes00); Stokes01 <- rbind(Stokes01, uncalStokes01)
load('/Volumes/SSD/PolariS/20150315/4933111618.uncalStokes.Rdata')
Stokes00 <- rbind(Stokes00, uncalStokes00); Stokes01 <- rbind(Stokes01, uncalStokes01)
load('/Volumes/SSD/PolariS/20150315/4933119750.uncalStokes.Rdata')
Stokes00 <- rbind(Stokes00, uncalStokes00); Stokes01 <- rbind(Stokes01, uncalStokes01)
}

plot(Stokes00$Q/Stokes00$I, Stokes00$U/Stokes00$I, xlim=c(-0.25, 0.25), ylim=c(-0.25,0.25), pch=20, cex=0.5)
#arrows( rep(0.0, length(Stokes00$I)), rep(0.0, length(Stokes00$I)), Stokes00$Q/Stokes00$I, Stokes00$U/Stokes00$I, length=0.0)
