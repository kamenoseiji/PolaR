lineLabel <- c('F=5-5','F=4-3','F=5-4','F=6-5','F=4-4')
restFreq  <- c(45488.8368, 45490.2614, 45490.3137, 45490.3373, 45492.1085)
lineStrength <- c(0.928, 18.053, 22.282, 27.432, 0.928)
Ta <- c(0.6, 3.0, 3.3, 3.6, 0.6)
HC3N <- data.frame( label=lineLabel, freq=restFreq, strength=lineStrength, Ta=Ta)
#vplot <- (args$restFreq - HC3N$freq) / args$restFreq * 299792.458 + args$trackVel
# aplot <- 2.5* HC3N$strength / max(HC3N$strength)
#aplot <- 0.1* HC3N$strength / max(HC3N$strength)
#arrows( vplot, rep(0.0, length(vplot)), vplot, aplot, angle=90, length=0, col='blue')
#text(vplot, aplot, lineLabel, srt=45, adj=0, cex=0.5)
plot(HC3N$strength, HC3N$Ta, pch=20, xlab='Line Strength', ylab='Ta* [K]', main='HC3N in TMC-1', xlim=c(0,40))
text(lineStrength, Ta, lineLabel, adj=-0.2, cex=0.5)
fit <- nls(formula = Ta ~ a * (1.0 - exp(-b* lineStrength)), start=list(a=5.8, b=0.1))
x <- 0:90
lines(x, predict(fit, data.frame(lineStrength = x)))
