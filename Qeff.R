erf <- function(x) 2* pnorm(x * sqrt(2)) - 1	# Error Function
erf2 <- function(x) 2* pnorm(x) - 1				# erf(x/sqrt(2))

probNbit <- function( a, levelNum=256 ){
	# Input vector a as a[1]=threshold and a[2]=bias scaled by sd.
	volt <- ((1:(levelNum-1)) - levelNum/2)* a[1]	# Threshold voltages
	prob <- 0.5* c((erf2(volt[1] - a[2]) + 1), (erf2(volt[2:(levelNum-1)] - a[2]) - erf2(volt[1:(levelNum-2)] - a[2])), (1 - erf2(volt[levelNum-1] - a[2])))
	return(prob)
}

probThresh <- function( a, Thresh ){
	levelNum <- length(Thresh) + 1
	# Input vector a as a[1]=threshold and a[2]=bias scaled by sd.
	volt <- Thresh* a[1]	# Threshold voltages
	prob <- 0.5* c((erf2(volt[1] - a[2]) + 1), (erf2(volt[2:(levelNum-1)] - a[2]) - erf2(volt[1:(levelNum-2)] - a[2])), (1 - erf2(volt[levelNum-1] - a[2])))
	return(prob)
}

#-------- Initial paramters of Gaussian fit
initGaussNbit <- function( prob, levelNum = 256 ){
	Vweight <- ((1:levelNum) - levelNum/2 - 0.5)		#	Voltage for each level
	Pweight <- Vweight^2								#	Power for each level
	Average <- Vweight %*% prob
	Variance <- Pweight %*%prob - Average^2
	return( c(1/sqrt(Variance), Average/sqrt(Variance)))
}	

#-------- Gaussian parameters determined by bit distribution
gaussNbit <- function( nsample, levelNum ){
	if(length(nsample) < levelNum){	return(NA); warning(message=sprintf("Need %d levels!", levelNum))}
	totalSample <- sum( as.numeric(nsample[1:levelNum]) )		# Total Number of samples
	prob <- nsample[1:levelNum] / totalSample					# Measures probabilities
	weight <- nsample / (1 - prob)^2					# Weight
	W <- diag( weight[1:levelNum] )							# Weight matrix
	SqPi2 <- 1/sqrt(2.0* pi)							# Constant
	a <- initGaussNbit(prob, levelNum)							# Initial Value
	
	niter <- 0
	while(niter < 10){									# Max iterations = 10
		resid <- prob - probNbit( a, levelNum )					# Residual from trial
		
		#-------- Factors used in the partial matrix
		erfDeriv <- c(0.0, exp(-0.5* (((1:(levelNum-1)) - levelNum/2)* a[1] - a[2])^2 ), 0.0)
		WerfDeriv <- ((0:levelNum) - levelNum/2)* erfDeriv
		
		#-------- Partial Matrix
		p <- matrix(nrow=levelNum, ncol=2, SqPi2*c(diff(WerfDeriv), -diff(erfDeriv)))
		t_p_W <- t(p) %*% W
		
		#-------- Solution
		solution <- solve( t_p_W %*% p )
		correction <- solution %*% (t_p_W %*% resid)
		
		#-------- Correction
		a <- a + correction
		
		#-------- Converged?
		if( sum(correction^2) < 1.0e-20 ) break
		niter <- niter + 1
	}
	return( c(a, diag(solution)))	# Return best-fit values and standard error
}

#-------- Gaussian parameters determined by bit distribution
gaussThresh <- function( nsample, Thresh ){
	levelNum <- length(Thresh) + 1
	if(length(nsample) < levelNum){	return(NA); warning(message=sprintf("Need %d levels!", levelNum))}
	totalSample <- sum( as.numeric(nsample[1:levelNum]) )		# Total Number of samples
	prob <- nsample[1:levelNum] / totalSample					# Measures probabilities
	weight <- nsample / (1 - prob)^2					# Weight
	W <- diag( weight[1:levelNum] )							# Weight matrix
	SqPi2 <- 1/sqrt(2.0* pi)							# Constant
	a <- initGaussNbit(prob, levelNum)							# Initial Value
	
	niter <- 0
	while(niter < 10){									# Max iterations = 10
		resid <- prob - probThresh( a, Thresh )			# Residual from trial
		
		#-------- Factors used in the partial matrix
		erfDeriv <- c(0.0, exp(-0.5* (Thresh* a[1] - a[2])^2 ), 0.0)
		WerfDeriv <- c(min(Thresh)-1, Thresh, max(Thresh)+1)* erfDeriv
		
		#-------- Partial Matrix
		p <- matrix(nrow=levelNum, ncol=2, SqPi2*c(diff(WerfDeriv), -diff(erfDeriv)))
		t_p_W <- t(p) %*% W
		
		#-------- Solution
		solution <- solve( t_p_W %*% p )
		correction <- solution %*% (t_p_W %*% resid)
		
		#-------- Correction
		a <- a + correction
		
		#-------- Converged?
		if( sum(correction^2) < 1.0e-20 ) break
		niter <- niter + 1
	}
	return( c(a, diag(solution)))	# Return best-fit values and standard error
}

#-------- Estimation of threshold voltages
threshLevel <- function( nsample ){
	# nlevel <- 256
	# if(length(nsample) < nlevel){	return(-1)}		# Check input sample
    nlevel <- length(nsample)
	gaussParam <- gaussNbit(nsample, nlevel)[1:2]	# threshold and bias
	probReal <- nsample / sum(nsample)				# Fraction at each level
	#-------- Cumulative fraction
	cumReal <- probReal[1:(nlevel-1)]
	for(level_index in 1:(nlevel-2)){
		cumReal <- cumReal + c( rep(0, level_index), probReal[1:(nlevel - level_index - 1)])
	}
    thresh <- rep(Inf, (nlevel - 1))
	tryCatch(
		{ thresh <- (qnorm(cumReal) + gaussParam[2]) / gaussParam[1] },
		warning = function(e){  },
		silent = TRUE
	)
	return(thresh)
}

#-------- biGauss function
biGauss <- function(x, y, rho){
	arg <- (x*x + y*y - 2.0* rho* x* y) / (2.0* (1 - rho* rho))
	fact <- 1.0 / (2* pi* sqrt(1 - rho* rho))
	return( fact* exp(-arg) )
}

#-------- Van Vleck correction
vanvleck4 <- function(thresh, rho){
	hfn <- function(rho){
		result <- 0.0
		for(x_index in 1:15){
			x_thresh <- thresh* (x_index - 8)
			for(y_index in 1:15){
				y_thresh <- thresh* (y_index - 8)
				result <- result + biGauss( x_thresh, y_thresh, rho)
			}
		}
		return(4.0* result)
	}
	
	rho_4 <- integrate( hfn, 0, rho )
	return(rho_4$value)
}

vanvleckN <- function(thresh, rho, levelNum = 256){
	hfn <- function(rho){
		result <- 0.0
		for(x_index in 1:(levelNum - 1)){
			x_thresh <- thresh* (x_index - levelNum/2)
			for(y_index in 1:(levelNum - 1)){
				y_thresh <- thresh* (y_index - levelNum/2)
				result <- result + biGauss( x_thresh, y_thresh, rho)
			}
		}
		return(4.0* result)
	}
	
	rho_N <- integrate( hfn, 0, rho )
	return(rho_N$value)
}
