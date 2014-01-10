erf <- function(x) 2* pnorm(x * sqrt(2)) - 1	# Error Function
erf2 <- function(x) 2* pnorm(x) - 1				# erf(x/sqrt(2))

#-------- Probabilities in 16 levels
prob4bit <- function( a ){
	# Input vector a as a[1]=threshold and a[2]=bias scaled by sd.
	volt <- (-7:7)* a[1]	# Threshold voltages
	prob <- 0.5* c((erf2(volt[1] - a[2]) + 1), (erf2(volt[2:15] - a[2]) - erf2(volt[1:14] - a[2])), (1 - erf2(volt[15] - a[2])))
	return(prob)
}

#-------- Initial paramters of Gaussian fit
initGauss4bit <- function( prob ){
	Vweight <- (-7:8 - 0.5)			#	Voltage for each level
	Pweight <- (-7:8 - 0.5)^2		#	Power for each level
	Average <- Vweight %*% prob
	Variance <- Pweight %*%prob - Average^2
	return( c(1/sqrt(Variance), Average/sqrt(Variance)))
}	

#-------- Gaussian parameters determined by bit distribution
gauss4bit <- function( nsample ){
	if(length(nsample) < 16){	return(NA); warning(message="Need 16 levels!")}
	totalSample <- sum( as.numeric(nsample[1:16]) )		# Total Number of samples
	prob <- nsample[1:16] / totalSample					# Measures probabilities
	weight <- nsample / (1 - prob)^2					# Weight
	W <- diag( weight[1:16] )							# Weight matrix
	SqPi2 <- 1/sqrt(2.0* pi)							# Constant
	a <- initGauss4bit(prob)							# Initial Value
	
	niter <- 0
	while(niter < 10){									# Max iterations = 10
		resid <- prob - prob4bit( a )					# Residual from trial
		
		#-------- Factors used in the partial matrix
		erfDeriv <- c(0.0, exp(-0.5* ((-7:7)* a[1] - a[2])^2 ), 0.0)
		WerfDeriv <- (-8:8)* erfDeriv
		
		#-------- Partial Matrix
		p <- matrix(nrow=16, ncol=2, SqPi2*c(diff(WerfDeriv), -diff(erfDeriv)))
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
