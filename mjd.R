doy2mjd <- function(year, doy, hh, mm, ss){
	# doy2mjd : Calculate MJD of given UTC
	# Year : Valid for 1901 - 2099
	# doy  : Day of year
	# hh   : Hour in UTC
	# mm   : Minutes
	# ss   : Second
	
	#-------- Check for validity
	if( year < 1901 || year > 2099){ return(-1)}
	
	#-------- Constants
	MJD_1901 <- 15384        # MJD at the beginning of 1901
	DAY_PER_4YEAR <- 1461
	DAY_PER_YEAR  <- 365
	SEC_PER_DAY   <- 86400
	
	#-------- MJD calculation
	mjd <- ((year - 1901) %/% 4)* DAY_PER_4YEAR + ((year - 1901)%%4)*DAY_PER_YEAR + doy + MJD_1901
	sec_day <- ((hh* 60 + mm)*  60 + ss) / SEC_PER_DAY
	return( mjd + sec_day)
}

mjd2doy <- function( mjdSec ){
	MJD_1901 <- 15384        # MJD at the beginning of 1901
	DAY_PER_4YEAR <- 1461
	DAY_PER_YEAR  <- 365
	SEC_PER_DAY   <- 86400
	Sec <- mjdSec %% 60
	Min <- floor(mjdSec / 60) %% 60
	Hour<- floor(mjdSec / 3600) %% 24
	Day <- floor(mjdSec / 86400)
	DaySinceLeapYear <- (Day -  MJD_1901) %% DAY_PER_4YEAR
	Year <- floor((Day - MJD_1901) / DAY_PER_4YEAR)* 4 + 1901 + floor(DaySinceLeapYear/DAY_PER_YEAR)
	Doy  <- DaySinceLeapYear %% DAY_PER_YEAR
	return(list(year=Year, doy=Doy, hour=Hour, min=Min, sec=Sec))
}	

mjd2gmst <- function(mjd, ut1utc = 0){
	# mjd2gmst : Calculate Greenwidge Mean Sidereal Time
	# mjd : Modified Julian Date
	# ut1utc : UT1 - UTC [sec]
	
	FACT <- c(24110.54841, 8640184.812866, 0.093104, 0.0000062)
	MJD_EPOCH <- 51544.5  # MJD at 2000 1/1 12:00:00 UT
	TU_UNIT <- 36525.0
	SEC_PER_DAY   <- 86400
	
	tu <- (mjd - MJD_EPOCH) / TU_UNIT
	ut1 <- (mjd%%1)* SEC_PER_DAY + ut1utc
	gmst <- (ut1 + FACT[1] + ((FACT[4]* tu + FACT[3])* tu + FACT[2])* tu) / SEC_PER_DAY
	return(2* pi* (gmst %% 1))
}

gst2ha <- function(gst, lambda=131.556644*pi/180, ra=pi){
	# gst2ha : Calculate Hour Angle
	# gst : Greenwidge Sidereal Time
	# lambda : Longtude of the site [rad]
	# ra  : Right ascension of the source [rad]
	
	lst <- gst + lambda  # Local Sidereal Time
	ha  <- lst - ra      # Hour Angle
	ha <- (ha + pi)%%(2*pi) - pi
	return(ha)
}

ha2azel <-function(ha, phi=34.217019*pi/180, dec=0){
	# ha2azel : Calculate Azimuth, Elevation, and Parallactic Angle
	# ha  : Hour Angle
	# phi : Latitude of the site [rad]
	# dec : Declination of the source [rad]

	sin_el <- sin(phi)* sin(dec) + cos(phi)* cos(dec)* cos(ha)
	el <- asin(sin_el)
	az <- atan2(cos(dec)*sin(ha), sin(phi)*cos(dec)*cos(ha) - cos(phi)*sin(dec)) + pi
	pa <- atan2(cos(phi)*sin(ha), sin(phi)*cos(dec) - cos(phi)*sin(dec)* cos(ha))
	return(data.frame(ha=ha, az=az, el=el, pa=pa))
}