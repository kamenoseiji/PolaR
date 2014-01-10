#-------- Constants for time/date system
DEGRAD <- pi/180.0
MJD_1901 <- 15384		# MJD at Y1901
DAY_PER_YEAR <- 365
DAY_PER_4YEAR <- 1461
SEC_PER_DAY <- 86400
DOY_MON <- c(0, 31, 59, 90, 120, 151, 181, 212, 242, 273, 303, 334)				# DOY at the 1st day of month
DOY_MON_LEAP <- c(0, 31, 60, 91, 121, 152, 182, 213, 243, 274, 304, 335)		# DOY at the 1st day of monthe (for leap year)
#-------- Calculate Day of year from Month and date
md2doy <- function(year, month, date){
	is_leap <- ((year%%4 == 0) && (month > 2))	# Leap Year Frag
	DOY_MON[month] + date + is_leap
}

#-------- Calculate Modified Julian Date from Year and Day of Year
doy2mjd <- function( year, doy ){
	if((year < 1901) || (year > 2099)){ return(-1)}
	year <- year - 1901
	mjd <- year %/%4 * DAY_PER_4YEAR + year%%4 * DAY_PER_YEAR + doy + MJD_1901
	return(mjd)
}

#-------- Calculate Modified Julian Date from Year, Month, and Day
day2mjd <- function(year, month, day){
	doy <- md2doy(year, month, day)
	doy2mjd(year, doy)
}

#-------- Calculate (fractional) Modified Julian Date
doy2fmjd <- function(year, doy, hour, min, sec){
	if((year < 1901) || (year > 2099)){ return(-1)}
	year <- year - 1901
	mjd <- year %/%4 * DAY_PER_4YEAR + year%%4 * DAY_PER_YEAR + doy + MJD_1901
	sod <- (hour*60 + min)*60 + sec
	mjd <- mjd + sod/SEC_PER_DAY
	return(mjd)
}

#-------- Calculate (fractional) Modified Julian Date in unit of second. This unit is used in CASA
doy2mjdSec <- function(year, doy, hour, min, sec){
	if((year < 1901) || (year > 2099)){ return(-1)}
	year <- year - 1901
	mjd <- year %/%4 * DAY_PER_4YEAR + year%%4 * DAY_PER_YEAR + doy + MJD_1901
	sod <- (hour*60 + min)*60 + sec
	return(mjd* SEC_PER_DAY + sod)
}

#-------- Convert Time String (JST) into MJD in second
timeString2MJD <- function( timeString ){
	year <- as.integer(substring(timeString, 1, 4))
	month <- as.integer(substring(timeString, 5, 6))
	day <- as.integer(substring(timeString, 7, 8))
	jst_hour <- as.integer(substring(timeString, 9, 10))
	minute <- as.integer(substring(timeString, 11, 12))
	second <- as.integer(substring(timeString, 13, 14))
	return(doy2mjdSec(year, md2doy(year, month, day), jst_hour-9, minute, second))
}
 
#-------- Convert Time String (YYYDOYHHMMSS) into MJD in second
prefix2MJDsec <- function(prefix){
	year <- as.integer(substring(prefix, 1, 4))
	doy <- as.integer(substring(prefix, 5, 7))
	hour <- as.integer(substring(prefix, 8, 9))
	minute <- as.integer(substring(prefix, 10, 11))
	second <- as.integer(substring(prefix, 12, 13))
	return(doy2mjdSec(year, doy, hour, minute, second))
}

#-------- Calculate Greenwich Mean Sidereal Time
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

#-------- Calculate Hour Angle
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

#-------- Calculate Azimuth, Elevation, and Parallactic angle
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
