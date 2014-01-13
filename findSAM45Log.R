library(downloader)
source_url("https://raw.github.com/kamenoseiji/PolaR/master/date.R", , prompt=F, quiet=T)
findSAM45Log <- function(dir, prefix){
	# dir: directory to search SAM45 logging files
	# prefix: PolariS file prefix (=YYYYDOYHHMMSS)
	PolariSmjd <- prefix2MJDsec(prefix)
	fileList <- list.files(dir)[grep("^SAM45", list.files(dir))]	# file list that start by "SAM45"
	SAM45mjd <- numeric(0)
	for( index in 1:length(fileList) ){
		SAM45mjd[index] <- timeString2MJD(unlist(strsplit(fileList[index], "[.]"))[5])	# Date in YYYYMMDDHHMMSS -> MJD
	}
	return(fileList[which.min(abs(SAM45mjd - PolariSmjd))])			# The nearest file
}
