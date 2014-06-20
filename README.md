This repository contains R scripts that are used in analysis of polarization spectra taken by PolariS.
To include R scripts in the GitHub repository, follow the instruction:

(0) If RCurl package is not installed, get it from Package Installer
(1) library(RCurl)
(2) eval(parse(text = getURL("https://raw.githubusercontent.com/kamenoseiji/PolaR/master/readSAM45.R", ssl.verifypeer = FALSE)))
 # Replace tha file name with the source file you want.

If you want to commit to the repository, send pull request.
