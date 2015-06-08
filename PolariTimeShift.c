// PolariBunch: A tool to bunch spectrum
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "shm_k5data.inc"
#define	MAX(a,b)	a>b?a:b	// Larger Value
#define	FNAME	1
#define	SHIFT	2
#define ARGNUM	3

int soy2dhms(
    unsigned long   sc,         // IN: Second of Year
    int   *doy_ptr,   // Out: Pointer to Day of Year
    int   *hh_ptr,    // Out: Pointer to Hour
    int   *mm_ptr,    // Out: Pointer to Minute
    int   *ss_ptr)    // Out: Pointer to Second
{
    *doy_ptr = sc / 86400 + 1;
    *hh_ptr  = (sc % 86400 ) / 3600;
    *mm_ptr  = (sc % 3600 ) / 60;
    *ss_ptr  =  sc % 60;
    if( sc != (((*doy_ptr - 1)* 24 + (*hh_ptr))* 60 + (*mm_ptr))* 60 + (*ss_ptr)){
        return(-1);
    }
    return(0);
}

unsigned long dhms2soy(
    unsigned long   doy,    // IN: Day of Year
    unsigned long   hh,     // IN: Hour
    unsigned long   mm,     // IN: Min
    unsigned long   ss)     // IN: Sec
{
    return((((doy - 1)* 24 + hh)* 60 + mm)* 60 + ss);
}

int	timeShift(
	char	*fname,			// IN: Input File Name
	int		shift)	    	// IN: Time shift in sec
{
	FILE	*file_ptr;				// File Pointer
	char	outName[24];			// Output File Name
	char	prefix[14];			    // Output File Name
	struct	SHM_PARAM	param;		// File Header
    unsigned long    startSoY, shiftSoY;     // Second of Year

	//-------- Open File
	if((file_ptr = fopen(fname, "rb+")) == NULL){ return(-1);}	// Open input file

	//-------- Read and Write Header
	fread(&param, sizeof(struct SHM_PARAM), 1, file_ptr);

	//-------- Time Shift in Header
	// startSod = hms2sod(param.hour, param.min, param.sec);
	// shiftSod = startSod + shift;
	// sod2hms(shiftSod, &(param.hour), &(param.min), &(param.sec));
	startSoY = dhms2soy(param.doy, param.hour, param.min, param.sec);
	shiftSoY = startSoY + shift;
	soy2dhms(shiftSoY, &(param.doy), &(param.hour), &(param.min), &(param.sec));

    strcpy(outName, fname);
    sprintf(prefix, "%04d%03d%02d%02d%02d", param.year, param.doy, param.hour, param.min, param.sec );
    strncpy(outName, prefix, 13);

	//-------- Skip to the start position and overwrite the header
	rewind(file_ptr);
	fwrite(&param, sizeof(struct SHM_PARAM), 1, file_ptr);

	//-------- Save to File
	fclose(file_ptr);
    rename(fname, outName);

	return(0);
}

int main(
	int		argc,		// Number of Arguments
	char	**argv)		// Pointer to Arguments
{
	if( argc < ARGNUM ){
		printf("USAGE: PolariBunch [file name] [BUNCH] [Start PP] [End PP] !!\n");
		printf("  file name : input file name (e.g. 2014016154932.C.00)\n");
		printf("  Shift     : Time shift in [sec]. Positive value makes time stamp greater.\n");
		exit(-1);
	}
	timeShift(argv[FNAME], atoi(argv[SHIFT]));
	return(0);
}
