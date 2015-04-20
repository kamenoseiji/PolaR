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

int	hms2sod(
	int		hour,		// IN: Hour
	int		min,		// IN: Minute
	int		sec)		// IN: Sec
{
	return( sec + 60*(min + 60*hour));
}

int	sod2hms(
	int		sod,		// IN: Second of the day
	int		*hour,		// Hour
	int		*min,		// Minute
	int		*sec)		// Sec
{
	*hour = sod / 3600;
	*min  = (sod % 3600) / 60;
	*sec  = sod % 60;
	return(sod - hms2sod(*hour, *min, *sec));
}

int	timeShift(
	char	*fname,			// IN: Input File Name
	int		shift)	    	// IN: Time shift in sec
{
	FILE	*file_ptr;				// File Pointer
	char	outName[24];			// Output File Name
	char	prefix[14];			    // Output File Name
	struct	SHM_PARAM	param;		// File Header
    long    startSod, shiftSod;     // Second of Year

	//-------- Open File
	if((file_ptr = fopen(fname, "rb+")) == NULL){ return(-1);}	// Open input file

	//-------- Read and Write Header
	fread(&param, sizeof(struct SHM_PARAM), 1, file_ptr);

	//-------- Time Shift in Header
	startSod = hms2sod(param.hour, param.min, param.sec);
	shiftSod = startSod + shift;
	sod2hms(shiftSod, &(param.hour), &(param.min), &(param.sec));

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
