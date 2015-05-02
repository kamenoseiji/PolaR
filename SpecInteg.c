// PolariBunch: A tool to bunch spectrum
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "shm_k5data.inc"
#define	MAX(a,b)	a>b?a:b	// Larger Value
#define	FNAME	1
#define STARTPP	2
#define ENDPP	3
#define ARGNUM	4

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

int	fileType( char	*fname)
{
	return(
		strstr(fname, ".A.")?1:
		strstr(fname, ".C.")?2:0);
}

int timeIntegReal(
    int     chNum,          // IN: Number of spectral channel
    int     timeNum,        // IN: Number of time records to be integrated
    float   *inData,        // IN: Pointer to input original data
    float   *outData)       // OUT:Pointer to output integrated data)
{
    int     timeIndex, chIndex;

    for(timeIndex=0; timeIndex<timeNum; timeIndex++){
        for(chIndex=0; chIndex<chNum; chIndex++){
            outData[chIndex] += inData[timeIndex* chNum + chIndex];
        }
    }
    return(timeIndex);
}

int	specTimeInteg(
	char	*fname,			// IN: Input File Name
	int		filetype,		// IN: 1->autocorr, 2->crosscorr
	int		StartPP,		// IN: Rec# to start extract
	int		EndPP,          // IN: Final Rec# to extract
    float   *outdata)		// OUT:Pointer to integrated spectral data
{
	FILE	*file_ptr;				// File Pointer
	struct	SHM_PARAM	param;		// File Header
	float	*readdata;				// Power spectrum data record
	int		readSize, outSize;		// Original Record size [bytes]
	int		startSod, endSod;		// Start and End second of day
	int		startH, startM, startS, endH, endM, endS;
	int		integNum;

	//-------- Open File
	if((file_ptr = fopen(fname, "r")) == NULL){ return(-1);}	// Open input file

	//-------- Read file header
	fread(&param, sizeof(struct SHM_PARAM), 1, file_ptr);
    integNum = EndPP - StartPP + 1; if(integNum < 1){   return(-1);}
    outSize  = filetype* param.num_ch* sizeof(float);
    readSize  = integNum* outSize;
	if( readSize == 0 ){	perror("Unavailable Data Format"); return(-1);}

	readdata = (float *)malloc(readSize);
	startSod = hms2sod(param.hour, param.min, param.sec);
	endSod = startSod + EndPP; startSod += StartPP;
	sod2hms(startSod, &startH, &startM, &startS);
	sod2hms(endSod,   &endH,   &endM,   &endS);
	printf("SpecInteg.c: %02d:%02d:%02d - %02d:%02d:%02d\n", startH, startM, startS, endH, endM, endS);

	//-------- Skip to the start position
	fseek(file_ptr, StartPP* outSize, SEEK_CUR);

	//-------- Read and Write Records
	if(fread(readdata, readSize, 1, file_ptr) != 1){
        printf("File Read Error [%s]\n", fname);
        return(-1);
    }
	fclose(file_ptr);

    //-------- Time-integration
    timeIntegReal(filetype* param.num_ch, integNum, readdata, outdata); 
	
	//-------- Close File
	free(readdata);
	return(integNum);
}

int fileInteg(
    char    *fname,
    char    *fsave,
    int     StartPP,
    int     EndPP)
{
    int     outSize;
    float   *outData;
    FILE    *file_ptr;
    FILE    *save_ptr;
    struct  SHM_PARAM   param;
    int     filetype;

    //-------- Read File Pointer
    if(( filetype = fileType(fname)) == 0){  return(-1);}
    if((file_ptr = fopen(fname, "r")) == NULL){ return(-1);}    // Open input file
    fread(&param, sizeof(struct SHM_PARAM), 1, file_ptr);
    fclose(file_ptr);

    outSize = filetype * param.num_ch * sizeof(float);
    if( outSize == 0){ printf("Invalid File Type!\n"); return(-1); }
    outData = (float *)malloc( outSize ); memset(outData, 0, outSize);
    specTimeInteg( fname, filetype, StartPP, EndPP, outData );

    save_ptr = fopen(fsave, "w");
    fwrite( outData, outSize, 1, save_ptr);
    fclose( save_ptr);
    free(outData);
    return(outSize); 
}

int main(
	int		argc,		// Number of Arguments
	char	**argv)		// Pointer to Arguments
{
	if( argc < ARGNUM ){
		printf("USAGE: SpecInteg [file name] [Start PP] [End PP] !!\n");
		printf("  file name : input file name (e.g. 2014016154932.C.00)\n");
		printf("  Start PP  : Start record number to extract\n");
		printf("  End   PP  : Final record number to extract\n");
		exit(-1);
	}
	fileInteg(argv[FNAME], "tmp.spec", atoi(argv[STARTPP]), atoi(argv[ENDPP]));
	return(0);
}
