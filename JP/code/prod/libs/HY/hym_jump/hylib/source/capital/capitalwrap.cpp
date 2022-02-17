// capitalwrap7.cpp : Defines the entry point for the application.
//

//#include "stdafx.h"
#include "stdio.h"
#include "iostream.h"
#include "string.h"
#include "stdlib.h"
#include "capital.h"
#include "cerror.h"
#include "time.h"

#define FNAME "testFile.txt"
#define FILESIZE 50
#define MAXLINE 2000
#define MAXARRAYSIZE 200
#define SEPARATOR ','


int closeFile(FILE *fname) ;
FILE* openReadOnlyFile ( FILE *fname ) ;
int openFile ( FILE *fname ) ;
int extractInputArray( char *ch1, char ch2[200][2000]) ;
int loadDoubleArray(double numberArray[], char *tempLine)  ;
int loadIntegerArray(long numberArray[], char *tempLine)  ;


FILE *theFile ;



int main()
{
//
	void GtoErrMsgOn();

	char  tempOutput[FILESIZE][MAXLINE] ;
	int	ctr ;
	char tempLine[MAXLINE];
//	char field[200][MAXLINE] ;


	ctr = 0 ;
	tempLine[0] = '\0';

	while (ctr < FILESIZE )
	{
			*tempOutput[ctr++] = '\0' ;
	}


   if( (theFile = fopen("testFile.txt","r")) == NULL )
      printf( "The file 'data' was not opened\n" );
   else
      printf( "The file 'data' was opened\n" );

   /* load up the file */


	ctr = 0;
	while (!feof (theFile))
	{
		fgets(tempLine,MAXLINE,theFile);
		strcpy(tempOutput[ctr++],tempLine);
	}

	closeFile (theFile);

	/*
	* Now that the file has been read, fill the arrays
	*/


	


	double spotPrice[MAXARRAYSIZE];
	loadDoubleArray(spotPrice, tempOutput[0]) ;

	long divDates[MAXARRAYSIZE];
	loadIntegerArray(divDates, tempOutput[1]) ;

	double divRates[MAXARRAYSIZE];
	loadDoubleArray(divRates, tempOutput[2]) ;

	double dps[MAXARRAYSIZE];
	loadDoubleArray(dps, tempOutput[3]) ;

	long repoDates[MAXARRAYSIZE];
	loadIntegerArray(repoDates, tempOutput[4]) ;

	double repoRates [MAXARRAYSIZE];
	loadDoubleArray(repoRates, tempOutput[5]) ;

	long swapDates[MAXARRAYSIZE];
	loadIntegerArray(swapDates, tempOutput[6]) ;

	double swapRates[MAXARRAYSIZE];
	loadDoubleArray(swapRates, tempOutput[7]) ;

	long volDates[MAXARRAYSIZE];
	loadIntegerArray(volDates, tempOutput[8]) ;

	double volRates[MAXARRAYSIZE];
	loadDoubleArray(volRates, tempOutput[9]) ;

	double volShift[MAXARRAYSIZE];
	loadDoubleArray(volShift, tempOutput[10]) ;

	long instType[MAXARRAYSIZE];
	loadIntegerArray(instType, tempOutput[11]) ;
 
	double notional[MAXARRAYSIZE];
	loadDoubleArray(notional, tempOutput[12]) ;

	double recoveryRate[MAXARRAYSIZE];
	loadDoubleArray(recoveryRate, tempOutput[13]) ;

	long cfAccStartDates[MAXARRAYSIZE];
	loadIntegerArray(cfAccStartDates, tempOutput[14]) ;

	long cfAccEndDates[MAXARRAYSIZE];
	loadIntegerArray(cfAccEndDates, tempOutput[15]) ;

	long cfDates[MAXARRAYSIZE];
	loadIntegerArray(cfDates, tempOutput[16]) ;

	double cfCoupons[MAXARRAYSIZE];
	loadDoubleArray(cfCoupons, tempOutput[17]) ;

	double cfAmort[MAXARRAYSIZE];
	loadDoubleArray(cfAmort, tempOutput[18]) ;

	long cfClaimDates[MAXARRAYSIZE];
	loadIntegerArray(cfClaimDates, tempOutput[19]) ;

	double cfClaimAmounts[MAXARRAYSIZE];
	loadDoubleArray(cfClaimAmounts, tempOutput[20]) ;

	long exerStartDates[MAXARRAYSIZE];
	loadIntegerArray(exerStartDates, tempOutput[21]) ;

	long exerEndDates[MAXARRAYSIZE];
	loadIntegerArray(exerEndDates, tempOutput[22]) ;

	double exerStartStrikes[MAXARRAYSIZE];
	loadDoubleArray(exerStartStrikes, tempOutput[23]) ;

	double exerEndStrikes[MAXARRAYSIZE];
	loadDoubleArray(exerEndStrikes, tempOutput[24]) ;

	long optionDir[MAXARRAYSIZE];
	loadIntegerArray(optionDir, tempOutput[25]) ;

	long optionType[MAXARRAYSIZE];
	loadIntegerArray(optionType, tempOutput[26]) ;

	double optionBarrier[MAXARRAYSIZE];
	loadDoubleArray(optionBarrier, tempOutput[27]) ;

	long exerciseType[MAXARRAYSIZE];
	loadIntegerArray(exerciseType, tempOutput[28]) ;

	double lim1[MAXARRAYSIZE];
	loadDoubleArray(lim1, tempOutput[29]) ;

	double lim2[MAXARRAYSIZE];
	loadDoubleArray(lim2, tempOutput[30]) ;

	double vollim[MAXARRAYSIZE];
	loadDoubleArray(vollim, tempOutput[31]) ;

	double x[MAXARRAYSIZE];
	loadDoubleArray(x, tempOutput[32]) ;

	double lim[MAXARRAYSIZE];
	loadDoubleArray(lim, tempOutput[33]) ;

	double beta[MAXARRAYSIZE];
	loadDoubleArray(beta, tempOutput[34]) ;

	long assetProcessType[MAXARRAYSIZE];
	loadIntegerArray(assetProcessType, tempOutput[35]) ;

	long ppy[MAXARRAYSIZE];
	loadIntegerArray(ppy, tempOutput[36]) ;

	long valueDates[MAXARRAYSIZE];
	loadIntegerArray(valueDates, tempOutput[37]) ;


	/************
	double spotPrice[MAXARRAYSIZE];
	loadDoubleArray(spotPrice, tempOutput[0]) ;

	double spotPrice[MAXARRAYSIZE];
	loadDoubleArray(spotPrice, tempOutput[0]) ;

	double spotPrice[MAXARRAYSIZE];
	loadDoubleArray(spotPrice, tempOutput[0]) ;

	double spotPrice[MAXARRAYSIZE];
	loadDoubleArray(spotPrice, tempOutput[0]) ;

	double spotPrice[MAXARRAYSIZE];
	loadDoubleArray(spotPrice, tempOutput[0]) ;
	************/

	char outChar[10] = "WQWQWQ\0";
	double outputs[40];
  double elapsedTime = 0.00;
	long* dummyTime = 0;
  long startTime;
  long endTime;
  char* hyVersion ;


     hyVersion = HY_version() ;

	startTime = time(dummyTime);


	int callReturn =
	HYMCapitalWrapper(spotPrice,
					  divDates,
					  divRates,
					  dps,
					  repoDates,
					  repoRates,
					  swapDates,
					  swapRates,
					  volDates,
					  volRates,
					  volShift,
					  instType,
						notional,
					  recoveryRate,
					  cfAccStartDates,
					  cfAccEndDates,
					  cfDates,
					  cfCoupons,
					  cfAmort,
					  cfClaimDates,
					  cfClaimAmounts,
					  exerStartDates,
					  exerEndDates,
					  exerStartStrikes,
					  exerEndStrikes,
					  optionDir,
					  optionType,
					  optionBarrier,
					  exerciseType,
					  lim1,
					  lim2,
					  vollim,
					  x,
					  lim,
					  beta,
					  assetProcessType,
					  ppy,
					  valueDates,
					  outChar,
					  outputs);

	/*
	printf("price = %f\n, delta = %f\n,gamma = %f\n,vega = %f\n,accr=%f\n,annuity = %f\n",
		outputs[1],outputs[3],outputs[4],outputs[5],outputs[7],outputs[13]);
    */

	endTime = time(dummyTime);
	elapsedTime = difftime(endTime, startTime) ;

	printf("pricerVersion = %s\n price = %f\n, delta = %f\n,gamma = %f\n,first = %f\n,callReturn=%d\n, elapsedTime=%f\n, text = %s\n",
		hyVersion, outputs[1],outputs[3],outputs[4],outputs[0],callReturn, elapsedTime, &outChar);
    
        int status = 1;


  






	return 0;
}

int loadIntegerArray(long numberArray[], char *tempLine)  
{

	int arraySize;
	char thisField[100];
	char field[200][MAXLINE] ;
//	char tempLine[MAXLINE]

	/* spot Price  */
//	double numberArray[MAXARRAYSIZE];
	int ctr = 0;
	while(ctr<MAXARRAYSIZE) numberArray[ctr++] = 0;
	arraySize = 0;
//	strcpy(tempLine, tempOutput[0]) ;
	extractInputArray(tempLine, field);
	arraySize =  atoi(field[1]) + 1;
	
	if (arraySize > 0)
	{
		int ctr = 0;
		while (ctr < arraySize)
		{
			strcpy(thisField,field[++ctr]);
			numberArray[ctr-1]=atoi(thisField);			
		}
		
	}

		return 0;
}


int loadDoubleArray(double numberArray[], char *tempLine)  
{

	int arraySize;
	char thisField[100];
	char field[200][MAXLINE] ;
//	char tempLine[MAXLINE]

	/* spot Price  */
//	double numberArray[MAXARRAYSIZE];
	int ctr = 0;
	while(ctr<MAXARRAYSIZE) numberArray[ctr++] = 0.0;
	arraySize = 0;
//	strcpy(tempLine, tempOutput[0]) ;
	extractInputArray(tempLine, field);
	arraySize =  atoi(field[1]) + 1;
	
	if (arraySize > 0)
	{
		int ctr = 0;
		while (ctr < arraySize)
		{
			strcpy(thisField,field[++ctr]);
			numberArray[ctr-1]=atof(thisField);			
		}
		
	}

		return 0;
}



int extractInputArray (char *inputLine, char outputArray[200][2000])
{

	int ctr = 0;
	int charNo = 0;
	int fieldNo = 0;
	char tempString[MAXLINE] ;
	char ch;

	ctr = 0;
	tempString[0] = '\0' ;

	while (ctr < 200 )
	{
			*outputArray[ctr++] = '\0';
	}

	ctr = 0;
	while ((ch = inputLine[ctr++]) != '\0')
	{
		if ((ch != SEPARATOR) && (ch != 10))
		{				
				tempString[charNo++]=ch ;				
		}
		else
		{
				tempString[charNo]='\0' ;
				strcpy(outputArray[fieldNo++],tempString) ;
				*tempString = '\0' ;
				charNo = 0;
		}

	}
	return 0;
}




int openFile(char fname)
{
	return 0;
}


FILE* openReadOnlyFile ( char *fname )
{

	return fopen(fname, "r");

}



int closeFile(FILE *fname)
{

		fclose(fname) ;
		return 0;

}


int closeAllFiles ()
{

//		_fcloseAll () ;
	
	return 0;
}


int writeCrapToFile (char fname, char inputText)
{
	return 0;
}
