//**************************************************************
//
//
//
//**************************************************************
				// DR C++ Library
#include "kstdinc.h"
#include "kutilios.h"
#include "kmodpar.h"

				// Products
#include "krate.h"
				// VTree Tools
#include "vtlbase.h"
#include "vtltest.h"
				// ALIb Tree
#include "kvtalib.h"



extern	"C" {
#include "date_sup.h"

#include "drlstr.h"		// Strtok
#include "drltime.h"
#include "drlio.h"
#include "drlmem.h"
#include "drlproc.h"
#include "drlts.h"
#include "drlstr.h"
#include "drlgetop.h"		// getopt() 

#include "dritkwrp.h"		// TDrWrapper
};

#define SUP_MAX_NUM_ZC 3


#define	__PROGRAM__	"alfrate_t"

static	void	usage(char *s, ...);

#define	__DEBUG__
#undef	__DEBUG__




//--------------------------------------------------------------
//

int
main(int argc, char **argv)
{
static	char			routine[] = __PROGRAM__;
	char			drwDir[1024];
	char			dealFnam[1024];
	char			xdlFnam[1024];	// instrument filename
	int			c;

	KVTreeAL		vt;

	TDrWrapperData		*drWrapData = NULL;

	int		numVolDates;
	TDate		*volDates = NULL;
	double		*volMat = NULL;
	int		*volFreq = NULL;
	double		*volRates = NULL;


	KVector(TDate)	resetDates;
	KVector(TDate)	payDates;
	KVector(double)	strikes;
	KRate		floatRate;
	string		indxCurveName;
	KZCurve		indxZcCurve;
	string		discCurveName;
	KZCurve		discZcCurve;






    try {

	KTreeModelParam		modParam[3];
	KDealModelParam		dealParam;



	// Enable logging
	//
	DppLogging(-1, "run.log", 1, NULL);


	//
	// Defaults
	//
	strcpy(drwDir, ".");
	//strcpy(drwDir, "tmp.drWrapper");
	strcpy(dealFnam, __PROGRAM__".dat");
	strcpy(xdlFnam, "");

	//
	// Parse options
	//
  	while ((c = getopt(argc, argv, "D:W:h")) != EOF)
  	switch (c) {
	case 'D':
		ASSERT_OR_THROW(sscanf(optarg, "%s", xdlFnam) == 1);
		dppLog << "Using xdl file: " << xdlFnam << endl;
		break;
	case 'W':
		ASSERT_OR_THROW(sscanf(optarg, "%s", drwDir) == 1);
		dppLog << "Using directory: " << drwDir << endl;
		break;
	default:
	case 'h':
		usage(NULL);
		break;
	}
	argv += optind - 1 ;
	/*if (argv[1] == NULL) usage(NULL);*/



	//
	// (1) Read data 
	//
	if (*(++argv) == NULL) {

	} else {
		strcpy(dealFnam, *argv);
		dppLog << "Using dat file: " << dealFnam << endl;
	}


	// Read Dr Wrapper data
	IF_FAILED_THROW( DriTDrWrapperDataGetFull(
		drwDir,
		DRI_DRW_TYPE2_2CURVES,
		&drWrapData));


	// Read model parameters
	KTreeModelParamKapitalDrWrapperRead(
		drwDir,
		modParam);



	//throw KFailure ("Done parsing (abort).\n");
	dppLog << "===========================================================================" << endl;

	//
	// Open deal data file
	//
	ifstream is(dealFnam);
	check(is, dealFnam);



	//
	// Read deal file
	//
#if defined(V10)
	if (strlen(xdlFnam) <= 0) {
		strcpy(xdlFnam, getString(is, "#DealInputFile"));
	} else {
		getString(is, "#DealInputFile");
	}
	IF_FAILED_THROW( DrlStrSubsEnv(xdlFnam));
#else
	strcpy(xdlFnam, dealFnam);
#endif


	//
	// Read test data
	//
	{
		int	numItems, idx;

		floatRate.Get(is, TRUE);
		numItems = getInt(is, "numDates");
		resetDates.resize(numItems);
		payDates.resize(numItems);
		for (idx=0; idx<numItems; idx++) {
			resetDates[idx] =  getTDate(is, "reset date");
			payDates[idx]   =  getTDate(is, "pay date");
		}
		numItems = getInt(is, "numStrikes");
		strikes.resize(numItems);
		for (idx=0; idx<numItems; idx++) {
			strikes[idx] = getDouble(is, "strikes");
		}

	}



	//
	// Get model parameters overwrite
	//
	dealParam.GetDrwOverwrite(is, 3, modParam);
	dppLog << "Deal Pricing Parameters:" << endl;
	dppLog << dealParam;
	dppLog << "===========================================================================" << endl;



	//
	// Get volatilities
	//

	IF_FAILED_THROW( DriTDrWrapperDataGetInterpVol(
		drWrapData,
		dealParam.mCalibFinal,
		(dealParam.mFinalDate ? dealParam.mMatDate : 0L),
		dealParam.mMatInt,
		&numVolDates,
		&volDates,
		&volMat,
		&volFreq,
		&volRates));

	TDate		todayDate = drWrapData->fToday;
	int		numZcCurves = 2;		// SUP_MAX_NUM_ZC;  
	TCurve		*zcCurves[SUP_MAX_NUM_ZC];
	KVector(string)	zcNames;

	zcCurves[0] = drWrapData->fZcCurve;
	zcNames.push_back("Curve1");        

	zcCurves[1] = drWrapData->fDiscZcCurve;
	zcNames.push_back("Curve2");

	zcCurves[2] = drWrapData->fRiskZcCurve;
	//zcNames.push_back("Curve3"); 

	string diffuseName("Curve1");


	vt.Initialize(
		//todayDate,
			drWrapData->fZcCurve->fBaseDate,
		numZcCurves,
		zcCurves,
		zcNames,
		diffuseName,
		dealParam.mNumFact,
		dealParam.mBeta,
		dealParam.mAlpha,
		dealParam.mRho,
		dealParam.mPpy,
		dealParam.mSmoothFact,
		dealParam.mNumStdevCut,
		dealParam.mQ1,
		dealParam.mQ2,
		0.,
        	numVolDates,
        	volDates,
        	volMat,
        	volFreq,
        	volRates);


	dppLog << "===========================================================================" << endl;



	//
	// Perform Test
	//
	indxCurveName = "Curve1";
	indxZcCurve   = drWrapData->fZcCurve;

	discCurveName = "Curve2";
	discZcCurve = drWrapData->fDiscZcCurve;

	ostrstream	ostOut;
	//ofstream	osOut(__PROGRAM__".out");


	KVTreeTestForwardRates1(
	        vt,
	        resetDates,
	        payDates,
	        strikes,
	        floatRate,
	        indxCurveName,
	        indxZcCurve,
	        discCurveName,
	        discZcCurve,
	        ostOut);





	cout.write(ostOut.str(), ostOut.pcount());




	// Free memory
	DriTDrWrapperDataFree(drWrapData);


    }
    catch (KFailure) {
	// Free memory
	DriTDrWrapperDataFree(drWrapData);

	dppErr << format("%s: failed.\n", routine);
	return(FAILURE);
    }
    catch (...) {
	dppErr << format("%s: failed (uncaught).\n", routine);
	return(FAILURE);
    }
    return(SUCCESS);
}







/*------------------------------------------------------
 *
 */

static	void
usage(char *s, ...)
{
	va_list	ap;
	char	buf[255];

	if (s != NULL) {
		va_start(ap, s);
		vsprintf(buf, s, ap);
		va_end(ap);

		fprintf(stderr, "ERROR: %s (-h for help)\n", buf);
		exit(1);
	} else {
		fprintf(stdout, "\
USAGE\n\
    "__PROGRAM__" [options] [<input dat file>]\n\
OPTIONS\n\
    -W <dir>    Read DrWrapper data files from <dir>\n\
    -D <file>   Use <file> for XDL input\n\
INPUT FILE FORMAT\n\
<template file>\n\
<old environment wrapper dir>\n\
<new environment wrapper dir>\n\
<portfolio>\n\
<factor file>\n\
\n", __PROGRAM__);
		fprintf(stdout, "    Created %s  %s\n",
			__DATE__,__TIME__);
		fflush(stdout);
		exit(0);
	}
}




