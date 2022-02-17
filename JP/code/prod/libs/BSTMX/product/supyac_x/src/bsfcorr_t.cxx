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
#include "vtbstest.h"
				// Basis Tree
#include "kbirtree.h"

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

extern 	int      debugLevel;


#define	__PROGRAM__	"bsfcorr_t"

static	void	usage(char *s, ...);

#define	__DEBUG__
#undef	__DEBUG__



//--------------------------------------------------------------
// Compares and prints forward rates computed by the tree
// and the implied tree correlation.
//

int
main(int argc, char **argv)
{
static	char			routine[] = __PROGRAM__;

	char			drwDir[1024];
	char			dealFnam[1024];
	int                     c;

	char			loggingLevel[1024];	// logging level

	KVector(TDate)  resetDates;
        KVector(TDate)  payDates;
        KVector(double) strikes;

        KRate           floatRate0;
        string          indxCurveName0;
        KZCurve         indxZcCurve0;

        KRate           floatRate1;
        string          indxCurveName1;
        KZCurve         indxZcCurve1;

	TDate           today;          // today's date
	KMarketCurves   mktCurves;      // curves and curve types
	KVolDiag        irVolDiag;      // IR volatility data.
	KMrParam        irMrParam;      // IR mr data.
	KSmileParam     irSmileParam;   // IR skew data.
	KVolDiag        bsVolDiag;      // Basis volatility data.
	KMrParam        bsMrParam;      // Basis mr data.
	KSmileParam     bsSmileParam;   // Basis skew data.
	double          irBsCorr;       // IR basis correlation.
 
        string          discCurveName;
	double		discZeroShift;
        KZCurve         discZcCurve;
 
	KMrParam        treeMrParam;    // full tree mr parameters
 
	KBirTree	vt;
	KResetBank	vpResetBank;
 
	bool		isBasis;	

        ostrstream      ostOut;

    try {



	// Enable logging
	//
	DppLogMsgSetFile("run.log");


	//
	// Defaults
	//
	strcpy(drwDir, ".");
	strcpy(dealFnam, __PROGRAM__".dat");
	strcpy(loggingLevel, "0");

	//
        // Parse options
        //
        while ((c = getopt(argc, argv, "D:W:l:h")) != EOF)
        switch (c) {
        case 'W':
                ASSERT_OR_THROW(sscanf(optarg, "%s", drwDir) == 1);
                dppLog << "Using directory: " << drwDir << endl;
                break;
	case 'd':
		ASSERT_OR_THROW(sscanf(optarg, "%s", loggingLevel) == 1);
		dppLog << "Logging level: " << loggingLevel << endl;
		break;
        default:
        case 'h':
                usage(NULL);
                break;
        }
        argv += optind - 1 ;
 
 
        //
        // (1) Read data
        //
        if (*(++argv) == NULL) {
 
        } else {
                strcpy(dealFnam, *argv);
                dppLog << "Using dat file: " << dealFnam << endl;
        }

	// Set debugging level
	debugLevel = atoi(loggingLevel);

	dppLog << "============================================="
		  "==============================" << endl;

	//
	// Open deal data file
	//
	
	ifstream is(dealFnam);
	check(is, dealFnam);


	//
        // Read test data
        //
        {
		int     numItems, idx;
 
		floatRate0.Get(is, TRUE);
		floatRate1.Get(is, TRUE);

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
        // Get the market curves
	// Return true if basis curves are required
        //
        isBasis = BasisTWrapperRead_Mod(
			is,
			drwDir,
			'B',
 
			mktCurves,
			irVolDiag,
			irMrParam,
			irSmileParam,
			bsVolDiag,
			bsMrParam,
			bsSmileParam,
			irBsCorr);


	// We need the full tree parameters
	// to compute the vnfm composite vols,
	treeMrParam = Correlate(irMrParam, bsMrParam, irBsCorr);
 
 
	// Initialize the tree
	//

	vt.Initialize( 
		mktCurves,
		irVolDiag,
		irMrParam,
		irSmileParam,
		bsVolDiag,
		bsMrParam,
		bsSmileParam,
		irBsCorr,
		vpResetBank);



	//
        // Perform Test
        //
	dppLog << "============================================="
		  "==============================" << endl;

	//
        // Perform Test
        //
	int cv0 = mktCurves.mCVTypes[0];
	int cv1 = mktCurves.mCVTypes[1];
	int cv2 = mktCurves.mCVTypes[2];
	
        discCurveName = mktCurves.mCVNames[cv0];
	KMap(int, KZCurve)::iterator itCV0=mktCurves.mCV.find(cv0);
	if(itCV0 != mktCurves.mCV.end())
        	discZcCurve = (*itCV0).second;
	else
		throw KFailure("%s: invalid curve type %d. "
			"Not found from the environment.\n", 
			routine, cv0);
 
        indxCurveName0 = mktCurves.mCVNames[cv1];
	KMap(int, KZCurve)::iterator itCV1=mktCurves.mCV.find(cv1);
	if(itCV1 != mktCurves.mCV.end())
        	indxZcCurve0 = (*itCV1).second;
	else
		throw KFailure("%s: invalid curve type %d. "
			"Not found from the environment.\n", 
			routine, cv1);
 
        indxCurveName1 = mktCurves.mCVNames[cv2];
	KMap(int, KZCurve)::iterator itCV2=mktCurves.mCV.find(cv2);
	if(itCV2 != mktCurves.mCV.end())
        	indxZcCurve1 = (*itCV2).second;
	else
		throw KFailure("%s: invalid curve type %d. "
			"Not found from the environment.\n", 
			routine, cv2);
 
 
	// Zero shift from value date to today for discount curve.
	//
	discZeroShift = mktCurves.ZeroShift(discCurveName);



	// Zero shift from value date to today for discount curve.
	//
	discZeroShift = mktCurves.ZeroShift(discCurveName);

	//
	// Perform the test
	//

	KVTreeTestBasisRatesVolAndCorr(
                &vt,

		irVolDiag,
		bsVolDiag,

		treeMrParam,		// full model parameters

                resetDates,
                payDates,
                strikes,

                floatRate0,
                indxCurveName0,
                indxZcCurve0,

                floatRate1,
                indxCurveName1,
                indxZcCurve1,

                discCurveName,
                discZcCurve,
                ostOut);


	cout.write(ostOut.str(), ostOut.pcount());

	delete [] ostOut.str();


    }
    catch (KFailure) {

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
    -d <dbg>    Specify debug levle\n\
\n", __PROGRAM__);
		fprintf(stdout, "    Created %s  %s\n",
			__DATE__,__TIME__);
		fflush(stdout);
		exit(0);
	}
}

