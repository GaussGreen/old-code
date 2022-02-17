//**************************************************************
//
//
//
//**************************************************************
//#include <cstdio>
				// DR C++ Library
#include "kstdinc.h"
#include "kutilios.h"
#include "kmodpar.h"

				// VTree Tools
#include "vpipars.h"
#include "vtlbase.h"
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


#define	__PROGRAM__	"cet_t"

static	void	usage(char *s, ...);

#define	__DEBUG__
#undef	__DEBUG__


//--------------------------------------------------------------
// Reads and loads the model/market data.
//

bool
BasisTWrapperRead(
	istream& is,			// (I) deal input stream
	char *pathDir,			// (I) directory to read from
	char drwType,			// (I) '2', 'B'asis

	TDate &today,			// (O) today's date 
	KMarketCurves&	mktCurves,	// (O) curves and curve types
	KVolDiag&	irVolDiag,	// (O) IR volatility data.
	KMrParam&	irMrParam,	// (O) IR mr data.
	KSmileParam&	irSmileParam,	// (O) IR skew data.
	KVolDiag&	bsVolDiag,	// (O) Basis volatility data.
	KMrParam&	bsMrParam,	// (O) Basis mr data.
	KSmileParam&	bsSmileParam,	// (O) Basis skew data.
	double&		irBsCorr)	// (O) IR basis correlation.

{
static	char	routine[] = "BasisTWrapperRead";
		
	TDrWrapperData	*drWrapData = NULL;
	long		drw_type = DRI_DRW_TYPE2_2CURVES;

			// --- zero curves

			// --- vol indices
	char		idx1s[128], idx2s[128];
	KVolCalibIndex	volIdx1, volIdx2;


	KMrParam	treeMrParam;

	int		bsDim;
	int		bsType;
	char		bsLiborCVName[256];
	char		bsDiscCVName[256];
	double		bsDelayShift,
			bsQ1, bsQ2, bsQF;




			// --- volatilities for calibration
	int	numVolDates;
	TDate	*volDates = NULL;
	double	*volMat = NULL;
	int	*volFreq = NULL;
	double	*volRates = NULL;

			// --- basis volatility curve
	bool	isBasis = false;	// default
	TCurve  *bsVolCurve = NULL;
        double  bsVolInterp;
        KVector(double) bsVolTmp;

    try {

	//----------------------------------------------
	// (1) Read Dr Wrapper data
	//----------------------------------------------

        IF_FAILED_THROW( DriTDrWrapperDataGetFull(
                pathDir,
                DRI_DRW_BASIS,
                &drWrapData));


	//----------------------------------------------
	// Today date
	//----------------------------------------------
	today = drWrapData->fToday;


	//----------------------------------------------
	// Load zero curves
	//----------------------------------------------

	isBasis = mktCurves.ReadDrw(is, drWrapData);
	

	//----------------------------------------------
	// Calibration indices
	//----------------------------------------------
	strcpy(idx1s, getString(is, "calibration index 1"));
	strcpy(idx2s, getString(is, "calibration index 2"));

	DppReadFromString(idx1s, volIdx1);
	DppReadFromString(idx1s, volIdx2);

	if (!(volIdx1 == volIdx2)) {
		throw KFailure("%s: different volatility calibration"
			" indices (%s and %s).\n",
			routine,
			DppWriteToString(volIdx1).c_str(),
			DppWriteToString(volIdx2).c_str());
	}

	//----------------------------------------------
	// Distribution type
	//----------------------------------------------
	irSmileParam.ReadDrw(is, drWrapData);



	//----------------------------------------------
	// Read full tree mr parameters
	//----------------------------------------------
	treeMrParam.ReadDrw(is, drWrapData);

	//----------------------------------------------
	// Read extra basis parameters and fill structures.
	//----------------------------------------------

	// Read wrapper data 
	if (isBasis) {

		bsDim  = getInt(is, "basis dimension");
		bsType = getInt(is, "basis type (0=Spd, 1=Percent)");
		strcpy(bsLiborCVName, getString(is, "basis libor curve name"));
		strcpy(bsDiscCVName, getString(is, "basis disc curve name"));
		bsDelayShift = getDouble(is, "basis delay shift");
		bsQ1 = getDouble(is, "basis smile Q1");
		bsQ2 = getDouble(is, "basis smile Q2");
		bsQF = getDouble(is, "basis smile QF");

		//
		// Fill market curve basis info
		//
		switch (bsType) {
		case 0:
			mktCurves.mBSType = SPREAD;
			break;
		case 1:
			mktCurves.mBSType = PERCENT;
			break;
		default:
			throw KFailure("%s: bad basis type %d (only 0 or 1).\n",
				routine, bsType);
		}
		mktCurves.mLiborCVName = string(bsLiborCVName);
		mktCurves.mBSDiscCVName = string(bsDiscCVName);

		if (bsDelayShift >= 0e0 && bsDelayShift <= 1e0)
			mktCurves.mBSDelayShift = bsDelayShift;
		else
			throw KFailure("%s: invalid basis delay shift (%f). "
				       "Must be within the range [0, 1].\n",
					routine, bsDelayShift);

		//
		// Fill basis mr parameters
		//
		bsMrParam.mNumFact = bsDim;
		switch (bsDim) {
		case 0:
			// No basis
			irMrParam = treeMrParam;

			bsMrParam    = 0;
			bsMrParam.mNumFact = bsDim;

			bsSmileParam = 0;

			break;
		case 1:
			switch (treeMrParam.mNumFact) {
			case 2:
				// Keep only first dim in IR
				irMrParam = treeMrParam;
				irMrParam.mNumFact = 1;

				// Get second dim for basis
				bsMrParam = treeMrParam;
				bsMrParam.mNumFact = 1;
				bsMrParam.mBeta[0]  = irMrParam.mBeta[1];
				bsMrParam.mAlpha[0] = irMrParam.mAlpha[1];
				bsMrParam.mBackboneQ = 0;
	
				// IR basis correlation.
				irBsCorr = treeMrParam.mRho[0];

				break;
			default:
				throw KFailure("%s: tree dim %d and bs dim %d "
					"(only support 1+1 mode).\n", 
					routine, treeMrParam.mNumFact,
					bsDim);
			}
			break;
		default:
			throw KFailure("%s: only supports basis dimension of 0 or 1.\n",
				routine);
		}

		//
		// Basis smile parameters
		//
		bsSmileParam.mQ1 = 1e0 - bsQ1;
		bsSmileParam.mQ2 = 1e0 - bsQ2;
		bsSmileParam.mQF = bsQF;
		bsSmileParam.mNumIter = 0;

	} else {
		//
		// No basis present (dummy params)
		//
		irMrParam = treeMrParam;

		bsMrParam    = 0;
		bsMrParam.mNumFact = 0;
		bsSmileParam = 0;
		irBsCorr = 0;
	}



	//----------------------------------------------
	// IR Volatilities: interpolate at calib index
	//----------------------------------------------

	irVolDiag = KVolDiag(volIdx1, drWrapData);


	//----------------------------------------------
	// Basis volatilities: use same timeline as IR
	//----------------------------------------------
	bsVolDiag = irVolDiag;
	bsVolDiag = 0e0;

	// Basis volatility curve
	//
	if (isBasis) {
	    bsVolCurve = drWrapData->fBSVolCurve;
	    if (bsVolCurve == NULL) {
		// If no basis curve (must be a type II)
		// set vols to 1 to use the weight.
		if (drwType == 'B') {
			throw KFailure("%s: basis volatility curve "
			       "(basisvol.dat) is NOT available "
			       "from the environment.\n",
				routine);
		}
		bsVolDiag = 1e0;

	    } else {
		// Basis vol curve available, interpolate
        	for (int i=0; i<bsVolDiag.mVolRates.size(); i++)
        	{
			// interpolate basis vol on these benchmark dates
			IF_FAILED_THROW(GtoInterpRate(
					irVolDiag.mVolDates[i],
					bsVolCurve,
					GTO_LINEAR_INTERP,
					&bsVolInterp));

			bsVolDiag.mVolRates[i] = bsVolInterp;
			bsVolDiag.mVolMats[i]  = 0e0;
			bsVolDiag.mVolFreqs[i] = 0;
			bsVolDiag.mSpotVolRates[i] = 0e0;
		}
	    }
	} else {
		// No basis
		bsVolDiag = 0e0;
	}



	//----------------------------------------------
	// Logging
	//----------------------------------------------

	dppLog << "===========================================================================" << endl;
	dppLog << routine << ": INPUT MARKET AND MODEL DATA " << endl;
	dppLog << "===========================================================================" << endl;

	dppLog << "---------------------------------------------------------------------------" << endl;
	dppLog << "DATA: today:\n" << DrlTDatePrint(NULL, today) << endl;
	dppLog << "---------------------------------------------------------------------------" << endl;
	dppLog << "DATA: mktCurves:\n" << mktCurves << endl;
	dppLog << "---------------------------------------------------------------------------" << endl;
	dppLog << "DATA: irVolDiag:\n" << irVolDiag << endl;
	dppLog << "---------------------------------------------------------------------------" << endl;
	dppLog << "DATA: irMrParam:\n" << irMrParam << endl;
	dppLog << "---------------------------------------------------------------------------" << endl;
	dppLog << "DATA: irSmileParam:\n" << irSmileParam << endl;
	dppLog << "---------------------------------------------------------------------------" << endl;
	dppLog << "DATA: bsVolDiag:\n" << bsVolDiag << endl;
	dppLog << "---------------------------------------------------------------------------" << endl;
	dppLog << "DATA: bsMrParam:\n" << bsMrParam << endl;
	dppLog << "---------------------------------------------------------------------------" << endl;
	dppLog << "DATA: bsSmileParam:\n" << bsSmileParam << endl;
	dppLog << "---------------------------------------------------------------------------" << endl;
	dppLog << "DATA: irBsCorr:\n" << irBsCorr << endl;
	dppLog << "===========================================================================" << endl;

	// Free memory
	DriTDrWrapperDataFree(drWrapData);

	return(isBasis);

    }
    catch (KFailure)
    {
	// Free memory
	DriTDrWrapperDataFree(drWrapData);

	throw KFailure("%s: failed (dir=`%s').\n",
		routine, (pathDir ? pathDir : "."));
		
    }
}





//--------------------------------------------------------------
//

int
main(int argc, char **argv)
{
static	char			routine[] = __PROGRAM__;

	char			drwDir[1024];
	char			dataFnam[1024];
	char			insFnam[1024];	// instrument filename
	int                     c;

	char			loggingLevel[1024];	// logging level


	TDate		today;		// today's date 
	KMarketCurves	mktCurves;	// curves and curve types
	KVolDiag	irVolDiag;	// IR volatility data.
	KMrParam	irMrParam;	// IR mr data.
	KSmileParam	irSmileParam;	// IR skew data.
	KVolDiag	bsVolDiag;	// Basis volatility data.
	KMrParam	bsMrParam;	// Basis mr data.
	KSmileParam	bsSmileParam;	// Basis skew data.
	double		irBsCorr;	// IR basis correlation.



	KBirTree		vt;


	bool			isBasis;	
	double			discZeroShift;
	string			discCurveName;	// product discount curve name 

	KVPAtom			*vpRoot = NULL;
	KVPToolAtom		*vpToolRoot = NULL;
	KResetBank		vpResetBank;

    try {



	// Enable logging
	//
	DppLogging(-1, "run.log", 1, NULL);


	//
	// Defaults
	//
	strcpy(drwDir, ".");
	strcpy(dataFnam, __PROGRAM__".dat");
	strcpy(insFnam, "");

	strcpy(loggingLevel, "0");

	//
        // Parse options
        //
        while ((c = getopt(argc, argv, "D:I:W:l:h")) != EOF)
        switch (c) {
        case 'D':
                ASSERT_OR_THROW(sscanf(optarg, "%s", dataFnam) == 1);
                dppLog << "Using xdl file: " << dataFnam << endl;
                break;
        case 'I':
                ASSERT_OR_THROW(sscanf(optarg, "%s", insFnam) == 1);
                dppLog << "Using instrument file: " << insFnam << endl;
                break;
        case 'W':
                ASSERT_OR_THROW(sscanf(optarg, "%s", drwDir) == 1);
                dppLog << "Using directory: " << drwDir << endl;
                break;
	case 'l':
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
	// Set debugging level
	debugLevel = atoi(loggingLevel);

	dppLog << "===========================================================================" << endl;

	//
	// Open deal data file
	//
	
	ifstream is(dataFnam);
	check(is, dataFnam);

	//
	// Read deal file
	//

	//----------------------------------------------
        // Parse instrument description file
	//----------------------------------------------
	vpRoot = VPInsParEvalFile(
		(insFnam[0] != '\0' ? insFnam : dataFnam),
		vpResetBank,
		discCurveName);
	ASSERT_OR_THROW(vpRoot != NULL);




        advanceIs(is, "END");

	//----------------------------------------------
        // Get the market curves
	// Return true if basis curves are required
	//----------------------------------------------
	isBasis = BasisTWrapperRead(
			is,
			drwDir,
			'B',

			today,
			mktCurves,
			irVolDiag,
			irMrParam,
			irSmileParam,
			bsVolDiag,
			bsMrParam,
			bsSmileParam,
			irBsCorr);


	// Zero shift from value date to today for discount curve.
	//
	discZeroShift = mktCurves.ZeroShift(discCurveName);


	dppLog << "---------------------------------------------------------------------------" << endl;
	dppLog << "Before CET: irVolDiag:\n" << irVolDiag << endl;
	dppLog << "---------------------------------------------------------------------------" << endl;
 
	//----------------------------------------------
        // Perform Pricing
	//----------------------------------------------

	// Calibrate PirTree vol
        //
        vt.CET( mktCurves,
		irVolDiag,
		irMrParam,
		irSmileParam,
                vpResetBank);
 

	dppLog << "---------------------------------------------------------------------------" << endl;
	dppLog << "After CET: irVolDiag:\n" << irVolDiag << endl;
	dppLog << "---------------------------------------------------------------------------" << endl;
 
        //
        //
        //
        //
        // Free memory
        delete vpToolRoot;
        delete vpRoot;


    }
    catch (KFailure) {
	// Free memory
	delete vpToolRoot;
	delete vpRoot;


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
    -I <file>   Use <file> for instrument (but reads params from data file)\n\
    -l <dbg>    Specify debug levle\n\
\n", __PROGRAM__);
		fprintf(stdout, "    Created %s  %s\n",
			__DATE__,__TIME__);
		fflush(stdout);
		exit(0);
	}
}

