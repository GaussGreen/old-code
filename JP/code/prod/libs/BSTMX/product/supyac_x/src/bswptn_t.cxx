//**************************************************************
//
//
//
//**************************************************************
//#include <cstdio>
				// DR C++ Library
#include "kstdinc.h"
#include "kutilios.h"
#include "kmktcrv.h"

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


#define	__PROGRAM__	"bswptn_t"

static	void	usage(char *s, ...);

#define	__DEBUG__
#undef	__DEBUG__



//--------------------------------------------------------------
// Reads and loads the model/market data.
// Basis zero curve is taken from file riskzero.dat
// Spread volatiltiy is read from the last section of deal file.
//
//

bool
BasisTWrapperRead_MOD(
	istream& is,			// (I) deal input stream
	char *pathDir,			// (I) directory to read from
	char drwType,			// (I) '2', 'B'asis

	KMarketCurves&	mktCurves,	// (O) curves and curve types
	KVolDiag&	irVolDiag,	// (O) IR volatility data.
	KMrParam&	irMrParam,	// (O) IR mr data.
	KSmileParam&	irSmileParam,	// (O) IR skew data.
	KVolDiag&	bsVolDiag,	// (O) Basis volatility data.
	KMrParam&	bsMrParam,	// (O) Basis mr data.
	KSmileParam&	bsSmileParam,	// (O) Basis skew data.
	double&		irBsCorr)	// (O) IR basis correlation.

{
static	char	routine[] = "BasisTWrapperRead_MOD";
		
	TDrWrapperData	*drWrapData = NULL;
	long		drw_type = DRI_DRW_TYPE2_2CURVES;

			// --- zero curves

			// --- vol indices
	char		idx1s[128], idx2s[128];
	KVolCalibIndex	volIdx1, volIdx2;


	KMrParam	treeMrParam;

	int		bsDim;
	char		bsType;
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
                DRI_DRW_TYPE2_3CURVES,
                &drWrapData));


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
	irMrParam.ReadDrw(is, drWrapData);

	//----------------------------------------------
	// Read extra basis parameters and fill structures.
	//----------------------------------------------

	// Read wrapper data 
	if (isBasis) {

		bsType = getChar(is, "basis type ('S'pd, 'P'ercent)");
		strcpy(bsLiborCVName, getString(is, "basis libor curve name"));
		strcpy(bsDiscCVName, getString(is, "basis disc curve name"));

		//
		// Fill market curve basis info
		//
		switch (toupper(bsType)) {
		case 'S':
			mktCurves.mBSType = SPREAD;
			break;
		case 'P':
			mktCurves.mBSType = PERCENT;
			break;
		default:
			throw KFailure("%s: bad basis type %c (only S or P).\n",
				routine, bsType);
		}
		mktCurves.mLiborCVName = string(bsLiborCVName);
		mktCurves.mBSDiscCVName = string(bsDiscCVName);

		bsSmileParam.ReadDrw(is, drWrapData);
		
		bsMrParam.ReadDrw(is, drWrapData);


		bsDelayShift = 0.0;

		if (bsDelayShift >= 0e0 && bsDelayShift <= 1e0)
			mktCurves.mBSDelayShift = bsDelayShift;
		else
			throw KFailure("%s: invalid basis delay shift (%f). "
				       "Must be within the range [0, 1].\n",
					routine, bsDelayShift);

		//
		// To allow SIMPLE/PAR stub in basis leg,
		// the delayShift can NOT be anything other than 0.
		//
		if (!IS_ALMOST_ZERO(bsDelayShift))
			throw KFailure("%s: basis deley shift (%f) != 0!\n",
					routine,
					 bsDelayShift);


		// Chech consistency
		if (bsMrParam.mPpy != irMrParam.mPpy)
			throw KFailure("%s: IR ppy (%d) != Basis ppy (%d).\n",
					routine, 
					irMrParam.mPpy, bsMrParam.mPpy);

		if (!IS_ALMOST_ZERO(bsMrParam.mSmoothFact-irMrParam.mSmoothFact))
			throw KFailure("%s: IR smooth factor (%f) != Basis smooth factor (%f).\n",
					routine, 
					irMrParam.mSmoothFact, bsMrParam.mSmoothFact);
 
		if (!IS_ALMOST_ZERO(bsMrParam.mNumStdevCut-irMrParam.mNumStdevCut)) 
			throw KFailure("%s: IR # of stdev (%f) != Basis # of stdev (%f).\n",
					routine, 
					irMrParam.mNumStdevCut, bsMrParam.mNumStdevCut);

	}



	if (!volIdx1.IsNil()) {

		//----------------------------------------------
		// IR Volatilities: interpolate at calib index
		//----------------------------------------------

		irVolDiag = KVolDiag(volIdx1, drWrapData);

	} else {
		//----------------------------------------------
		// "NIL" calibration: we set up the dates, but basence
		// of vol is represented by -1 (UGLY !!! MUST BE CHANGED)
		// $$$
		//----------------------------------------------
		KVolCalibIndex dummyVolIdx("3m");
		irVolDiag = KVolDiag(dummyVolIdx, drWrapData);
	    	irVolDiag = -1e0;
	}



	//----------------------------------------------
	// Basis volatilities: use same timeline as IR
	//----------------------------------------------
	bsVolDiag = irVolDiag;
	bsVolDiag = 0e0;

	// Basis volatility curve
	//
	if (isBasis) {
	    int    i;
	    int    numBSVols   = getInt(is, "number of spread vols");
	    TDate  *bsVolDates = new TDate[numBSVols];
	    double *bsVolRates = new double[numBSVols];

	    for (i=0; i<=numBSVols-1; i++)
	    {
		bsVolDates[i] = getTDate(is, "Spread Vol Date");	
		bsVolRates[i] = getDouble(is, "Spread Vol Rate")/100e0;
	    }

	    bsVolCurve = GtoMakeTCurve(mktCurves.mToday,
				       bsVolDates,
				       bsVolRates,
				       numBSVols,
				       (double) 4L,	// 3M
				       GTO_ACT_365F);	// ignored

	    ASSERT_OR_THROW(bsVolCurve != NULL);

	    // Basis vol curve available, interpolate
            for (i=0; i<bsVolDiag.mVolRates.size(); i++)
            {
		// interpolate basis vol on these benchmark dates
		IF_FAILED_THROW(GtoInterpRate(
				irVolDiag.mVolDates[i],
				bsVolCurve,
				GTO_LINEAR_INTERP,
				&bsVolInterp));

		bsVolDiag.mVolRates[i] = bsVolInterp;
		bsVolDiag.mVolMats[i]  = 1e0/(double)(bsVolCurve->fBasis);
		bsVolDiag.mVolFreqs[i] = 0;
	    }

	    delete [] bsVolDates;
	    delete [] bsVolRates;
	    GtoFreeTCurve(bsVolCurve);

	    irBsCorr = getDouble(is, "corr");	

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
	DppLogMsgSetFile("run.log");


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
	isBasis = BasisTWrapperRead_MOD(
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


	// Zero shift from value date to today for discount curve.
	//
	discZeroShift = mktCurves.ZeroShift(discCurveName);



	//----------------------------------------------
        // Perform Pricing
	//----------------------------------------------
	//
	// Initialize tree
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
 
	dppLog << "===========================================================================" << endl;

	//
	// Build instrument tree pricing tools
	//
        vpToolRoot = NewToolRecursive(*vpRoot, vt);
        ASSERT_OR_THROW(vpToolRoot != NULL);
        dppLog << "===========================================================================" << endl;


	//
	// Set tree time line
	//
	vt.SetUpTimeline();


	//----------------------------------------------
        // CET
        //----------------------------------------------
 
	KPirTreeCet(
		vt,
		mktCurves,
		irVolDiag,
		irMrParam,
		irSmileParam,
		vpResetBank);
 
	if (debugLevel > 0) {
		dppLog << "irVolDiag (AFTER CET):\n" << irVolDiag << endl;
	}


	//
	// Calibrate tree
	//
	vt.Calibrate();

	//
	// Initialize tool
	//
	vpToolRoot->Initialize();

	//
	// Rollback
	//
	for (int tpIdx=vt.TPNum(); tpIdx >= 0; tpIdx--) {
		vt.Update(tpIdx);
		vpToolRoot->Update();
	}


	//
        //
        //
        PutPV(*vpToolRoot, dppLog, 0);


	// Output price.
	// Forward value at value date.
	//
	double  npv = vpToolRoot->GetResults()["PV"] / discZeroShift;

	cout << format("NPV:  %15.8f", npv) << endl;
	//DppWriteToFile("price", npv);
	ofstream outFile("price");
	if(!outFile)
		cerr << "cannot open \"price\" for output\n";

	// Output in decimal format. Kapital has a problem of
	// ignoring the scientific format (e.g. e+6).
	outFile << format("%lf", npv);
		

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

