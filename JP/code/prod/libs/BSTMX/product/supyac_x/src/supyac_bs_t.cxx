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


#define	__PROGRAM__	"supyac_bs_t"

static	void	usage(char *s, ...);

#define	__DEBUG__
#undef	__DEBUG__

#define FIDR_LITERAL(x) FIDR_LITERAL_(x)
#define FIDR_LITERAL_(x) #x
const char* _VERSION_ = FIDR_LITERAL(U_VERSION);

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

	strcpy(loggingLevel, "1");

	//
        // Parse options
        //
        while ((c = getopt(argc, argv, "I:W:l:h:v")) != EOF)
        switch (c) {
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
        case 'v':
                cout << format("VERSION %s COMPILED on %s %s\n",
                             _VERSION_, __DATE__,  __TIME__);
                return(SUCCESS);
        default:
        case 'h':
                usage(NULL);
                break;
        }
        argv += optind - 1 ;
 

	if (*(++argv) != NULL)
		strcpy(dataFnam, *argv);


 
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

	cout << format("NPV:  %.12f", npv) << endl;
	//DppWriteToFile("price", npv);
	ofstream outFile("price");
	if(!outFile)
		cerr << "cannot open \"price\" for output\n";
 
	// Output in decimal format. Kapital has a problem of
	// ignoring the scientific format (e.g. e+6). 
	outFile << format("%.12f", npv);

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
    -I <file>   Use <file> for instrument (but reads params from data file)\n\
    -l <dbg>    Specify debug level\n\
    -v <dbg>    Print version\n\
\n", __PROGRAM__);
		fprintf(stdout, "    Created %s  %s\n",
			__DATE__,__TIME__);
		fflush(stdout);
		exit(0);
	}
}

