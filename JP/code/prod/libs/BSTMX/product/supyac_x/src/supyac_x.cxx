//**************************************************************
//
//
//
//**************************************************************
                // DR C++ Library
#include "kstdinc.h"
#include "kutilios.h"
#include "kmktcrv.h"

                // VTree Tools
#include "vpipars.h"
#include "vtlprice.h"

extern	"C" {
#include "date_sup.h"

#include "drlstr.h"     // Strtok
#include "drltime.h"
#include "drlio.h"
#include "drlmem.h"
#include "drlproc.h"
#include "drlts.h"
#include "drlstr.h"
#include "drlgetop.h"   // getopt()

#include "dritkwrp.h"   // TDrWrapper
};

#define SUP_MAX_NUM_ZC 3

extern  int      debugLevel;

#define __PROGRAM__ "supyac_x"

static  void    usage(char *s, ...);

#define __DEBUG__
#undef  __DEBUG__


#define FIDR_LITERAL(x) FIDR_LITERAL_(x)
#define FIDR_LITERAL_(x) #x
const char* _VERSION_ = FIDR_LITERAL(U_VERSION);

//--------------------------------------------------------------
//

using std::auto_ptr;

int
main(int argc, char **argv)

{
static	char			routine[] = __PROGRAM__;

    int             TMXonFlag; // TRUE: use TMX3, FALSE: use FIX3
    char            CalibIdx[128];    // calibration index.  This is used in
                                      // loading smile parameters.  
//    auto_ptr<Tmx3Object> myTmxObj (new Tmx3Object);
    
    char            drwDir[1024];
    char            dataFnam[1024];
    char            insFnam[1024];	// instrument filename
    int             c;

    char            loggingLevel[1024];	// logging level


    KMarketCurves   mktCurves;      // curves and curve types
    KVolDiag        irVolDiag;      // IR volatility data.
    KMrParam        irMrParam;      // IR mr data.
    KSmileParam     irSmileParam;   // IR skew data.
    KVolDiag        bsVolDiag;      // Basis volatility data.
    KMrParam        bsMrParam;      // Basis mr data.
    KSmileParam     bsSmileParam;   // Basis skew data.
    double          irBsCorr;       // IR basis correlation.



    bool            isBasis;
    double          discZeroShift;
    String          discCurveName;  // product discount curve name 

    SharedPointer<KVPInstr> instrRoot;

    KResetBank      vpResetBank;

    KMap(String, double)    results;  


    try {


    // Enable logging
    //
    DppLogMsgSetFile("run.log");

    //
    // Enable error logging to stdout
    //
    DppErrMsgSet(DPP_ERR_MSG_COUT);


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
    instrRoot = VPInsParEvalFile(
        (insFnam[0] != '\0' ? insFnam : dataFnam),
        vpResetBank,
        discCurveName);
    ASSERT_OR_THROW(instrRoot != NULL);


    advanceIs(is, "END");


    //----------------------------------------------
    // Get the market curves
    // Return true if basis curves are required
    //----------------------------------------------
    isBasis = TmxBasisTWrapperRead(
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
            irBsCorr,
            &TMXonFlag,
            CalibIdx);


    // ZeroShift
    discZeroShift = mktCurves.ZeroShift(discCurveName);

    cout << "CalibIndex " << CalibIdx << endl;
    if (TMXonFlag == TRUE)
    {
        cout << "TMX" << endl;

        KVPRootPrice_TMX(   results,
                            instrRoot,
                            mktCurves,
                            irVolDiag,
                            irMrParam,
                            irSmileParam,
                            bsVolDiag,
                            bsMrParam,
                            bsSmileParam,
                            irBsCorr,
                            vpResetBank,
                            debugLevel);
    }
    //----------------------------------------------
    // Perform Pricing
    //----------------------------------------------
    else
    {
        cout << "FIX" << endl;
 
        KVPRootPrice(
            results,
            instrRoot,
            mktCurves,
            irVolDiag,
            irMrParam,
            irSmileParam,
            bsVolDiag,
            bsMrParam,
            bsSmileParam,
            irBsCorr,
            vpResetBank,
            debugLevel);
    }


        //
        //
        //PutPV(vpToolRoot, dppLog, 0);


    // Output price.
    // Forward value at value date.
    //
    double  npv = results["PV"] / discZeroShift;

    cout << format("NPV:  %15.8f", npv) << endl;
    ofstream outFile("price");
    if(!outFile)
        cerr << "cannot open \"price\" for output\n";
 
    // Output in decimal format. Kapital has a problem of
    // ignoring the scientific format (e.g. e+6). 
    outFile << format("%lf", npv);


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

