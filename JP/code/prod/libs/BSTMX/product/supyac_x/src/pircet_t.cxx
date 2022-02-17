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
				// Basis Tree
#include "kbirtree.h"

extern	"C" {
#include "drlstr.h"		// Strtok
#include "drltime.h"
#include "drlio.h"
#include "drlmem.h"
#include "drlproc.h"
#include "drlts.h"
#include "drlstr.h"
#include "drlgetop.h"		// getopt()

};


extern 	int      debugLevel;
static	void	usage(char *s, ...);


#define	__PROGRAM__	"pircet_t"



#define	__DEBUG__
#undef	__DEBUG__


//--------------------------------------------------------------
// Compares and prints the forward rates computed by the tree
// and the implied tree BS smile.
//

int
main(int argc, char **argv)
{
static	char			routine[] = __PROGRAM__;

	char			drwDir[1024];
	char			dealFnam[1024];
	int                     c;

	char			loggingLevel[1024];	// logging level

	TDate           today;          // today's date
	KMarketCurves   mktCurves;      // curves and curve types
	KVolDiag        irVolDiag;      // IR volatility data.
	KMrParam        irMrParam;      // IR mr data.
	KSmileParam     irSmileParam;   // IR skew data.
	KVolDiag        bsVolDiag;      // Basis volatility data.
	KMrParam        bsMrParam;      // Basis mr data.
	KSmileParam     bsSmileParam;   // Basis skew data.
	double          irBsCorr;       // IR basis correlation.
 
	KResetBank	resetBank;
 
	bool		isBasis;	

	KBirTree	vt;

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
        while ((c = getopt(argc, argv, "W:d:h")) != EOF)
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
	//debugLevel = 13;

	dppLog << "============================================="
		  "==============================" << endl;

	//
	// Open deal data file
	//
	
	ifstream is(dealFnam);
	check(is, dealFnam);


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


	cout << "INITIAL VOL" << endl;
	cout << irVolDiag;

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
                resetBank);

	//
	// Need a last product date on the tree to tell
	// CET to compute all the vol points
	// 
	vt.Insert(irVolDiag.mVolDates.back());

	vt.SetUpTimeline();

	KPirTreeCet(
		vt,
		mktCurves,
		irVolDiag,
		irMrParam,
		irSmileParam,
		resetBank);

	cout << "FINAL VOL" << endl;
	cout << irVolDiag;


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

