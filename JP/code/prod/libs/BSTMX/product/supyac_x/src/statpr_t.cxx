/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	dinstrum.h
 * Function:	
 * Author:	C. Daher
 * Revision:	$Header$
 ***************************************************************/
#include "kstdinc.h"
#include <errno.h>
#include "kutilios.h"
#include "kmrntree.h"
#include "kvtspar.h"

extern	"C" {

#include "drlmem.h"
#include "drlio.h"
#include "drltime.h"
#include "drlgetop.h"		// getopt() 
};

		char	buf[32];


#define	__PROGRAM__	"test"

static	void	usage(char *s, ...);

#define	__DEBUG__
#undef	__DEBUG__

inline	void	printf(KTSlice &z)
{
	fprintf(stdout, "%-30s = %12.8f\n", z.GetSliceName().c_str(),
			z.GetCenter());
}



/*--------------------------------------------------------------
 * This test just computes the EV of a formula of the
 * tree normal variables X_i.
 */

int
main(int argc, char **argv)
{
static	char	routine[] = "main";
	int	status = 0;
	int	c;
	FILE 	*fp = NULL;
	char	inputFnam[256];


	KMrNTree	mrnTree;

        int	ppy;
	double	smoothFact;
        char	EoI[32];
        int	numStdev;

        int	numFact;
        double	*factMr = NULL;
        double	*factWeight = NULL;
        double	*factCorr = NULL;

	double	tMat = 1e0;

        TDate				todayDate;
        int				numVolDates;
	KVector(TDate)			volDates;
	KVector(KVector(double))	factVol;

	KMap(String,double)		constTable;

#define	MAXFORM	256
	char	form[MAXFORM][1024];
	TDate	formDate;
	int	idxF,			// fact
		idxD;			// vol dates
	int	idxS, nf;		// formulas

	int	tpIdx;
	double	value0[MAXFORM];
	double	value1[MAXFORM];

	//
	//


    try {


	GtoErrMsgOn();
	//GtoLoggingSet(DrlGlobFlagGet(DRL_PROCFLAG_ID_GTO_LOGGING));
	GtoErrMsgFilePointer(stderr);
	setvbuf(stdout, NULL, _IONBF, 0);


	debugLevel = 0;

	/*
	 * Parse options
	 */
	ppy = 4;
	strcpy(EoI, "E");
	numStdev = 5;
	numFact = 3;

	while ((c = getopt(argc, argv, "d:f:p:s:h?")) != EOF)
 	switch (c) {
	case 'd':
		ASSERT_OR_THROW(sscanf(optarg, "%d", &debugLevel) == 1);
		break;
	case 'f':
		ASSERT_OR_THROW(sscanf(optarg, "%d", &numFact) == 1);
		break;
	case 'p':
		ASSERT_OR_THROW(sscanf(optarg, "%d", &ppy) == 1);
		break;
	case 's':
		ASSERT_OR_THROW(sscanf(optarg, "%s", EoI) == 1);
		break;
	case '?':
	case 'h':
		usage(NULL);
		break;
	default:
		break;
	}
	argv += optind - 1 ;
	/*if (argv[1] == NULL) usage(NULL);*/

	if (*(++argv) != NULL)
		strcpy(inputFnam, *argv);
	else
		strcpy(inputFnam, "statpr_t.dat");


	//
	// Read data
	//
	{
		ifstream is(inputFnam);
		check(is, inputFnam);

		cout << format("INPUT:\n");

		ppy = getInt(is, "ppy");
		strcpy(EoI, getString(is, "EoI"));
		smoothFact = 0e0;
		numStdev = getInt(is, "numStdev");
		numFact = getInt(is, "numFact");

		cout << format("#ppy\n\t%d\n", ppy);
		cout << format("#EoI\n\t%s\n", EoI);
		cout << format("#numStdev\n\t%d\n", numStdev);
		cout << format("#numFact\n\t%d\n", numFact);


		//
		// Read factor data
		//
        	factMr = DrlDoubleVectAlloc(0, numFact-1);
        	factWeight = DrlDoubleVectAlloc(0, numFact-1);
		factCorr   = DrlDoubleVectAlloc(0, MAX(numFact*(numFact-1)/2-1,0));

		for (idxF=0; idxF<numFact; idxF++)
			factMr[idxF] = getDouble(is, "factMr");
		for (idxF=0; idxF<numFact; idxF++)
			factWeight[idxF] = getDouble(is, "factWeight");
		for (idxF=0; idxF<numFact*(numFact-1)/2; idxF++)
			factCorr[idxF] = getDouble(is, "factCorr");

		// Print 
		//
		cout << "#factMr" << endl;
		for (idxF=0; idxF<numFact; idxF++)
			cout << '\t' << factMr[idxF];
		cout << endl;

		cout << "#factWeight" << endl;
		for (idxF=0; idxF<numFact; idxF++)
			cout << '\t' << factWeight[idxF];
		cout << endl;

		cout << "#factCorr" << endl;
		for (idxF=0; idxF<numFact*(numFact-1)/2; idxF++)
			cout << '\t' << factCorr[idxF];
		cout << endl;

		//
		// Read  factor volatilities
		//

		todayDate = getTDate(is, "todayDate");
		cout << format("#todayDate\n\t%s\n",
			DrlTDatePrint(NULL, todayDate));

	        numVolDates = getInt(is, "numVolDates");

		volDates.resize(numVolDates);

		KVector(double)	empty;	//  We had to do this on NT !
		factVol.resize(numFact, empty);

		for (idxF=0; idxF<numFact; idxF++) 
			factVol[idxF].resize(numVolDates);

		for (idxD=0; idxD<numVolDates; idxD++) {
			volDates[idxD] = getTDate(is, "volDates");
			for (idxF=0; idxF<numFact; idxF++) {
				factVol[idxF][idxD] = getDouble(is, "factVol");
			}
		}

		cout << "#volDates" << endl;
		for (idxD=0; idxD<numVolDates; idxD++) {
			cout << format("\t%10s",
				DrlTDatePrint(NULL, volDates[idxD]));
			for (idxF=0; idxF<numFact; idxF++) {
				cout << '\t' << factVol[idxF][idxD];
			}
			cout << endl;
		}

		// Read date
		formDate = getTDate(is, "formDate");


		// Read formulas
		nf = getInt(is, "numFormula");
		for (idxS=0; idxS<=nf-1; idxS++) {
			strcpy(form[idxS], getString(is, "form"));
			//cout << format("%3d  %s  \n", idxS, form[idxS]);
		}



	}



	//mrnTree.SetDebug(debugLevel);

	//
	// Initialize the tree
	//
	mrnTree.Initialize(
		todayDate,
        	ppy,
		smoothFact,
        	EoI[0],
        	numStdev,
        	numFact,
        	factMr,
        	factWeight,
        	factCorr,
        	volDates,
        	factVol);

	// We need to do that
	mrnTree.MapCurveName(K_DEFAULT_NAME, K_DEFAULT_IDX);


	//
	// Insert critical dates
	//
	mrnTree.Insert(formDate);


	//
	// Setup tree timeline and calibrate tree
	//
	mrnTree.SetUpTimeline();
	mrnTree.Calibrate();



	//
	// Create slices
	//
	KTSlice	*probaTs = new KTSlice(mrnTree);

	KTSlice	*u = KTSliceNewVector(nf, mrnTree);
	for (idxS=0; idxS<=nf-1; idxS++) u[idxS].SetSliceName(form[idxS]);

	KTSlice	*v = KTSliceNewVector(numFact, mrnTree);
	for (idxF=0; idxF<=numFact-1; idxF++)
		v[idxF].SetSliceName(format("VAR%d", idxF+1));

	//
	// Forward 
	//
	mrnTree.TSliceSetCenter(*probaTs, 1e0);

	for (tpIdx = 0; tpIdx <= mrnTree.TPNum(); tpIdx++) {
#ifdef	__DEBUG__
		cout << format("%s:TP %4d  %4d %4d\n", routine,
			tpIdx, mrnTree.TPIdxCurrent(), mrnTree.TPNum());
#endif

		if (formDate == mrnTree.TPDateCurrent())
			break;

		//
		// Update tree proba
		//
		mrnTree.Update(tpIdx);


		mrnTree.TSliceFwd(*probaTs);

#ifdef	__DEBUG__
		cout << *probaTs;
#endif
	}




	//
	// Backward
	//
	for (tpIdx = mrnTree.TPNum(); tpIdx >= 0; tpIdx--) {
#ifdef	__DEBUG__
		cout << format("%s:TP %4d  %4d %4d\n", routine,
			tpIdx, mrnTree.TPIdxCurrent(), mrnTree.TPNum());
#endif

		//
		// Update tree proba
		//
		mrnTree.Update(tpIdx);

		//
		// EV slices
		//
		for (idxS=0; idxS<=nf-1; idxS++) 
			mrnTree.TSliceEv(u[idxS]);



		/*for (idxS=0; idxS<=nf-1; idxS++) 
			u[idxS].Dev("");*/

		// Get gaussian variables (X_i)
		// 
		for (idxF=0; idxF<numFact; idxF++) {
		    mrnTree.TSliceScalarOper(v[idxF], (double)idxF, STVAR);
		}



		//
		// Compute formulas
		//
		for (idxS=0; idxS<=nf-1; idxS++) {

		    if (formDate == mrnTree.TPDateCurrent()) {
			// Evaluate the formula 
			KTSliceParEval(
				mrnTree,
				numFact,
				v,
				form[idxS],
				constTable,
				&u[idxS]);


			// Compute prob*value
			mrnTree.TSliceSpecialOper(u[idxS], "COPRODUCT",
				&u[idxS], probaTs, &value0[idxS]);

		    }
		}


	}

	// Collect ev
	for (idxS=0; idxS<=nf-1; idxS++)  {
		value1[idxS] = u[idxS].GetCenter();
	}

	//
	// Print results (slice EVs).
	//
	cout << "OUTPUT:" << endl;
	for (idxS=0; idxS<=nf-1; idxS++)  {
		cout << format("%3d  v0=%12.8f v1=%12.8f  diff=%18.14f `%s'\n",
			idxS,
			value0[idxS],
			value1[idxS],
			value0[idxS]-value1[idxS],
			form[idxS]);
	}



	// Free allocated memory
	//

        DrlDoubleVectFree(factMr, 0, numFact-1);
        DrlDoubleVectFree(factWeight, 0, numFact-1);
        DrlDoubleVectFree(factCorr, 0, MAX(numFact*(numFact-1)/2-1, 0));


	status = SUCCESS;
    }
    catch (KFailure) {
	// Free allocated memory
	//

        DrlDoubleVectFree(factMr, 0, numFact-1);
        DrlDoubleVectFree(factWeight, 0, numFact-1);
        DrlDoubleVectFree(factCorr, 0, numFact*(numFact-1)/2-1);

	GtoErrMsg("%s: failure caugth.\n", routine);
	status = FAILURE;
    }

    return(status);
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
    %s [options] \n\
OPTIONS\n\
    -p<ppy>     use <ppy> for periods per year\n\
    -d<level>   set debug level (0-20)\n\
    -s<E|I>     equal or increasing time steps\n\
\n", __PROGRAM__);
		fprintf(stdout, "    Created %s  %s\n",
			__DATE__,__TIME__);
		fflush(stdout);
		exit(0);
	}
}




