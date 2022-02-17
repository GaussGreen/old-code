/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	
 * Function:	
 * Author:	Th original model has been developped by 
 *	        Lionnel Pradier & London DR
 * 		This version has been adapted/modified C Daher.
 ***************************************************************/
//#include "dvtalib.h"
#include "kstdinc.h"    /* Standard definitions & error hadling */
#include "kutilios.h"

#define	_kmrntree_SRC
#include "kpirtree.h"
#include "kmodpar.h"

extern	"C" {
#include "date_sup.h"
#include "drloptio.h"		// Black 
#include "cashflow.h"
#include "stub.h"               // GTO_STUB_NONE 

#include "drlstr.h"		// Strtok
#include "drltime.h"
#include "drlio.h"
#include "drlmem.h"
#include "drlproc.h"
#include "drlts.h"

};

#define	__EDEV__	// Use express dev tool instead of dev 

const String	MASKLENGTH("24M");	// match CET and product dates up to 2y	

//--------------------------------------------------------------
// Data structure containing the CET (Calibration Enhancement Tool)
// information (calibrated benchmark options definitions).
//

class	KCetData {
public:
	TDate		todayDate;	// option start date
	TDate		exerDate;	// option exercise date
	TDate		matDate;	// swap maturity date
	KVector(TDate)	cfDates;	// cash flow dates
	KVector(double) cfAmounts;	// cash flow amounts
	int		isOption;	// 0=fwd, 1=put


	double		parYield;	// benchmark par yields
	double		bsPrice;	// BS option price
	double		bsVega;		// BS option vega (* 100)


	double		treePrice;	// Tree prices for the benchmarks 
	double		treePricePrev;	// Tree prices for previous iteration
	
	double		priceDiff;	// Differences between B&S and tree 
	double		priceDiffInVega;
	
	double		volPrev;	// Previous iteration vol
	double		volNext;	// Next iteration vol (for printing)

	double		fOVega;		// Price/Vol approximation of vega
};



//--------------------------------------------------------------
// This routine initializes a KCetData structure by computing
// the forward yield of the reference swap underlying the option,
// the BS price and vega, and the set of cash flows to be used by
// the tree for the valuation.
//

static	void
KCetDataSetOne(
	TDate todayDate,		// (I) todayDate
	TDate optExerDate,		// (I) option exercise date
	int isOption,			// (I) 0=fwd, 1=put
	TDate swapStartDate,		// (I) swap start date
	TDate swapMatDate,		// (I) swap maturity date
	int swapFreq,			// (I) swap frequency
	KDayCc swapDcc,		// (I) swap day count conv

	KZCurve &discCurve,		// (I) zero coupon curve
	double vol,		// (I) Black/yield volatility
	KVolType volType,               // (I) Yield or bp vol

	KCetData &cetData)		// (O) calculated cet data 
{
        static	char	routine[] = "KCetDataSetOne";
    
	// Vol tweak sizes for vega calculation.
	double          lognVegaTwk = 0.0010;
	double          normVegaTwk = 0.0001;

	TCashFlowList	*underCashFlows = NULL;
	TCashFlowList	*underPrincipals = NULL;
	TDateInterval	payInterval;
	double		tExp,
	                dvol,	// vega twk size
			pv, pvUp,
			principalPv, bondYtm, annuityPv;
	double		cpnRate;
	int		idxC;

	TCurve		*zcCurve = (TCurve*) discCurve;

   try {

	// 
	IF_FAILED_THROW( GtoFreq2TDateInterval(
		(long) swapFreq,
		&payInterval));

	// create annuity CFL 
	ASSERT_OR_THROW((underCashFlows = GtoMakeCFL(
		1e0,
		swapStartDate,
		&payInterval,
		swapMatDate,
		(long) swapDcc,
		GTO_STUB_SIMPLE_ACT_ACT,
		FALSE,
		0L, //GTO_PRESTART_ZERO_PAYMENT, GTO_SUBTRACT_INITIAL,
		GTO_BAD_DAY_NONE,
		GTO_BAD_DAY_NONE,
		"NONE")) != NULL);


	// create bond CFL 
	ASSERT_OR_THROW((underPrincipals = GtoMakeCFL(
		0e0,
		swapStartDate,
		&payInterval,
		swapMatDate,
		(long) swapDcc,
		GTO_STUB_SIMPLE_ACT_ACT,
		FALSE,
		GTO_ADD_FINAL | GTO_SUBTRACT_INITIAL,
		GTO_BAD_DAY_NONE,
		GTO_BAD_DAY_NONE,
		"NONE")) != NULL);


	// compute pv of annuity 
	IF_FAILED_THROW(GtoCashFlowPV(
		underCashFlows,
		zcCurve,
		discCurve.ZeroInterpType(),
		&annuityPv));


	// compute zc pv of bond 
	IF_FAILED_THROW(GtoCashFlowPV(
		underPrincipals,
		zcCurve,
		discCurve.ZeroInterpType(),
		&principalPv));


	bondYtm = - principalPv / annuityPv;


#ifdef	__DEBUG__
	dppLog << format( "%s: annuityPv=%lf principaPv=%lf " \
		"bondYtm = %lf\n", \
		routine, annuityPv, principalPv, bondYtm);
#endif


	// compte time to expiration 
	IF_FAILED_THROW(GtoDayCountFraction(
		todayDate,
		optExerDate,
		GTO_ACT_365F,
		&tExp));


	// ATM option
	cpnRate = bondYtm;

	// compute Black  & vega
	//
	// log normal 
	if (isOption) {

	    dvol = normVegaTwk;

	    if (volType == LOGVOL) {

		IF_FAILED_THROW(DrlBlack(
		    tExp,
		    bondYtm,
		    vol,
		    cpnRate,
		    "P",	// Put
		    "P",	// Premium
		    &pv));
		
		IF_FAILED_THROW(DrlBlack(
		    tExp,
		    bondYtm,
		    vol+dvol,
		    cpnRate,
		    "P",	// Put
		    "P",	// Premium
		    &pvUp));

	    } else if (volType == NORMVOL) {

		dvol = lognVegaTwk;

		IF_FAILED_THROW(DrlNormOption(
		    tExp,
		    bondYtm,
		    vol,
		    cpnRate,
		    "P",        // Put
		    "P",        // Premium
		    &pv));

		IF_FAILED_THROW(DrlNormOption(
		    tExp,
		    bondYtm,
		    vol+dvol,
		    cpnRate,
		    "P",        // Put
		    "P",        // Premium
		    &pvUp));
		
	    } else {

		throw KFailure("%s: Only lognormal and normal vol allowed.\n", routine);

	    }

	    cetData.bsPrice  = pv * annuityPv ;
	    cetData.bsVega   = (pvUp - pv) * annuityPv / dvol;
	    
	} else {

		cetData.bsPrice  = 0e0;
		cetData.bsVega   = 0e0;

	}

	cetData.parYield = bondYtm;
	cetData.isOption = isOption;
	cetData.exerDate = optExerDate;


	//
	// Store cash flows
	//
	for (idxC=0; idxC<underCashFlows->fNumItems; idxC++) {
	    if (!IS_ALMOST_ZERO(cpnRate*underCashFlows->fArray[idxC].fAmount)) {
		cetData.cfDates.insert(cetData.cfDates.end(),
			underCashFlows->fArray[idxC].fDate);
		cetData.cfAmounts.insert(cetData.cfAmounts.end(),
			(cpnRate)*underCashFlows->fArray[idxC].fAmount);
	    }
	}

	for (idxC=0; idxC<underPrincipals->fNumItems; idxC++) {
	    if (!IS_ALMOST_ZERO(underPrincipals->fArray[idxC].fAmount)) {
		cetData.cfDates.insert(cetData.cfDates.end(),
			underPrincipals->fArray[idxC].fDate);
		cetData.cfAmounts.insert(cetData.cfAmounts.end(),
			underPrincipals->fArray[idxC].fAmount);
	    }
	}



	/*GtoPrintCFL(underCashFlows,  routine);
	GtoPrintCFL(underPrincipals, routine);*/

	if(debugLevel > DEBUG_LEVEL_CET) {
	    dppLog << format( "%20s: %7.4f%% %d %8s (%10s->%10s) %10s :" \
		"BSVOL=%8.4f %3s PV=%12.8f%% VEGA=%12.8f%%\n", \
		routine, \
		cetData.parYield*1e2, \
		swapFreq, \
		DrlTDayCountPrint(NULL, ((long) swapDcc)), \
		GtoFormatDate(swapStartDate), \
		GtoFormatDate(swapMatDate), \
		GtoFormatDate(optExerDate), \
		vol*1e2,
		(isOption ? "OPT" : "FWD"),
		cetData.bsPrice*1e2, cetData.bsVega);
#ifdef	__DEBUG__
	    dppLog << "Cash Flows" << endl;
	    for (idxC=0; idxC<cetData.cfDates.size(); idxC++) {
		dppLog << format("\t%10s\t%14.8f", 
			DrlTDatePrint(NULL, cetData.cfDates[idxC]),
			cetData.cfAmounts[idxC]) << endl;
	    }
#endif
	}



	GtoFreeCFL(underCashFlows);
	GtoFreeCFL(underPrincipals);
    }
    catch (KFailure) {
	GtoFreeCFL(underCashFlows);
	GtoFreeCFL(underPrincipals);
	throw KFailure("%s: failed.\n", routine);
    }
}


//--------------------------------------------------------------
// This routine creates and initializes an array of KCetData 
// from a volatility diagonal (and the market curves).
// The isOption flag if FALSE computes the price of 
//

static	void
KCetDataSet(
	TDate	lastDate,		// (I) last tree date
	KMarketCurves& marketCurves,	// (I) curves and curve types
	KVolDiag& irVolDiag,		// (I) IR volatility data.
	int isOption,			// (I) 0=fwd, 1=put
	KVector(KCetData) &cetData)	// (O) array of cet data
{
static	char	routine[] = "KCetDataSet";
	int	status = FAILURE;

	int		numBM, idxB;		// instruments
	TDate		todayDate;
	TDate		optExerDate;		// option exercise date
	TDate		swapStartDate;		// swap start date
	TDate		swapMatDate;		// swap maturity date
	int		swapFreq;		// swap frequency
	KDayCc		swapDcc;		// swap day count conv
	double		bsVol;			// volatility
	KVolType        volType;                // yield or bp vol

   try {
	// Check vol OK
	//
	if (!irVolDiag.IsValid()) 
		throw KFailure("%s: invalid vol diag.\n", routine);


	todayDate = marketCurves.mToday;

	//
	// Find the smallest vol date >= lastDate
	//
	KVector(TDate)::iterator iterVolDate = 
		lower_bound(irVolDiag.mVolDates.begin(),
			    irVolDiag.mVolDates.end(),
			    lastDate);

	if (iterVolDate == irVolDiag.mVolDates.end())
		numBM = irVolDiag.mVolDates.size();
	else
		numBM = iterVolDate - irVolDiag.mVolDates.begin() + 1;

	// Get diffusion crve
	//
	KZCurve	&discCurve = marketCurves.GetDiffuse();

	// calc and load data
	cetData.clear();
	cetData.resize(numBM);

	for (idxB=0; idxB<numBM; idxB++)  {
		//
		// The following defaults are used.
		// exer = start
		//
		optExerDate = irVolDiag.mVolDates[idxB];
		swapStartDate = optExerDate;
		swapMatDate = swapStartDate + 
			KDateInterval(irVolDiag.mVolMats[idxB]);
		swapFreq = irVolDiag.mVolFreqs[idxB];
		swapDcc = GTO_B30_360;

		bsVol   = irVolDiag.mVolRates[idxB];
		volType = irVolDiag.mVolType; 

		cetData[idxB].todayDate = todayDate ;
		cetData[idxB].matDate = swapMatDate;

		// Initialize
		cetData[idxB].parYield = 0e0;
		cetData[idxB].bsPrice = 0e0;
		cetData[idxB].bsVega = 0e0;
		cetData[idxB].treePrice = 0e0;
		cetData[idxB].treePricePrev = 0e0;
		cetData[idxB].priceDiff = 0e0;
		cetData[idxB].priceDiffInVega = 0e0;
		cetData[idxB].volPrev = 0e0;
		cetData[idxB].volNext = 0e0;
		cetData[idxB].fOVega = 0e0;

		// Calc cash flows and option values
		KCetDataSetOne(
			todayDate,
			optExerDate,
			isOption,
			swapStartDate,
			swapMatDate,
			swapFreq,
			swapDcc,
			discCurve,
			bsVol,
			volType,
			cetData[idxB]);
	}



    }
    catch (KFailure) {
         throw KFailure ("%s: failed.\n", routine);
    }
}



//--------------------------------------------------------------
// Prints an array of CET data.
//

void
KCetDataPrint(
	ostream &os,			// (I) 
	KVector(KCetData) &cetData)	// (I) array of cet data
{
static	char	routine[] = "KCetDataPrint";
	int	numBM, idxB;


	numBM = cetData.size();

	os << "CET DATA:" << endl;
	os << " #  Opt Start  Expir Date  Mat Date    ( yrs) "
		" ParYield  BS Price BS Vega |"
		" TR Price TR Pr Old PriceDif  InVega |"
		"  Vol Prev Vol Nxt" << endl;


	for (idxB=0; idxB<numBM; idxB++)  {
		os << format(
		"%3d %10s %10s %10s (%5.2f) "
		"%8.4f %8.4f %8.4f | "
		"%8.4f %8.4f %8.4f %8.4f | "
		"%8.4f %8.4f\n",

		idxB+1,
		DrlTDatePrint(NULL, cetData[idxB].todayDate),
		DrlTDatePrint(NULL, cetData[idxB].exerDate),
		DrlTDatePrint(NULL, cetData[idxB].matDate),
		(cetData[idxB].matDate - cetData[idxB].exerDate)/365.25e0,

		cetData[idxB].parYield*1e2,
		cetData[idxB].bsPrice*1e2,
		cetData[idxB].bsVega*1e2,


		cetData[idxB].treePrice*1e2,
		cetData[idxB].treePricePrev*1e2,
		cetData[idxB].priceDiff*1e2,
		cetData[idxB].priceDiffInVega*1e2,
	
		cetData[idxB].volPrev*1e2,
		cetData[idxB].volNext*1e2,
		cetData[idxB].fOVega*1e2);

	}



}



//--------------------------------------------------------------
// Computes the tree prices of a series of CET benchmarks.
//

void
KPirTreeCalcCet(
	KPirTree &cetTree,		// (I) CET tree
	KResetBank &resetBank,		// (I) rate reset bank

	KVector(KCetData) &cetData)	// (I/O) array of cet data
{
static  char    routine[] = "KPirTreeCalcCet";

	int			numBM, idxB,		// instruments
				numCF, idxC;		// cashflows

	KVector(KTSlice*)	swapTs;		// to store swap value
	KTSlice			*valueTs;
	double			pv;

	TDate			currDate;

 try{

	//
	// Calibrate drift
	//
	cetTree.CalibrateDrift();

	String	discCurveName = cetTree.GetCurveName(KV_DIFF);
	numBM = cetData.size();
	swapTs.resize(numBM);

	//
	// Insert swap slices
	//
	for (idxB=0; idxB<numBM; idxB++) 
		swapTs[idxB] = NULL;


	
	// Roll-back loop
	//
	for (int tpIdx=cetTree.TPNum(); tpIdx >= 0; tpIdx--)
	{

	    // Update tree (probas)
            cetTree.Update(tpIdx);

	    // Get current date
	    currDate = cetTree.TPDateCurrent();


	    // For each instrument
	    //
	    for (idxB=0; idxB<numBM; idxB++) {
		valueTs = swapTs[idxB];

		// 
		// (1) Dev
		//
		if (swapTs[idxB] != NULL) {
			swapTs[idxB]->Dev(discCurveName);
#ifdef	__DEBUG__
			dppLog << format("%10s Dev    BM%03d     \n",
				DrlTDatePrint(NULL, currDate),
				idxB);
#endif
			//cetTree.TSlicePut(*swapTs[idxB], dppLog, FALSE);

		}



		// 
		// (2) Check cash flow dates
		//
		KVector(TDate)  &cfDates  = cetData[idxB].cfDates;
		KVector(double) &cfAmounts= cetData[idxB].cfAmounts;
		numCF = cfDates.size();
		for (idxC=0; idxC<numCF; idxC++) {
		    if (currDate == cfDates[idxC]) {
	
			// Check slice allocated
			//
			if (swapTs[idxB] == NULL) {
			    swapTs[idxB] = new KTSlice(
				cetTree,
				format("Benchmark %d", idxB),
				discCurveName);
			    (*swapTs[idxB]) = 0e0;
			}

			// Add cash flow
			(*swapTs[idxB]) += cfAmounts[idxC];

#ifdef	__DEBUG__
			dppLog << format("%10s Adding BM%03d CF %3d "
				"  %12.8f\n", DrlTDatePrint(NULL, currDate),
				idxB, idxC, cfAmounts[idxC]); 
#endif
			//cetTree.TSlicePut(*swapTs[idxB], dppLog, FALSE);

		    }
		}


		//
		// (3) if expiration, express DEV
		//
		if (currDate == cetData[idxB].exerDate) {


#ifdef	__DEBUG__
			dppLog << format("%10s Exerc  BM%03d \n",
				DrlTDatePrint(NULL, currDate),
				idxB);
#endif

			// Evaluate payoff of option
			//
			if (cetData[idxB].isOption) {
				cetTree.sliceScalarOper(
					(double*)(valueTs->mData),
					valueTs->GetSliceDim(),
					0.0,
					MAX);
			}
#if defined(__EDEV__)
			// Express DEV: get the state prices
			// multiply state pr * slice and take sum
			//

			cetTree.sliceUnaryOper((double*)(valueTs->mData),
					valueTs->GetSliceDim(),
					cetTree.GetStatePrice(currDate),
					MULT);

			cetTree.sliceSpecialOper((double*)(valueTs->mData),
                                	 valueTs->GetSliceDim(),
					 "sum",
					 &pv);

			cetData[idxB].treePrice = pv;

			// Free unused time slice

			delete swapTs[idxB];
			swapTs[idxB] = NULL;

			// Free the obsolete state price
			//
			cetTree.FreeStatePrice(currDate);
#endif


		}
	    }
	}


#if !defined(__EDEV__)

	for (idxB=0; idxB<numBM; idxB++) {
		cetData[idxB].treePrice = swapTs[idxB]->GetCenter();
		delete swapTs[idxB];
		swapTs[idxB] = NULL;
	}

#endif

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}




//--------------------------------------------------------------
// Master CET calibration routine.
//
// Adjust the input "irVolDiag" volatility diagonal to fit
// the corresponding benchmarks options using the CET algorithm.
//

void
KPirTreeCet(
	KPirTree& prodTree,		// (B) product tree
  	KMarketCurves &marketCurves,	// (I) curves and curve types
	KVolDiag& irVolDiag,		// (B) IR volatility data.
  	KMrParam& irMrParam,		// (I) IR mr data.
  	KSmileParam& irSmileParam,	// (I) IR skew data.
  	KResetBank &resetBank)		// (I) rate reset bank
{
static  char    routine[] = "KPirTreeCet";

	int	iterCet,				//
		iterCetMax = irSmileParam.mNumIter;	// num iter
	int	numBM, idxB, 		// instruments
		numCF, idxC;		// cash flows

	double	slopeFO;		// 1st order 1st derv in NRaphson
	double	MaxDiff=0e0;		// Max price diff between tree & BS

	KVector(KCetData) cetData;	// array of cet data

	KMrParam        irMrParamCet = irMrParam; // Initialize to be the
						  // same as product parameters

	KPirTree	cetTree;	// CET tree		

	int		nDim;


	KVector(TDate)		 factVolDates;		// vol dates
	KVector(KVector(double)) factVolCurves;	// factor vol curves

    try {

	//
	// Check if anything to do.
	//
	if (iterCetMax <= 0 ||
		IS_ALMOST_ZERO(irVolDiag.mVolRates[0]+1e0)) 
	{
		return ;
	}


	//
        // Adjust ppy for factor 1 and 2.
	// Only apply to timeline beyond MASKLENGTH (2y).
        //
        switch (irMrParamCet.mNumFact) {
        case 1:
                irMrParamCet.mPpy = 24;
                break;
        case 2:
                irMrParamCet.mPpy = 4;
                break;
        }


	TDate todayDate = marketCurves.Today();

	//
	// Initialize the tree ONCE in the CET
	//
	nDim  = irMrParamCet.mNumFact;
 
        // Calibrate IR spot vols in irDim dimensions
        DppCalibIRVolDiag(
                irMrParamCet,
                irSmileParam,
                nDim,
                irVolDiag,
                marketCurves.GetDiffuse(),
                factVolDates,
                factVolCurves);
 
 
        cetTree.Initialize(
                todayDate,
                irMrParamCet,
                irSmileParam,
                nDim,
                marketCurves.mCVTypes,
                marketCurves.mCV,
                marketCurves.mCVNames,
                marketCurves.mValueDates,
                factVolDates,
                factVolCurves,
                resetBank);


	//
	// Get the product critical dates before MASKLENGTH
	//
	KVector(TDate) prodDates = prodTree.GetTDatesRange(
				todayDate,
				todayDate + KDateInterval(MASKLENGTH.c_str()));


	// 
	// Insert product dates << MASKLENGTH
	//
	for (KVector(TDate)::iterator iter=prodDates.begin();
	     iter != prodDates.end(); ++iter)
		cetTree.Insert(*iter);



	//
	// Get the last tree date.  Only
	// compute CET benchmark vols up to the tree last date
	//
	TDate lastDate = prodTree.GetLastDate();

	//
	// (1) Set up the CET data
	//
	KCetDataSet(
		lastDate,
		marketCurves,
		irVolDiag,
		1,			// option
		cetData);

	numBM = cetData.size();

	//
	// Insert benchmark swap cash flow dates
	//
	String	&discCurveName = marketCurves.mCVNames[KV_DIFF];
	for (idxB=0; idxB<numBM; idxB++) 
	{
		// Insert cash flow dates
		//
		KVector(TDate) &cfDates = cetData[idxB].cfDates;
		numCF = cfDates.size();
		for (idxC=0; idxC<numCF; idxC++) {
			cetTree.Insert(cfDates[idxC]);
		}

		cetTree.Insert(cetData[idxB].exerDate);

#if defined(__EDEV__)
		// Insert exp date for express dev
		//
		cetTree.InsertStatePrice(cetData[idxB].exerDate);
#endif

	}


	//
	// Tree set up
	//
	cetTree.SetUpTimeline();
	cetTree.TreeSetUp();
	cetTree.CheckTreeValid();



	//
	// !!!DONE with tree Init
	// Now begin the CET iteration
	//



	// For debugging: dry run
	if (iterCetMax == 0) {
	    //
	    // calc tree option prices
	    //
	    KPirTreeCalcCet(
		cetTree,	// (I) tree
		resetBank,
		cetData);

	    MaxDiff = 0.;
	    for (idxB=0; idxB<numBM; idxB++)
	    {
		slopeFO  = cetData[idxB].treePrice;
		slopeFO /= irVolDiag.mVolRates[idxB];

		cetData[idxB].priceDiff  = cetData[idxB].bsPrice 
                                              - cetData[idxB].treePrice;
		cetData[idxB].fOVega  = slopeFO;
		cetData[idxB].volPrev = irVolDiag.mVolRates[idxB];

		/*// Update inputs market volatilities
		//
		irVolDiag.mVolRates[idxB] +=
				cetData[idxB].priceDiff / slopeFO;*/

		// Store last step
		cetData[idxB].treePricePrev = cetData[idxB].treePrice;
		cetData[idxB].priceDiffInVega  = cetData[idxB].priceDiff/
                                                   cetData[idxB].bsVega;
		MaxDiff = MAX(fabs(cetData[idxB].priceDiffInVega),MaxDiff);         
	    }

 	    if (debugLevel > DEBUG_LEVEL_CET) {
	    	dppLog << format("CET DATA: NO ITERATIONS") << endl;
		KCetDataPrint(dppLog, cetData);
	    }

	    return;
	}

	if (debugLevel > DEBUG_LEVEL_CET) {
		dppLog << "INITIAL CET DATA:" << endl;
		KCetDataPrint(dppLog, cetData);
	}

	//
	// (2) Perform the CET iteration
	//
	for (iterCet=1; iterCet<=iterCetMax; iterCet++)
	{
	    //
            // Calibrate IR spot vols in irDim dimensions
	    //
	    if (iterCet > 1)
	    {
	    	//
	    	// Free tree memories that depends on the factor vols,
	    	// such as tree limits and probabilities, for
	    	// next CET iteration.
	    	//
	    	cetTree.ClearTreeVolMem();
	  
		//
		// clear tree vol and recalibrate
		//
		factVolDates.clear();
		factVolCurves.clear();

        	DppCalibIRVolDiag(
                	irMrParamCet,
                	irSmileParam,
             	        nDim,
                	irVolDiag,
                	marketCurves.GetDiffuse(),
                	factVolDates,
                	factVolCurves);

		// Update tree vol and geometry and re-calibrate the drift
		//
		cetTree.UpdateFactorVol(factVolDates,
			      		factVolCurves);

	    } 


	    //
	    // (2.1) update tree option prices
	    //
	    KPirTreeCalcCet(
		cetTree,	// (I) tree
		resetBank,
		cetData);


	    //
	    // (2.2) Update target diffs and improve vols on used vols
	    //
	    MaxDiff = 0.;
	    for (idxB=0; idxB<numBM; idxB++)
	    {
		slopeFO  = cetData[idxB].treePrice;
		slopeFO /= irVolDiag.mVolRates[idxB];

		cetData[idxB].priceDiff  = cetData[idxB].bsPrice 
                                              - cetData[idxB].treePrice;
		cetData[idxB].fOVega  = slopeFO;
		cetData[idxB].volPrev = irVolDiag.mVolRates[idxB];

		// Update inputs market volatilities
		//
		irVolDiag.mVolRates[idxB] +=
				cetData[idxB].priceDiff / slopeFO;
		cetData[idxB].volNext = irVolDiag.mVolRates[idxB];

		// Store last step
		cetData[idxB].treePricePrev = cetData[idxB].treePrice;
		cetData[idxB].priceDiffInVega  = cetData[idxB].priceDiff/
                                                   cetData[idxB].bsVega;
		MaxDiff = MAX(fabs(cetData[idxB].priceDiffInVega),MaxDiff);         
	    }


	    //
	    // Printout (if required)
	    //
	    if (debugLevel > DEBUG_LEVEL_CET) {
	   	 dppLog << format("CET DATA: ITERATION %d MAXDIFF=%5.2f%%",
			iterCet, MaxDiff*1e2) << endl;
		KCetDataPrint(dppLog, cetData);
	    }

	}


	//
	// Compute the final tree vol 
	//
	factVolDates.clear();
	factVolCurves.clear();

        DppCalibIRVolDiag(
               	irMrParamCet,
               	irSmileParam,
                nDim,
               	irVolDiag,
               	marketCurves.GetDiffuse(),
               	factVolDates,
               	factVolCurves);


	prodTree.InitializeFactorVol(factVolDates,
			  	     factVolCurves);

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}



