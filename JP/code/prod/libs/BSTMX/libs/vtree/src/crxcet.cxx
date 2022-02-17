/***************************************************************
 * Module:	Credit CET calibration module for credit tree
 * Submodule:	
 * File:	crxcet.cxx
 * Function:	
 * Author:	Charles Morcom, adapted and extended from IR CET code 
 *          by Lionnel Pradier, Christian Daher, and David Liu
 ***************************************************************/
//#include "dvtalib.h"
#include "kstdinc.h"    /* Standard definitions & error handling */
#include "kutilios.h"

#define	_prcet_SRC
#include "kcrxtree.h"
#include "kmodpar.h"




extern	"C" {
#include "date_sup.h"
#include "drloptio.h"		// Black 
#include "cashflow.h"
#include "stub.h"               // GTO_STUB_NONE 
#include "crcrv.h" // crxflow protection and fee legs for BM CDS

#include "drlstr.h"		// Strtok
#include "drltime.h"
#include "drlio.h"
#include "drlmem.h"
#include "drlproc.h"
#include "drlts.h"
    #include "stdio.h"

};

#define	__EDEV__	// Use express dev tool instead of dev 

const String	MASKLENGTH("24M");	// match CET and product dates up to 2y	

//--------------------------------------------------------------
// Data structure containing the CET (Calibration Enhancement Tool)
// information (calibrated benchmark options definitions).
//

class	KCrxCetData {
public:
	TDate		todayDate;	// option start date
	TDate		exerDate;	// option exercise date
	TDate		matDate;	// cds maturity date
	KVector(TDate)	cfDates;	// cash flow dates
	KVector(double) cfAmounts;	// cash flow amounts
    double      recovery;   // recovery R for protection leg (1-R)
    KZeroReset  protZero; // zero reset for protection leg calc.
	int		isOption;	// 0=fwd, 1=put

	double		parSpread;	// benchmark par spread
	double		bsPrice;	// BS option price
	double		bsVega;		// BS option vega (* 100)
    double      bsVol;      // BM Black-Scholes implied vol.

	double		treePrice;	// Tree prices for the benchmarks 
	double		treePricePrev;	// Tree prices for previous iteration
	
	double		priceDiff;	// Differences between B&S and tree 
	double		priceDiffInVega;
	
	double		volPrev;	// Previous iteration vol
	double		volNext;	// Next iteration vol (for printing)

	double		fOVega;		// Price/Vol approximation of vega
    KCrxCetData() {}; // MSVCC needs this for vector, but shouldn't be necessary
};



//--------------------------------------------------------------
// This routine initializes a KCrxCetData structure by computing
// the forward par spread of the reference CDS underlying the 
// option,
// the BS price and vega, and the set of cash flows to be used by
// the tree for the valuation.
//

static	void
KCrxCetDataSetOne(
	TDate todayDate,		// (I) todayDate
	TDate optExerDate,		// (I) option exercise date
	int isOption,			// (I)4 0=fwd, 1=put
	TDate cdsStartDate,		// (I) swap start date
	TDate cdsMatDate,		// (I) swap maturity date
	int cdsFreq,			// (I) swap frequency
	KDayCc cdsDcc,		    // (I) swap day count conv
	KZCurve &discCurve,		// (I) zero coupon discount curve
    KZCurve &credCurve,     // (I) credit default curve
    double recoveryRate,    // (I) CDS recovery rate
	double vol,		// (I) Black/spread volatility
	KVolType volType,               // (I) Spread or bp vol
	KCrxCetData &cetData)		// (O) calculated cet data 
{
    static	char	routine[] = "KCrxCetDataSetOne";
    
	// Vol tweak sizes for vega calculation.
	double          lognVegaTwk = 0.0010;
	double          normVegaTwk = 0.0001;
    double discountToPayoff;

	TDateInterval	payInterval, integrate3M, zeroDelay;
	double		tExp,
	                dvol,	// vega twk size
			pv, pvUp,
			protectionPv, annuityPv;
	double		fwdSpread, spdStrike;
	int		idxC;

	TCurve		*zcCurve = (TCurve*) discCurve;
    TCurve      *crCurve = (TCurve*) credCurve;
    KFeeLeg_D* feeLeg = 0; // fee leg
    KProtLeg_D* protLeg = 0; // protection leg
    
    try {

	    // 
	    IF_FAILED_THROW( GtoFreq2TDateInterval(
		    (long) cdsFreq,
		    &payInterval));

        // 3M integration frequency for CDS protection
        IF_FAILED_THROW(
            GtoMakeDateInterval(
		    1,'Q',
		    &integrate3M));

        // zero payment delay on default
        IF_FAILED_THROW(
            GtoMakeDateInterval(
		    0,'Q',
		    &zeroDelay));

        /* Create fee leg for fwd annuity value and cash-flows */
        feeLeg = CrxFeeLegCreateFromFreq(
            cdsStartDate,       
            cdsMatDate,
            payInterval,
            SHORT_FRONT, /* Stub Location */
            1.0, /* Notional */
            1.0, /* Coupon */
            (long)cdsDcc, // converts to long, then implicitly to ALIB
            ACCRUAL_PAY_NONE, /* All accrued on default */
	        integrate3M);
        if (feeLeg==0) {
            throw KFailure("%s: Failed to create fee leg.\n", routine);
        }

        // compute fv of fee annuity
        IF_FAILED_THROW(RiskyFeePV_O(
            &annuityPv,
            optExerDate,
	        cdsStartDate,
            feeLeg,
            SIMPLE, /* Accrual type */
            zcCurve,
            crCurve
        ));

        double protNtl[1] = {1.0};
        /* Create protection leg */
        protLeg = ProtectionCreate(
            cdsStartDate,
            1, // only one protection period
            &cdsMatDate,
            protNtl,
            &recoveryRate,
            PAY_DEF, /* Pay protection on default */
	        zeroDelay, /* Delay from default to payment */
	        integrate3M
        );

        // compute fv of protection leg
    
        IF_FAILED_THROW(ProtectionPV_O(
            &protectionPv,
            optExerDate,
            cdsStartDate,
            protLeg,
            zcCurve,
            crCurve
        ));

	    fwdSpread = protectionPv / annuityPv;


        #ifdef	__DEBUG__
	        dppLog << format( "%s: feePv=%lf protectionPv=%lf " \
		        "fwdSpread = %lf\n", \
		        routine, feePv, protectionPv, fwdSpread);
        #endif


	    // compte time to expiration 
	    IF_FAILED_THROW(GtoDayCountFraction(
		    todayDate,
		    optExerDate,
		    GTO_ACT_365F,
		    &tExp));


	    // ATM option
	    spdStrike = fwdSpread;

        discountToPayoff = RiskyDiscountFactor(todayDate+1,optExerDate,zcCurve,crCurve);

	    // compute Black  & vega
	    if (isOption) {

	        dvol = normVegaTwk;

	        if (volType == LOGVOL) {

		        IF_FAILED_THROW(DrlBlack(
		            tExp,
		            fwdSpread,
		            vol,
		            spdStrike,
		            "P",	// Put
		            "P",	// Premium
		            &pv));
		        
		        IF_FAILED_THROW(DrlBlack(
		            tExp,
		            fwdSpread,
		            vol+dvol,
		            spdStrike,
		            "P",	// Put
		            "P",	// Premium
		            &pvUp));

	        } else if (volType == NORMVOL) {

		        dvol = lognVegaTwk;

		        IF_FAILED_THROW(DrlNormOption(
		            tExp,
		            fwdSpread,
		            vol,
		            spdStrike,
		            "P",        // Put
		            "P",        // Premium
		            &pv));

		        IF_FAILED_THROW(DrlNormOption(
		            tExp,
		            fwdSpread,
		            vol+dvol,
		            spdStrike,
		            "P",        // Put
		            "P",        // Premium
		            &pvUp));
		    
	        } else {

		        throw KFailure("%s: Only lognormal and normal vol allowed.\n", routine);

	        }

	        cetData.bsPrice  = discountToPayoff * pv * annuityPv ;
	        cetData.bsVega   = discountToPayoff * (pvUp - pv) * annuityPv / dvol;
	        
	    } else {

		    cetData.bsPrice  = 0e0;
		    cetData.bsVega   = 0e0;

	    }

    cetData.bsVol = vol;
	cetData.parSpread = fwdSpread;
	cetData.isOption = isOption;
	cetData.exerDate = optExerDate;
    cetData.recovery = recoveryRate;
	//
	// Store cash flows in cetData structure
	//
	for (idxC=0; idxC<feeLeg->mNbCF; idxC++) {
        double dcf;
        GtoDayCountFraction(feeLeg->mAccStDates[idxC], feeLeg->mAccEndDates[idxC], feeLeg->mDCC, &dcf);

	    if (!IS_ALMOST_ZERO(fwdSpread*feeLeg->mCoupons[idxC])) {
		    cetData.cfDates.insert(cetData.cfDates.end(),
			    feeLeg->mPayDates[idxC]);
		    cetData.cfAmounts.insert(cetData.cfAmounts.end(),
			    fwdSpread * feeLeg->mCoupons[idxC]*dcf);
	    }
	}

	if(debugLevel > DEBUG_LEVEL_CET) {
	    dppLog << format( "%20s: %7.4f%% %d %8s (%10s->%10s) %10s :" \
		"BSVOL=%8.4f %3s PV=%12.8f%% VEGA=%12.8f%%\n", \
		routine, \
		cetData.parSpread*1e4, \
		cdsFreq, \
		DrlTDayCountPrint(NULL, ((long) cdsDcc)), \
		GtoFormatDate(cdsStartDate), \
		GtoFormatDate(cdsMatDate), \
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

        CrxProtectionFree(protLeg);
        CrxFeeLegFree(feeLeg);

    } catch (KFailure) {

        CrxProtectionFree(protLeg);
        CrxFeeLegFree(feeLeg);
	    throw KFailure("%s: failed.\n", routine);
    }

    
} /* KCrxCetDataSetOne */


//--------------------------------------------------------------
// This routine creates and initializes an array of KCetData 
// from a volatility diagonal (and the market curves).
// The isOption flag if FALSE computes the price of 
//

static	void
KCrxCetDataSet(
	TDate	lastDate,		// (I) last tree date
	KMarketCurves& marketCurves,	// (I) curves and curve types
	KVolDiag& crVolDiag,		// (I) IR volatility data.
	int isOption,			// (I) 0=fwd, 1=put
	KVector(KCrxCetData) &cetData)	// (O) array of cet data
{
static	char	routine[] = "KCrxCetDataSet";
	int	status = FAILURE;

	int		numBM, idxB;		// instruments
	TDate		todayDate;
	TDate		optExerDate;		// option exercise date
	TDate		cdsStartDate;		// swap start date
	TDate		cdsMatDate;		// swap maturity date
	int		cdsFreq;		// swap frequency
	KDayCc		cdsDcc;		// swap day count conv
	double		bsVol;			// volatility
	KVolType        volType;                // yield or bp vol

   try {
	// Check vol OK
	//
	if (!crVolDiag.IsValid()) 
		throw KFailure("%s: invalid vol diag.\n", routine);

	todayDate = marketCurves.mToday;

	//
	// Find the smallest vol date >= lastDate
	//
	KVector(TDate)::iterator iterVolDate = 
		lower_bound(crVolDiag.mVolDates.begin(),
			    crVolDiag.mVolDates.end(),
			    lastDate);

	if (iterVolDate == crVolDiag.mVolDates.end())
		numBM = crVolDiag.mVolDates.size();
	else
		numBM = iterVolDate - crVolDiag.mVolDates.begin() + 1;

	// Get diffusion crve
    KZCurve	&discCurve = marketCurves.mCV[marketCurves.GetCurveType(marketCurves.mIRDiscCVName)];
    KZCurve &credCurve = marketCurves.mCV[KV_CREDIT_RISKY];
    double recovery = marketCurves.mRecovery;

	// calc and load data
	cetData.clear();
	cetData.resize(numBM);

	for (idxB=0; idxB<numBM; idxB++)  {
		//
		// The following defaults are used.
		// exer = start
		//
		optExerDate = crVolDiag.mVolDates[idxB];
		cdsStartDate = optExerDate;
		cdsMatDate = cdsStartDate + 
			KDateInterval(crVolDiag.mVolMats[idxB]);
		cdsFreq = crVolDiag.mVolFreqs[idxB];
		cdsDcc = GTO_ACT_360;

		bsVol   = crVolDiag.mVolRates[idxB];
		volType = crVolDiag.mVolType; 

		cetData[idxB].todayDate = todayDate ;
		cetData[idxB].matDate = cdsMatDate;

		// Initialize
		cetData[idxB].parSpread = 0e0;
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
		KCrxCetDataSetOne(
			todayDate,
			optExerDate,
			isOption,
			cdsStartDate,
			cdsMatDate,
			cdsFreq,
			cdsDcc,
			discCurve,
            credCurve,
            recovery,
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
KCrxCetDataPrint(
	ostream &os,			// (I) 
	KVector(KCrxCetData) &cetData)	// (I) array of cet data
{
static	char	routine[] = "KCrxCetDataPrint";
	int	numBM, idxB;


	numBM = cetData.size();

	os << "CET CREDIT DATA:" << endl;
	os << " #  Opt Start  Expir Date  Mat Date    ( yrs) "
		" ParSpread  BS Price BS Vega |"
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

		cetData[idxB].parSpread*1e4,
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
// Computes the tree prices of a series of Credit CET benchmarks.
// NOTE THAT THIS CALIBRATES TO CDS OPTIONS WITH NO DEFAULT ACCRUAL
// IN PRACTICE, IT MAKES VERY LITTLE DIFFERENCE, AND IS SIGNIFICANTLY
// FASTER.
void
KCrxTreeCalcCet(
	KCrxTree &cetTree,		// (I) CET tree
	KResetBank &resetBank,		// (I) rate reset bank
	KVector(KCrxCetData) &cetData)	// (I/O) array of cet data
{
    static  char    routine[] = "KCrxTreeCalcCet";

	int			numBM, idxB,		// instruments
				numCF, idxC;		// cashflows

	KVector(KTSlice*)	feeTs;		// to store swap fee leg value
    KVector(KTSlice*)   protTs;     // protection leg
	KTSlice			*valueTs;       // current slice  
    

	double			pv;
	TDate			currDate;


    try {

	    //
	    // Calibrate drift
	    //
	    cetTree.CalibrateDrift();
        

	    String discCurveName = cetTree.GetCurveName(KV_CREDIT_RISKY); // want Z[r]Z[l] discounting in DEV
        String protCurveName = cetTree.GetCurveName(KV_CREDIT_PROT_DEFRECOV); // protection leg "curve"
        String cpnDefCurveName = cetTree.GetCurveName(KV_CREDIT_PROT_BINRECOV); // recover 100%
	    numBM = cetData.size();
	    feeTs.resize(numBM);
        protTs.resize(numBM);

	    //
	    // Initialize all swap slices to NULL
	    //
        for (idxB=0; idxB<numBM; idxB++) {
            feeTs[idxB] = NULL;
            protTs[idxB] = NULL;
        }

	
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

		        // (1) Dev to current node's date, if already created
		        if (feeTs[idxB] != NULL) {
                    feeTs[idxB]->Dev(discCurveName); // risky discounting
                    #ifdef	__DEBUG__
			                    dppLog << format("%10s Dev    BMFee%03d     \n",
				                    DrlTDatePrint(NULL, currDate),
				                    idxB);
                    #endif

                }
                if (protTs[idxB] != NULL) {
                    protTs[idxB]->Dev(protCurveName); // protection leg discounting
                    #ifdef	__DEBUG__
			                    dppLog << format("%10s Dev    BMProt%03d     \n",
				                    DrlTDatePrint(NULL, currDate),
				                    idxB);
                    #endif

                } 

                /* Create protection leg slice, if at maturity date */
                if (currDate == cetData[idxB].matDate) {
                    if (protTs[idxB] == NULL) 
                        protTs[idxB] = new KTSlice(cetTree, format("BM Protection %d", idxB), protCurveName);
                   *protTs[idxB] = 0.0; // start protection leg at zero value at maturity
                    //*protTs[idxB] *= (1. - cetData[idxB].recovery);    
                }

		        // (2) Check cash flow dates
		        KVector(TDate)  &cfDates  = cetData[idxB].cfDates;      // dates
		        KVector(double) &cfAmounts= cetData[idxB].cfAmounts;    // payments
		        numCF = cfDates.size();

                // add cash-flows at time points where they happen and create slices
		        for (idxC=0; idxC<numCF; idxC++) {
		            if (currDate == cfDates[idxC]) {
			            // Check slice allocated - if not, create it
			            if (feeTs[idxB] == NULL) {
			                feeTs[idxB] = new KTSlice(
				                cetTree,
				                format("BM Fee %d", idxB),
				                discCurveName);
			                (*feeTs[idxB]) = 0e0; // initial value is zero before add flow below
			            }

			            // Add cash flow at all nodes on payment date
			            (*feeTs[idxB]) += cfAmounts[idxC];

                        #ifdef	__DEBUG__
			                        dppLog << format("%10s Adding BM%03d CF %3d "
				                        "  %12.8f\n", DrlTDatePrint(NULL, currDate),
				                        idxB, idxC, cfAmounts[idxC]); 
                        #endif
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

                    valueTs = feeTs[idxB]; 
                    *valueTs -= *protTs[idxB];

			        // Evaluate payoff of option (call payoff - max(V,0))
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
                        // this discounts each node to now
			            cetTree.sliceUnaryOper((double*)(valueTs->mData),
					           valueTs->GetSliceDim(),
					            cetTree.GetStatePrice(currDate), // N.B. these are risky
					            MULT);

			            cetTree.sliceSpecialOper(
                            (double*)(valueTs->mData),
                             valueTs->GetSliceDim(),
					         "sum",
					         &pv);
			            cetData[idxB].treePrice = pv;

			            // Free unused time slice

			            delete feeTs[idxB];
			            feeTs[idxB] = NULL;
                        delete protTs[idxB];
                        protTs[idxB] = NULL;

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

    } catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
} /* KCrxTreeCalcCet */


//--------------------------------------------------------------
// Master CET calibration routine.
//
// Adjust the input "crVolDiag" volatility diagonal to fit
// the corresponding benchmarks options using the CET algorithm.
// this assumes that the IR vol has already been calibrated

void
KCrxTreeCet(
	KCrxTree& prodTree,		// (B) product tree
  	KMarketCurves &marketCurves,	// (I) curves and curve types
    KVolDiag& irVolDiag,        // (I) Calibrated IR vol curve
    KMrParam& irMrParam,		// (I) IR mr data.
  	KSmileParam& irSmileParam,	// (I) IR skew data.
    KVolDiag& crVolDiag,		// (B) CR volatility data.
  	KMrParam& crMrParam,		// (I) CR mr data.
  	KSmileParam& crSmileParam,	// (I) CR skew data.
    double irCrCorr,            // (I) credit/rates correlation
  	KResetBank &resetBank)		// (I) rate reset bank
{
    static char routine[] = "KCrxTreeCet";
	int	iterCet,				//
		iterCetMax = crSmileParam.mNumIter;	// num iter
	int	numBM, idxB, 		// instruments
		numCF, idxC;		// cash flows
	double	MaxDiff=0e0;		// Max price diff between tree & BS
	KVector(KCrxCetData) cetData;	// array of cet data
	KMrParam        crMrParamCet = crMrParam; // Initialize to be the
						  // same as product parameters
	KCrxTree	cetTree;	// CET tree		
	int		nDim, crDim, irDim;
	KVector(TDate)		 factVolDates;		// vol dates
	KVector(KVector(double)) factVolCurves;	// factor vol curves
    KVector(KVector(double)) vegaDifferences;
    int idx;
    //FILE* el;

    try {

	    if (iterCetMax <= 0 || IS_ALMOST_ZERO(crVolDiag.mVolRates[0]+1e0)) {
            /*=================================================================
             * NO ITERATIONS - prodTree INITIALIZED WITH VNFM OK OR LOW VOL
             *===============================================================*/
		    return ;
        } else if (!marketCurves.IsCredit()) {
            /*=================================================================
             * NO CREDIT DIFFUSION - NOTHING TO DO
             *===============================================================*/
            return;
        }

        int i;

	    //
        // Adjust ppy for factor 1 and 2.
	    // Only apply to timeline beyond MASKLENGTH (2y).
        // This is same as IR CET - check to see what is suitable for credit
        // clearly doesn't matter, yet, since only 1 factor possible
        switch (crMrParamCet.mNumFact) {
        case 1:
                crMrParamCet.mPpy = 24;
                break;
        case 2:
                crMrParamCet.mPpy = 4;
                break;
        }

	    TDate todayDate = marketCurves.Today();

	    

        /*=====================================================================
         * INITIALIZE THE cetTree AND DO IR AND CR VNFM APPROX AS APPROPRIATE 
         *===================================================================*/
        crDim = crMrParam.mNumFact;
        irDim = irMrParam.mNumFact;
        nDim = crDim + irDim;
        cetTree.Initialize(
                marketCurves,     // (I) curves and curve types
                irVolDiag,        // (I) IR volatility data.
                irMrParam,        // (I) IR mr data.
                irSmileParam,     // (I) IR skew data.
                crVolDiag,        // (I) Credit volatility data.
                crMrParam,        // (I) Credit mr data.
                crSmileParam,     // (I) Credit skew data.
                irCrCorr,         // (I) IR crdit correlation.
                resetBank);       // (I) rate reset bank

        /*=====================================================================
         * SET-UP TREE CRITICAL DATES, TIMES AND BM SWAPTIONS
         *===================================================================*/
	    // Get the product critical dates before MASKLENGTH
	    KVector(TDate) prodDates = prodTree.GetTDatesRange(
				    todayDate,
				    todayDate + KDateInterval(MASKLENGTH.c_str()));

	    // Insert product dates << MASKLENGTH
	    for (KVector(TDate)::iterator iter=prodDates.begin();
	         iter != prodDates.end(); ++iter)
		    cetTree.Insert(*iter);

	    // Get the last tree date.  Only
	    // compute CET benchmark vols up to the tree last date
	    TDate lastDate = prodTree.GetLastDate();

	    // (1) Set up the CET benchmark data and clc target swaption prices
	    KCrxCetDataSet(
		    lastDate,
		    marketCurves,
		    crVolDiag,
		    1,			// options, not forwards
		    cetData);

	    numBM = cetData.size();

	    // Insert benchmark CDS cash flow dates and option exercise dates in tree
	    for (idxB=0; idxB<numBM; idxB++) {
		    // Insert cash flow dates
		    KVector(TDate) &cfDates = cetData[idxB].cfDates;
		    numCF = cfDates.size();
		    for (idxC=0; idxC<numCF; idxC++) {
			    cetTree.Insert(cfDates[idxC]);
		    }

            // insert protection leg date ZeroReset in cetData and in cetTree
            cetData[idxB].protZero = KZeroReset(cetTree.GetCurveName(KV_CREDIT_PROT_DEFRECOV), 
                cetData[idxB].exerDate, cetData[idxB].matDate);
            cetTree.Insert(cetData[idxB].protZero, true);

            // insert option exercise date and maturity
		    cetTree.Insert(cetData[idxB].exerDate);
            cetTree.Insert(cetData[idxB].matDate);

            #if defined(__EDEV__)
		            // Insert exp date discount factor for express dev
		            cetTree.InsertStatePrice(cetData[idxB].exerDate);
            #endif
        }
	    // Final tree set up
	    cetTree.SetUpTimeline();
	    cetTree.TreeSetUp();
	    cetTree.CheckTreeValid();
        /*=====================================================================
         * NOW FINISHED TREE AND BM SWAPTION SET-UP
         *===================================================================*/

        /*=====================================================================
         * GET VOL DATES AND SPOT VOLS FROM VNFM TREE
         *===================================================================*/
	    factVolDates = cetTree.SpotVolDates();
        factVolCurves = cetTree.SpotVolatilities();
        //el = fopen("C:\\ceterr.log","a");
        //fprintf(el,"Spot vols\n");
        //for (idx=0; idx<factVolDates.size(); idx++) {
        //    fprintf(el,"\t%ld\t\t%lf\t\t%lf\n",factVolDates[idx]-todayDate,
        //        (factVolCurves[0].size()>idx ? factVolCurves[0][idx] : 0.0),
        //        (factVolCurves[1].size()>idx ? factVolCurves[1][idx] : 0.0));
        //    fflush(el);
        //}
        //fclose(el);

        /*=====================================================================
         * IR VNFM ON CET-CALIBRATED irVolDiag SO THIS INCLUDES IR CET
         *===================================================================*/
        //DppCalibIRVolDiag(
        //    irMrParam,
        //    irSmileParam,
        //    irDim,
        //    irVolDiag,
        //    marketCurves.GetIRDiffuse(),
        //    factVolDates,
        //    factVolCurves);

        /*=====================================================================
         * CR VFNM APPROXIMATION
         *===================================================================*/
        //DppCalibCRSpreadVolDiag(
        //    crMrParam,
        //    irMrParam,
        //    irCrCorr,
        //    crSmileParam,
        //    marketCurves.GetDiffuse().BaseDate(),
        //    marketCurves.mCV[marketCurves.GetCurveType(marketCurves.mIRDiscCVName)],
        //    marketCurves.mCV[KV_CREDIT_RISKY],
        //    crVolDiag,
        //    marketCurves.mRecovery,
        //    factVolDates,
        //    factVolCurves);


        /*=====================================================================
         * SET UP PARTITION ON SPOT VOLS TO MATCH NUMBER OF BM INSTRUMENTS
         * FOR ITERATIVE SOLUTION. ALSO INITIALIZE INVERSE JACOBIAN
         *===================================================================*/
        KVector(int) volDatePartition;
        KVector(double) origVNFMVols(factVolCurves[1]); // keep copy
        KVector(double) fLast(numBM); // last tree price diffs from target
        KVector(double) xLast(numBM); // spot vol diffs from originals
        KVector(double) minVol(numBM);
        vegaDifferences.resize(numBM);
        KVector(KVector(double)) invJLast(numBM); // inverse Jacobian
        int j;
        idx = 0;
        for (idxB=0; idxB<numBM; idxB++) {
            /* INITIALIZE INVERSE JACOBIAN TO DIAGONAL OF 1/BS VEGA OF THE BM INSTRUMENTS */
            invJLast[idxB].assign(numBM, 0.0);
            invJLast[idxB][idxB] = (cetData[idxB].bsVega!=0 ? 1.0/cetData[idxB].bsVega : 1.0);
            minVol[idxB] = 1000.0;
            while (idx<factVolDates.size() && factVolDates[idx]<=cetData[idxB].exerDate) {
                minVol[idxB] = min(minVol[idxB], factVolCurves[1][idx]);
                idx++;
            }
            // add the index of the last voldate just before the benchmark expiry date.
            if (idxB>0 && idx==volDatePartition.back()) idx++; // no duplicates
            if (idx>=factVolDates.size()) {
                throw KFailure("%s: Vol date/benchmark mismatch.\n", routine);
            }
            volDatePartition.push_back(idx);
            xLast.push_back(0.0);
        }
        /* MAKE THE LAST ONE THE END, NO MATTER WHAT */
        volDatePartition[volDatePartition.size()-1] = factVolDates.size();
       
        /*=====================================================================
         * COMPUTE INITIAL TREE PRICE DIFFERENCES
         *===================================================================*/
        cetTree.UpdateFactorVol(factVolDates, factVolCurves);
        // Compute BM option prices
	    KCrxTreeCalcCet(
		    cetTree,
		    resetBank,
		    cetData);
        for (idxB=0; idxB<numBM; idxB++) {
            fLast[idxB] = cetData[idxB].treePrice - cetData[idxB].bsPrice;
            vegaDifferences[idxB].push_back(fLast[idxB]/cetData[idxB].bsVega);
        }

	    /*=====================================================================
         * DO THE ITERATIONS
         *===================================================================*/
	    for (iterCet=1; iterCet<=iterCetMax; iterCet++) {
            KVector(double) xDelta(numBM);
            KVector(double) fDelta(numBM);
            KVector(double) xThis(numBM);
            KVector(double) fThis(numBM);
            //fprintf(el,"\t\tCET ITERATION %d\n",iterCet);

            /*=================================================================
             * Calculate BM vector sized vol diffs: xDelta = -invJLast * fLast 
             * and xThis = xLast + xDelta
             *===============================================================*/
            for (i=0; i<numBM; i++) {
                xDelta[i] = 0.0;
                for (j=0; j<numBM; j++) {
                    xDelta[i] += -invJLast[i][j]*fLast[j];
                }
                xThis[i] = xLast[i] + xDelta[i];
                /* DON'T LET VOLS GO NEGATIVE */
                if (xThis[i]+minVol[i]<DBL_EPSILON) {
                    xThis[i] = 0.5*xLast[i] - 0.5*minVol[i];
                    xDelta[i] = xThis[i] - xLast[i];
                }

            }

            /*=================================================================
             * Set actual CR spot vol curve using vector vol differences
             *===============================================================*/
            for (idxB=0; idxB<numBM; idxB++) {
                for (idx=(idxB>0 ? volDatePartition[idxB-1] : 0); idx<volDatePartition[idxB]; idx++) {

                    factVolCurves[1][idx] = origVNFMVols[idx] + xThis[idxB];
                }
            }
            
            /*=================================================================
             * Recompute tree prices and differences
             *===============================================================*/
            cetTree.ClearTreeVolMem();
            cetTree.UpdateFactorVol(factVolDates,factVolCurves);
	        KCrxTreeCalcCet(
		        cetTree,
		        resetBank,
		        cetData);
            for (idxB=0; idxB<numBM; idxB++) {
                fThis[idxB] = cetData[idxB].treePrice - cetData[idxB].bsPrice;
                fDelta[idxB] = fThis[idxB]-fLast[idxB];
                vegaDifferences[idxB].push_back(fThis[idxB]/cetData[idxB].bsVega);
            }

            /*=================================================================
             * Recalculate inverse Jacobian estimate using Broyden's method
             * and Sherman-Morrison
             *===============================================================*/
            KVector(double) invJDeltaF(numBM);
            for (i=0; i<numBM; i++) {
                invJDeltaF[i] = 0.0;
                for (j=0; j<numBM; j++) {
                    invJDeltaF[i] += invJLast[i][j]*fDelta[j];
                }
            }
            KVector(double) deltaXInvJ(numBM);
            double p = 0;
            for (j=0; j<numBM; j++) {
                deltaXInvJ[j] = 0.0;
                for (i=0; i<numBM; i++) {
                    deltaXInvJ[j] += xDelta[i]*invJLast[i][j];
                }
                p += xDelta[j]*invJDeltaF[j];
            }
            for (i=0; i<numBM; i++) {
                for (j=0; j<numBM; j++) {
                    invJLast[i][j] = invJLast[i][j] + ((xDelta[i] - invJDeltaF[i])*deltaXInvJ[j])/p;
                }
              
                xLast[i] = xThis[i];
                fLast[i] = fThis[i];
            }
        } /* END OF ITERATION LOOP */

        

        
            
        /*=====================================================================
         * SET THE FACTOR VOLS FOR THE PRODUCTION TREE
         *===================================================================*/
        if (iterCet>1) {
            for (idxB=0; idxB<numBM; idxB++) {
                for (idx=(idxB>0 ? volDatePartition[idxB-1] : 0); idx<volDatePartition[idxB]; idx++) {
                    factVolCurves[1][idx] = origVNFMVols[idx] + xLast[idxB];
                }
            }
        }
        
        prodTree.InitializeFactorVol(factVolDates, factVolCurves);
    
        // log the convergence
        

	    // Printout (if required)
	    if (debugLevel > DEBUG_LEVEL_CET) {
            dppLog << endl << 
            format("CREDIT CET CONVERGENCE: CR mean-reversion=%lf CR leftQ=%lf, "
                "CR rightQ=%lf, IR/CR correlation=%lf\n"
                "BM\tBM Expiry\tBM Maturity\tBS Impl Vol\tBS Price\tPS "
                "Vega\tVega Differences (1.0=1%) at Iteration Number",
                crMrParam.mBeta[0], crSmileParam.mQ1, crSmileParam.mQ2, irCrCorr);
            for (idx=0; idx<vegaDifferences[0].size(); idx++) {
                dppLog << "\t" << idx;
            }
            dppLog << endl;
            for (idxB=0; idxB<numBM; idxB++) {
                dppLog << format("%d\t%s\t%s\t%lf\t%lf\t%lf",
                    idxB, GtoFormatDate(cetData[idxB].exerDate), GtoFormatDate(cetData[idxB].matDate),
                    cetData[idxB].bsVol, cetData[idxB].bsPrice, cetData[idxB].bsVega);
                for (idx=0; idx<vegaDifferences[idxB].size(); idx++) {
                    dppLog << "\t" << 100*vegaDifferences[idxB][idx];
                }
                dppLog << endl;
            }

            dppLog << endl << "SPOT VOLS:\nNumber\tDate\tRates\tCredit" << endl;
            for (idx=0; idx<factVolCurves[0].size(); idx++) {
                dppLog <<  format("%d\t%s\t%lf\t%lf", idx, GtoFormatDate(factVolDates[idx]), 
                    factVolCurves[0][idx],factVolCurves[1][idx]) << endl;
            }
            dppLog.flush();
	    }

    } catch (KFailure) {

        // Printout (if required)
	    if (debugLevel > DEBUG_LEVEL_CET) {
            dppLog << endl << 
            format("CREDIT CET CONVERGENCE: CR mean-reversion=%lf CR leftQ=%lf, "
                "CR rightQ=%lf, IR/CR correlation=%lf\n"
                "BM\tBM Expiry\tBM Maturity\tBS Impl Vol\tBS Price\tPS "
                "Vega\tVega Differences (1.0=1%) at Iteration Number",
                crMrParam.mBeta[0], crSmileParam.mQ1, crSmileParam.mQ2, irCrCorr);
            for (idx=0; idx<vegaDifferences[0].size(); idx++) {
                dppLog << "\t" << idx;
            }
            dppLog << endl;
            for (idxB=0; idxB<numBM; idxB++) {
                dppLog << format("%d\t%s\t%s\t%lf\t%lf\t%lf",
                    idxB, GtoFormatDate(cetData[idxB].exerDate), GtoFormatDate(cetData[idxB].matDate),
                    cetData[idxB].bsVol, cetData[idxB].bsPrice, cetData[idxB].bsVega);
                for (idx=0; idx<vegaDifferences[idxB].size(); idx++) {
                    dppLog << "\t" << 100*vegaDifferences[idxB][idx];
                }
                dppLog << endl;
            }

            dppLog << endl << "SPOT VOLS:\nNumber\tDate\tRates\tCredit" << endl;
            for (idx=0; idx<factVolCurves[0].size(); idx++) {
                dppLog <<  format("%d\t%s\t%lf\t%lf", idx, GtoFormatDate(factVolDates[idx]), 
                    factVolCurves[0][idx],factVolCurves[1][idx]) << endl;
            }
            dppLog.flush();
	    }

        throw KFailure("%s: failed.\n", routine);
    }
}



