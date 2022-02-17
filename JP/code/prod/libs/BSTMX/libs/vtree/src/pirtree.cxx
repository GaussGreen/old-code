/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	
 * Function:	
 * Author:	Th original model has been developped by 
 *	        Lionnel Pradier & London DR
 * 		This version has been adapted/modified NY DR.
 ***************************************************************/
#include "kstdinc.h"    /* Standard definitions & error hadling */
#include "kutilios.h"

extern  "C" {
#include "drlroot.h"
#include "drlsort.h"
}

#define	_kmrntree_SRC
#include "kpirtree.h"
#include "kmodpar.h"

#include "kvtspar.h"        // Slice parser
//
// Q cutoff
//
static	const	double	QCUTOFF = 1e-4;

#ifndef IS_Q
#define IS_Q(q) (fabs((q)) > QCUTOFF)
#endif

#ifndef SHIFT_ZERO
#define SHIFT_ZERO(x) ((x) = (IS_ALMOST_ZERO(x) ? 1e-5 : (x)))
#endif

extern  KMap(String,double)     globVpConstTable;        // Constant tables

//#define	__NO_CALIB__		// do not perform drift calibration

//--------------------------------------------------------------
//

KPirTree::KPirTree() : KMrNTree()
{
	mIRDim = 1;

	// Lognormal
	mQLo   = 1e0; 
	mQHi   = 1e0;
	mFSh   = 0e0;
	mTpZCenter = NULL;
}

//--------------------------------------------------------------
//

KPirTree::~KPirTree()
{
	for(KMap(TDate,double*)::iterator iterDev=mDevStatePr.begin();
            iterDev != mDevStatePr.end(); ++iterDev) {
                sliceDelete((*iterDev).second);
        }

	for(KVector(int)::iterator iterCV  = mCVTypes.begin();
				   iterCV != mCVTypes.end(); ++iterCV) 
	{
		sliceDelete(mDiscountIR[*iterCV]);
		delete [] mTpZPrices[*iterCV];
		delete [] mTpFRates[*iterCV];
	}

	if (mTpZCenter)
		delete [] mTpZCenter;

}




//--------------------------------------------------------------
//

void 
KPirTree::DeleteMemory()
{
	KMrNTree::DeleteMemory();

	for(KMap(TDate,double*)::iterator iterDev=mDevStatePr.begin();
            iterDev != mDevStatePr.end(); ++iterDev) {
                sliceDelete((*iterDev).second);
        }
	

	for(KVector(int)::iterator iterCV  = mCVTypes.begin();
				   iterCV != mCVTypes.end(); ++iterCV) 
	{
		sliceDelete(mDiscountIR[*iterCV]);
		delete [] mTpZPrices[*iterCV];
		delete [] mTpFRates[*iterCV];
	}

	delete [] mTpZCenter;

	mDiscountIR.clear();
	mTpZPrices.clear();
	mTpFRates.clear();

	mCVTypes.clear();
	mZBanks.clear();
	mKZCurves.clear();
	mDevDates.clear();
	mDevStatePr.clear();

	
	// Re-initialize the tree member
	mTpZCenter = NULL;

}





//--------------------------------------------------------------
// Free tree memories that depends on the factor vols, including
// tree limits, transitional probabilities. discount slicess, etc.
//
void
KPirTree::ClearTreeVolMem()
{
static  char	routine[] = "KPirTree::ClearTreeVolMem";

	KMrNTree::ClearTreeVolMem();

	//
	// Free discount slices
	//
	for(KVector(int)::iterator iterCV  = mCVTypes.begin();
				   iterCV != mCVTypes.end(); ++iterCV) 
		sliceDelete(mDiscountIR[*iterCV]);

	mDiscountIR.clear();
}



//--------------------------------------------------------------
// Insert crtical zero dates in the tree and zero bank
void
KPirTree::Insert(const KZeroReset &zeroReset, bool isCrit)
{
static	char	routine[] = "KPirTree::Insert(ZeroReset)";

   try {
	
	// Check if matures in the future 
	//
	if (zeroReset.mMaturityDate >= mTodayDate)
	{
		// Check valid dates
		zeroReset.IsValid();

		// Add to the tree if are critcal dates
		if (isCrit)
		{
			Insert(zeroReset.mMaturityDate);
			Insert(zeroReset.mEarliestDate);
		}

		// Insert in the corresponding zero bank.
		int curveIdx = GetCurveIdx(zeroReset.mCurveName);
		GetZeroBank(curveIdx).InsertDates(zeroReset.mMaturityDate, 
					  	  zeroReset.mEarliestDate);

	}

   }    // try block
   catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
   }
}




//--------------------------------------------------------------
// Insert rate cash flow dates in the tree and zero bank
TDate
KPirTree::Insert(const KRateReset &rateReset, bool isCrit)
{
static	char	routine[] = "KPirTree::Insert(RateReset)";

	TDate			resetDate = rateReset.ResetDate();
	KVector(KZeroReset)	zeroResetList;
	long            	curveIdx;
	
    try {
	

	// Reset date in the future
	if (resetDate >= mTodayDate)
	{
	    Insert(resetDate);	// Always add reset date

	    if (rateReset.Rate().IsFloating()) {

		curveIdx = GetCurveIdx(rateReset.Rate().CurveName());

		// TRUE floating rate
		//
		// zero reset dates start from effDate.
		rateReset.GetZeroResetList(zeroResetList);

                if(zeroResetList.size() < 2)
			throw KFailure( "%s: invalid input rateReset, "
				        "Number of zero resets is %d < 2.\n"
					"Minimum requires effective date and "
					"maturity date for a simple rate.\n",
					routine,
					zeroResetList.size());

		// Add to the tree and zero bank
		for (KVector(KZeroReset)::iterator 
				iterZero=zeroResetList.begin();
				iterZero!=zeroResetList.end();
				++iterZero)
		{

			// Add to the tree if is critcal date
			if (isCrit) 
				Insert(iterZero->mMaturityDate);

			// Insert in the corresponding zero bank.
			GetZeroBank(curveIdx).InsertDates(
						iterZero->mMaturityDate,
						resetDate);
		}

            }  // if floating

            return resetDate;

	}  // if reset in the future
	else
            return resetDate;

    }
    catch (KFailure) {
        throw KFailure("%s: failed inserting at %s.\n", 
			routine,
			GtoFormatDate(resetDate));
    }

}

//--------------------------------------------------------------
// Insert complex rate cash flow dates in the tree and zero bank
// Same reset date for all instruments.
TDate
KPirTree::Insert(const KCplxRateReset &cplxRateReset, bool isCrit)
{
static	char	routine[] = "KPirTree::Insert(CplxRateReset)";
    TDate       resetDate;
    int         i, nbInstr;
    try
    {
        nbInstr = cplxRateReset.CplxRate().nbInstr();
        for ( i = 0; i < nbInstr; i++)
            resetDate = Insert(cplxRateReset.RateResetIndex(i), isCrit);
        return resetDate;
    }
    catch (KFailure) {
        throw KFailure("%s: failed inserting at %s.\n", 
			routine,
			GtoFormatDate(resetDate));
    }
    	
}



//--------------------------------------------------------------
// Insert rate cash flow dates in the tree and zero bank
TDate
KPirTree::Insert(const KRateReset &rateReset, TDate endDate, bool isCrit)
{
static	char	routine[] = "KPirTree::Insert(RateReset, endDate)";

	TDate	resetDate = rateReset.ResetDate();
	TDate	lastZeroDate;
	
    try {
	

	// Reset date in the future
	if (resetDate >= mTodayDate)
	{
	    Insert(resetDate);	// Always add reset date

	    if (rateReset.Rate().IsFloating()) {

		Insert(rateReset, isCrit);

		//
		// Added extra zero matured at "endDate + rateReset.SpotOffset
		// + rateReset.Maturity()" to zero list all reset at 
		// original resetDate.  This represents the last possible
		// rate reset for american type of excersize,
		// and allow us to compute the rate at any
		// reset date between resetDate and endDate,
		// with interpolation if needed.
		//
		lastZeroDate = endDate 
			     + rateReset.Rate().SpotOffset()
			     + rateReset.Rate().Maturity();

		Insert( KZeroReset(rateReset.Rate().CurveName(),
				   resetDate,
				   lastZeroDate), 
			isCrit);
 
            }  // if floating

            return resetDate;

	}  // if reset in the future
	else
            return resetDate;
		

    }
    catch (KFailure) {
        throw KFailure("%s: failed inserting between %s and %s.\n",
			 routine,
			 GtoFormatDate(resetDate),
			 GtoFormatDate(endDate));
    }

}

//--------------------------------------------------------------
// Insert complex rate cash flow dates in the tree and zero bank
// Same reset date for all instruments.
TDate
KPirTree::Insert(
                 const KCplxRateReset &cplxRateReset, 
                 TDate endDate, 
                 bool isCrit)
{
static	char	routine[] = "KPirTree::Insert(CplxRateReset, endDate)";
    TDate       resetDate;
    int         i, nbInstr;
    try
    {
        nbInstr = cplxRateReset.CplxRate().nbInstr();
        for ( i = 0; i < nbInstr; i++)
           resetDate = Insert(cplxRateReset.RateResetIndex(i), 
                              endDate,
                              isCrit);
        
        return resetDate;
    }
    catch (KFailure) {
        throw KFailure("%s: failed inserting at %s.\n", 
			routine,
			GtoFormatDate(resetDate));
    }
    	
}

//---------------------------------------------------------------
// Compute a zero reset.  The pure interest tree can
// deal with two types of zeros: 
// 1. rate out of diffuse curve.
// 2. rate with deterministic spread over diffuse curve.
//
void    
KPirTree::Get(KTSlice &zeroTS, const KZeroReset &zeroReset)
{
static	char	routine[] = "KPirTree::Get(KZeroReset)";

	int	curveIdx;

	String	zeroTSCVName;


	TDate	resetDate = zeroReset.mEarliestDate;
	TDate	matDate   = zeroReset.mMaturityDate;

	int	currTP   = TPIdxCurrent();
	TDate	currDate = TPDateCurrent();

try {
	// Ensure target ts is allocated
        TSliceCreate(zeroTS);

	// Save the curve name if available
	zeroTSCVName = zeroTS.GetCurveName();

	curveIdx = GetCurveIdx(zeroReset.mCurveName);

	// Check that reset date >= current date
	//
	if (resetDate < currDate)	
		throw KFailure( "%s: reset date(%s) < current date(%s).\n",
				routine,
				GtoFormatDate(resetDate),
				GtoFormatDate(currDate));


	// Compute the index rate from zero bank and copy to zeroTS
	GetZeroBank(curveIdx).GetZeroSlice(zeroTS, matDate);

	// Set the time point of the slice
	zeroTS.SetTpIdx(currTP);

	// Restore the slice curve name if is NON-default
    // But keep the curveIdx for dimension changes.
    //
	if (zeroTSCVName != K_DEFAULT_NAME)
		zeroTS.SetSliceName(zeroTSCVName);

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }

}




//---------------------------------------------------------------
// Compute a reset index rate.  The pure interest tree can
// deal with two types of rate: 
// 1. rate out of diffuse curve.
// 2. rate with deterministic spread over diffuse curve.
//
void
KPirTree::Get(KTSlice &rateTS, const KRateReset &rateReset)
{
static	char	routine[] = "KPirTree::Get(KRateReset)";

	int	curveIdx;

	TDate	resetDate    = rateReset.ResetDate();
	TDate   resetEffDate = rateReset.EffDate();


	String	rateTSCVName;	// curve name for rateTS

	int	currTP   = TPIdxCurrent();
	TDate	currDate = TPDateCurrent();

	double	resetValue;

	double	resetFwd, currFwd;


try {

	// Ensure that slice is allocated.
	TSliceCreate(rateTS);

	// Save the curve name if available
	rateTSCVName = rateTS.GetCurveName();

	// Check if is fixed rate or floating rate
	//
	if (rateReset.Rate().IsFloating())
	{
	    curveIdx = GetCurveIdx(rateReset.Rate().CurveName());

	    // 
	    // Check if has already been manually reset
	    //
	    // We strip the spread out to check reset bank
	    // and will add it back
            KRateReset rateResetNS = rateReset;
	    rateResetNS.Rate().SetSpread(0e0);

	    if(mResetBank->Get(rateResetNS, &resetValue))
	    {
		rateTS = resetValue + rateReset.Rate().Spread();
	    }
	    else	// no reset value
	    {
		if (resetDate < mTodayDate) {
			throw KFailure("%s: (today %s) no available reset "
				"value on %s.\n", routine,
				GtoFormatDate(mTodayDate),
				GtoFormatDate(resetDate));

		}

		// Temp slice for annuity
		KTSlice annuity(*this, 
				format("Annuity reset on %s", 
					    GtoFormatDate(resetDate)),
				rateReset.Rate().CurveName());


		//
	    	// Reset strictly in the past, need to adjust by the forward
		// ratio.
		// If resetDate < currDate <= resetEffDate, 
		// use (currDate, resetEFfDate) as reset and reset effective
		// dates.  This would use correct zeros to produce
		// correct rate, rather than approximation.  
		// 
		if (currDate <= resetDate) 
		    GetZeroBank(curveIdx).ComputeParYieldAndAnnuity(currTP,
						        	rateReset,
						        	rateTS,  
						        	annuity);
		else if (currDate <= resetEffDate)
		{
		    // Same KRate reset at (currDate, resetEffDate)
		    //
		    KRateReset currRateReset(rateReset.Rate().CurveName(),
					     currDate,
					     resetEffDate, 
					     rateReset.Rate());

		    GetZeroBank(curveIdx).ComputeParYieldAndAnnuity(currTP,
						        	currRateReset,
						        	rateTS,  
						        	annuity);
                   
		}
	    	else 	// (resetEffDate < currDate). truely reset in the past
		{
		    // Same KRate reset at current date
		    //
		    KRateReset approxRateReset(
					rateReset.Rate().CurveName(),
					currDate,
					currDate+rateReset.Rate().SpotOffset(), 
					rateReset.Rate());

		    GetZeroBank(curveIdx).ComputeParYieldAndAnnuity(currTP,
						        	approxRateReset,
						        	rateTS,  
						        	annuity);

		    resetFwd = rateReset.Rate().Forward(GetKZCurve(curveIdx),
						 resetDate);

		    currFwd  = rateReset.Rate().Forward(GetKZCurve(curveIdx),
						 currDate);
 
		    rateTS *= (resetFwd / currFwd);
		}
	    }
	}
	else	// fixed rate
	    rateTS = rateReset.Rate().Spread();


	// Set the time point of the slice
	rateTS.SetTpIdx(currTP);

	// Restore the slice curve name if is NON-default
	if (rateTSCVName != K_DEFAULT_NAME)
		rateTS.SetCurveIdx(rateTSCVName);

    dppLog << routine << '\t' << GtoFormatDate(TPDateCurrent()) << '\t' << currTP<< endl;
    this->slicePrint(
                    (double*) rateTS.mData,
                    rateTS.GetSliceDim(),
                    currTP,
                    FALSE,
                    dppLog);
                    
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }

}

//---------------------------------------------------------------
// Compute a reset index rate.  The pure interest tree can
// deal with complex rate via payment formula
// 
void
KPirTree::Get(KTSlice &rateTS, const KCplxRateReset &cplxRateReset)
{
static	char	routine[] = "KPirTree::Get(KCplxRateReset)";
    KTSlice     *rateArrayTs = NULL;
    KRateReset  *rateReset = NULL;
    int         i, nbInstr;
    TDate       resetDate;
    String      rateTSCVName;	// curve name for rateTS
    int	        currTP   = TPIdxCurrent();
    try
    {
      // Ensure that slice is allocated.
        TSliceCreate(rateTS);

    	// Save the curve name if available
	    rateTSCVName = rateTS.GetCurveName();
        
        nbInstr = cplxRateReset.CplxRate().nbInstr();
        resetDate = cplxRateReset.ResetDate();

        rateArrayTs = KTSliceNewVector(
                        nbInstr,
                        *this,
                        cplxRateReset.CplxRate().CurveName(), 
                        format("rate reset slices at %s",
                            DrlTDatePrint(NULL, resetDate)));
        for (i = 0; i < nbInstr; i++)
        {
            if (cplxRateReset.RateIndex(i).IsFloating())
                this->Get(rateArrayTs[i], 
                      cplxRateReset.RateResetIndex(i));
            else
                // fixed rate
                rateArrayTs[i] = cplxRateReset.RateIndex(i).Spread();
                            
        }
        KTSliceParEval(
                        *this,
                        nbInstr,
                        rateArrayTs,
                        cplxRateReset.CplxRate().getFormula().c_str(),
                        globVpConstTable,        // Constant tables
                        &rateTS);
	
        // Set the time point of the slice
	    rateTS.SetTpIdx(currTP);

	    // Restore the slice curve name if is NON-default
	    if (rateTSCVName != K_DEFAULT_NAME)
		    rateTS.SetCurveIdx(rateTSCVName);

/*        dppLog << "currentDate " << GtoFormatDate(TPDateCurrent()) << endl;
        dppLog << cplxRateReset.RateResetIndex(0) << endl;
        dppLog << rateArrayTs[0];*/

        delete[] rateArrayTs;

    }
    catch (KFailure) 
    {
        delete[] rateArrayTs;
	    throw KFailure("%s: failed.\n", routine);
    }
}
//--------------------------------------------------------------
//

TDate
KPirTree::GetValueDate(int curveIdx)
{
static	char	routine[] = "KPirTree::GetValueDate";

	KMap(int, TDate)::iterator itDate = mValueDates.find(curveIdx);

	if (itDate == mValueDates.end())
		throw KFailure("%s: invalid curve index (%d). "
			       "Value date for curve [%d] does not exist.\n",
				routine, curveIdx, curveIdx);

	return (*itDate).second;

}



//--------------------------------------------------------------
//

KZCurve&
KPirTree::GetKZCurve(int curveIdx)
{
static	char	routine[] = "KPirTree::GetKZCurve";

	KMap(int, KZCurve)::iterator itCV = mKZCurves.find(curveIdx);

	if (itCV == mKZCurves.end())
		throw KFailure("%s: invalid curve index (%d). "
			       "Zero curve [%d] does not exist.\n",
				routine, curveIdx, curveIdx);

	return (*itCV).second;

}


//--------------------------------------------------------------
//

KZeroBank&
KPirTree::GetZeroBank(int curveIdx)
{
static	char	routine[] = "KPirTree::GetZeroBank";

	KMap(int, KZeroBank)::iterator itZB = mZBanks.find(curveIdx);

	if (itZB == mZBanks.end())
		throw KFailure("%s: invalid curve index (%d). "
			       "Zero bank [%d] does not exist.\n",
				routine, curveIdx, curveIdx);

	return (*itZB).second;

}



//--------------------------------------------------------------
//

double*
KPirTree::GetDiscount(int curveIdx)
{
static	char	routine[] = "KPirTree::GetDiscount";

	KMap(int, double*)::iterator it = mDiscountIR.find(curveIdx);

	if (it == mDiscountIR.end())
		throw KFailure("%s: invalid curve index (%d). Discount "
			       "factor for curve [%d] does not exist.\n",
				routine, curveIdx, curveIdx);

	return (*it).second;

}





//--------------------------------------------------------------
//

double*
KPirTree::GetZeroPrice(int curveIdx)
{
static	char	routine[] = "KPirTree::GetZeroPrice";

	KMap(int, double*)::iterator it = mTpZPrices.find(curveIdx);

	if (it == mTpZPrices.end())
		throw KFailure("%s: invalid curve index (%d). "
			       "Zero price for curve [%d] does not exist.\n",
				routine, curveIdx, curveIdx);

	return (*it).second;

}



//--------------------------------------------------------------
//

double
KPirTree::GetZeroPrice(int curveIdx, TDate dt)
{
static	char	routine[] = "KPirTree::GetZeroPrice(TDate)";

        int     tpL;
        double  result;
        double  dccFactor, offsetTime, zeroRate;
        double  *zeroPrice = NULL;
        TCurve  *tCurve = NULL;

	KMap(int, double*)::iterator it = mTpZPrices.find(curveIdx);

	if (it == mTpZPrices.end())
		throw KFailure("%s: invalid curve index (%d). "
			       "Zero price for curve [%d] does not exist.\n",
				routine, curveIdx, curveIdx);

	zeroPrice = (*it).second;

        IF_FAILED_THROW (DrlTDateArrayFloorIdx(TPDates, TPNum(),
                                               dt,
                                               &tpL));  

        if (TPDates[tpL] == dt)
            result = zeroPrice[tpL];
        else
        {
            tCurve = (TCurve*)GetKZCurve(curveIdx);

            // Check zero curve day count is OK
            switch (tCurve->fDayCountConv)
            {
            case GTO_ACT_365F:
                dccFactor = 1e0;
                break;
            case GTO_ACT_360:
                dccFactor = 365e0/360e0;
                break;
            default:
                throw KFailure("%s: Day count conv for zero curve %d must be"
                        "\n either Act/360 or Act/365F, not %s.\n",
                        routine, curveIdx,
                        GtoFormatDayCountConv(tCurve->fDayCountConv));
            }

            offsetTime = ((double)dt - (double)TPToday()) / 365e0;

            // Interpolate zero rate on input date
            IF_FAILED_THROW( GtoInterpRateExact(
                        dt,
                        tCurve,
                        GetKZCurve(curveIdx).ZeroInterpType(),
                        &zeroRate));

            // Convert the rate to a discount factor.
            IF_FAILED_THROW( GtoRateToDiscountYearFrac(
                        zeroRate,
                        offsetTime * dccFactor,
                        tCurve->fBasis,
                        &result));
        }

        return result;
}




//--------------------------------------------------------------
//

double*
KPirTree::GetForwardRate(int curveIdx)
{
static	char	routine[] = "KPirTree::GetForwardRate";

	KMap(int, double*)::iterator it = mTpFRates.find(curveIdx);

	if (it == mTpFRates.end())
		throw KFailure("%s: invalid curve index (%d). "
			       "Forward rate for curve [%d] does not exist.\n",
				routine, curveIdx, curveIdx);

	return (*it).second;

}




//--------------------------------------------------------------
//

double*
KPirTree::GetStatePrice(TDate devDate)
{
static	char	routine[] = "KPirTree::GetStatePrice";

	KMap(TDate, double*)::iterator it = mDevStatePr.find(devDate);

	if (it == mDevStatePr.end())
		throw KFailure("%s: invalid dev date (%s). "
			       "State price on %s does not exist.\n",
				routine, 
				GtoFormatDate(devDate),
				GtoFormatDate(devDate));

	return (*it).second;

}



//--------------------------------------------------------------
//
double          
KPirTree::GetVolBbq(int t)
{
static	char	routine[] = "KPirTree::GetVolBbq";

	double  QLeft    = this->mQLo;
	double  QRight   = this->mQHi;
	double  FwdShift = this->mFSh;

	double	QMid;

	double  FwdRate;	// Fwd rate 
	double  InsRate;	// Instantaneous rate  

	double	VolBbq;


	//
	// Calc needed constants
	//
	FwdRate = GetForwardRate(KV_DIFF)[t];

	// Avoid zero forward rate, which is allowed in normal case but
        // canceled out after the mapping.  Set the minimum as 0.1bp
	SHIFT_ZERO(FwdRate);

	QMid = (QLeft + QRight) / 2;
  

        InsRate   = FwdRate / Length[t];    

        VolBbq  = mVolLogn * mBbq  
                + mVolNorm * (1. - mBbq) / InsRate;

        VolBbq *= (1. + FwdShift) / (1. + QMid * FwdShift);


	if (debugLevel >= DEBUG_LEVEL_DRIFT) {
		dppLog << format("%s: TPIDX=%4d VolBbq=%14.10f\n",
			routine, t, VolBbq);
	}

	return VolBbq;
}


//--------------------------------------------------------------
// Create a slice of dimension mIRDim.  

KTSlice&
KPirTree::TSliceCreate(KTSlice &ts)
{
static	char	routine[] = "KPirTree::TSliceCreate";

try{

	return TSliceDimCreate(ts, mIRDim);

   } 
   catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
   }

}




//--------------------------------------------------------------
// Slice Dev.  Should be called after tree being updated, i.e.
// the discount factor and probabality are calculated. 

KTSlice&
KPirTree::TSliceDev(KTSlice &ts, const String &discCurveName)
{
static	char	routine[] = "KPirTree::TSliceDev";
	
	int	discCurveIdx;
	
try {
	// Dev empty slices does nothing
	if (ts.mData == NULL) return ts;

	// Don't do anything at last timept
	if (this->tpIdxCurrent == TPNum()) return ts;

	// Check timepoint consistency
	if (ts.GetTpIdx() != tpIdxCurrent+1)
		throw KFailure("%s: attempting to Ev slice `%s' "
			"with TP %d on tree at TP %d.\n",
			routine, ts.GetSliceName().c_str(), 
			ts.GetTpIdx(), tpIdxCurrent);

	// Check slice dimension
	if (ts.GetSliceDim() != mIRDim)
		throw KFailure("%s: slice dimension (%d) does not equal to "
			       "the IR tree dimension (%d).\n",
			       routine, ts.GetSliceDim(), mIRDim);
	
	// Discount curve index 
	//
	discCurveIdx = this->GetCurveIdx(discCurveName);

	// Perform EV
	KMrNTree::sliceEv( 
			  (double*) ts.mData,
			  (int) ts.mSliceDim,
			  GetDiscount(discCurveIdx),
			  tpIdxCurrent);

	// Set new TP
	ts.SetTpIdx(tpIdxCurrent);

	return (ts);

   }    // try block
   catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
   }

}



//--------------------------------------------------------------
// Slice Forward.  Should be called after tree being updated, i.e.
// the discount factor and probabality are calculated. 

KTSlice&
KPirTree::TSliceFwd(KTSlice &ts)
{
static	char	routine[] = "KPirTree::TSliceFwd";
	
	int	curveIdx;
	int	sliceDim;
	
try {
	// Fwd empty slices does nothing
	if (ts.mData == NULL) return ts;

	// Don't do anything at last timept
	if (this->tpIdxCurrent == TPNum()) return ts;

	// Check timepoint consistency
	if (ts.GetTpIdx() != tpIdxCurrent)
		throw KFailure("%s: attempting to forward slice `%s' "
			"with TP %d on tree at TP %d.\n",
			routine, ts.GetSliceName().c_str(), 
			ts.GetTpIdx(), tpIdxCurrent);

	curveIdx = ts.GetCurveIdx();
	sliceDim = ts.GetSliceDim();

	// Perform FWD based on slice dimension
	if (sliceDim == mIRDim)
		sliceFw( 
		  	(double*) ts.mData,
			sliceDim,
			GetDiscount(curveIdx),
			tpIdxCurrent);
	else if (sliceDim == mNbFactor)
	{
		//
		// Expand DiscountIR in nDim dimensions, 
		// which only varies in the first nIRDim 
		// dimensions and is constant in the 
		// higher nDim-nIRDim dimensions
		//
		double *discount = sliceNew(mNbFactor);

		sliceExpand(tpIdxCurrent, mNbFactor, mIRDim, 
			    GetDiscount(curveIdx),
			    discount);
		
		sliceFw( 
		  	(double*) ts.mData,
			sliceDim,
			discount,
			tpIdxCurrent);

		sliceDelete(discount);
	}
	else
		throw KFailure("%s: invalid slice dimension (%d) "
			       "the IR tree dimension (%d).\n",
			       routine, ts.GetSliceDim(), mIRDim);
	

	// Set new TP
	ts.SetTpIdx(tpIdxCurrent+1);

	return (ts);

   }    // try block
   catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
   }

}



//--------------------------------------------------------------
// Update the tree 

void
KPirTree::Update(int tpIdx)
{
static	char	routine[] = "KPirTree::Update";

	int	t = tpIdx;			// Time step index

	double	Zt;


    try {

	// Don't do anything at last timept
        if (tpIdx != TPNum()) 
	{
	    // Compute the transition probability between t and t+1
	    KMrNTree::Update(t);

	    // Calibrated center offset 
	    Zt = mTpZCenter[t];

	    // Calculate the discount between t and t+1
	    //
	    CalcDiscount(t, mIRDim, Zt, mDiscountIR);
	}


	// Update the zero bank
	for(KMap(int, KZeroBank)::iterator itZB  = mZBanks.begin();
					   itZB != mZBanks.end(); ++itZB)
		(*itZB).second.Update(t);

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}





//---------------------------------------------------------------
// Check the tree.
void    
KPirTree::CheckTreeValid()
{
static  char    routine[] = "KPirTree::CheckTreeValid";

	// Check if KV_DIFF curve exists
	if(find(mCVTypes.begin(), mCVTypes.end(), KV_DIFF) == mCVTypes.end())	
		throw KFailure("%s: No diffuse curve specified in the tree.\n",
				routine);


	if (debugLevel > DEBUG_LEVEL_TIMELINE) {
		dppLog << format("%s:\n", routine);
		dppLog << format("mIRDim = %d\n", 	mIRDim);
		dppLog << format("mQLo = %lf\n", 	mQLo);
		dppLog << format("mQHi = %lf\n", 	mQHi);
		dppLog << format("mFSh = %lf\n", 	mFSh);
		dppLog << format("mBbq = %lf\n", 	mBbq);
		dppLog << format("mVolNorm = %lf\n", 	mVolNorm);
		dppLog << format("mVolLogn = %lf\n", 	mVolLogn);
	}
}






//---------------------------------------------------------------
// Initialize the tree.

void    
KPirTree::Initialize(
	TDate		   todayDt,		// (I) today's date
	KMrParam           &irMrParam,		// (I) full dimension mr info
	KSmileParam        &irSmilePar,		// (I) ir smile info

	int		   nIRDim,		// (I) IR dimension
	
	KVector(int)	   &cvTypes,		// (I) curve types (KV_DIFF,.)
	KMap(int, KZCurve) &cv,    		// (I) array of curves
	KMap(int, String)  &cvNames,		// (I) array of curve names
	KMap(int, TDate)   &cvValueDates,	// (I) array of cv value dates
 
	KVector(TDate) 	   &volDates,   	// (I) volatility dates
	KVector(KVector(double)) &factVol,	// (I) spot vol
	KResetBank	   &resetBank)		// (I) rate reset bank
{
static	char	routine[] = "KPirTree::Initialize";

	char	EoI = 'I';
	double	smoothFact;
	
	int	nDim;

 try{
	if (!irMrParam.IsValid() ||
	    !irSmilePar.IsValid())
		throw KFailure("%s: Invalid model inputs.\n", routine);

	// Check dimensions
	nDim = irMrParam.mNumFact;
	if (nIRDim > nDim ||
	    nDim < 1      || nDim > 3   ||
	    nIRDim < 1    || nIRDim > 3 )
		throw KFailure("%s: invalid factor dimensions. "
			       "nDim = %d; nIRDim = %d.\n",
				routine, nDim, nIRDim);

	// Use smoothing factor from parameters
	smoothFact = irMrParam.mSmoothFact;

	//
	// Initialize the model parameters 
	//
	KMrNTree::Initialize(
			todayDt,
			irMrParam.mPpy,             
			smoothFact,
			EoI,            
			irMrParam.mNumStdevCut,      
			irMrParam.mNumFact,      
			irMrParam.mBeta,      
			irMrParam.mAlpha, 
			irMrParam.mRho,
			volDates,      
			factVol);


	// Smile parameters
	mIRDim = nIRDim;
	mQLo = irSmilePar.mQ1; 
	mQHi = irSmilePar.mQ2; 
	mFSh = irSmilePar.mQF;


	//
	// Backbone coefficients
	//
	mBbq = irMrParam.mBackboneQ;

	// Compute reference vol constants 
	double norm = 0e0;
	for (int idxF = 0; idxF<nIRDim; idxF++)
		norm += irMrParam.mAlpha[idxF] * irMrParam.mAlpha[idxF];
	norm = sqrt(norm);
	if (fabs(norm) < DBL_EPSILON)
		throw KFailure("%s: total alpha too small %g.\n",
			routine, norm);

	if (IS_ALMOST_ZERO(mBbq - 1e0)) {
		mVolNorm = 0.;
		mVolLogn = norm;
	} else if (IS_ALMOST_ZERO(mBbq - 0e0)) {
		mVolNorm = norm;
		mVolLogn = 0.;
	} else {
		throw KFailure("%s: backbone must be 0 or 1 (%g).\n",
			routine, mBbq);
	}
 

	//
	// Initialize market data
	//
	InitializeCurves(cvTypes,
			 cv,
			 cvNames,
			 cvValueDates);

	// Rate reset bank
	mResetBank = &resetBank;

   }    // try block
   catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
   }
}




//--------------------------------------------------------------
//
void
KPirTree::InitializeCurves(
	KVector(int)	   &cvTypes,	// (I) array of curve types(KV_DIFF...)
	KMap(int, KZCurve) &cv,    	// (I) array of curves
	KMap(int, String)  &cvNames,	// (I) array of curve names
	KMap(int, TDate)   &cvValueDates) // (I) array of cv value dates
{
static	char	routine[] = "KPirTree::InitializeCurves";

	int	curveIdx;

 try{
	if (!mZBanks.empty() ||
	    !mKZCurves.empty()){
		throw KFailure("%s: zero curves or banks are not empty.\n",
			       routine);
	}
	else{	
		this->mCVTypes  = cvTypes;
		this->mKZCurves = cv;
		this->mValueDates = cvValueDates;

		for(KMap(int,TDate)::iterator itDate = mValueDates.begin();
		    itDate != mValueDates.end(); ++itDate) 
			Insert((*itDate).second);

		for(KVector(int)::iterator iterCV = cvTypes.begin(); 
					   iterCV != cvTypes.end(); ++iterCV) {
			curveIdx = *iterCV;

			KMap(int, String)::iterator itStr 
						= cvNames.find(curveIdx);
			KMap(int, KZCurve)::iterator itZCV 
						= cv.find(curveIdx);
			if(itStr  != cvNames.end() && 
			   itZCV  != cv.end()) {
				MapCurveName((*itStr).second, curveIdx);

				mZBanks.insert(KMap(int, KZeroBank)::value_type(
				    curveIdx,
				    KZeroBank(*this,
					      format("Zero Bank %d", 
						     curveIdx),
					      (*itStr).second,
					      (*itZCV).second.ZeroInterpType())
				    ));
			}
			else
				throw KFailure("%s: incorrect mapping of "
					       "cvNames[%d], or cv[%d].\n",
						routine, curveIdx, curveIdx);

		}
	}	

   }    // try block
   catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
   }
}




//--------------------------------------------------------------
//  Compute zero prices and forward rates at each time step
void
KPirTree::ComputeZeroAndForward()
{
static	char	routine[] = "KPirTree::ComputeZeroAndForward";

	TCurve  *tCurve = NULL;

	double	todayDT,		// in days
		offsetDT,
		exactDT;
	double	offsetTime,		// in yrs
		dccFactor,
		zeroRate,
		zeroPrice;

	int	idxT;		 // Time step index

	int	curveIdx;

	int	nbZCurves = mCVTypes.size();

try{

    // Allocate arrays for zero curves
    //
    ASSERT_OR_THROW((this->mTpZCenter = new double [NbTP+2]) != NULL);


    for(KVector(int)::iterator iterCV = mCVTypes.begin();
	 		       iterCV != mCVTypes.end(); ++iterCV) 
    {
	curveIdx = *iterCV;

	mTpZPrices.insert(KMap(int, double*)::value_type(
				curveIdx,
				new double [NbTP+2]));

	mTpFRates.insert(KMap(int, double*)::value_type(
				curveIdx,
				new double [NbTP+2]));

	// Compute zero prices: if non integer timepoint, use
	// fractional interpolation
	//

	tCurve = (TCurve*)GetKZCurve(curveIdx);

	// Check zero curve day count is OK
        switch (tCurve->fDayCountConv)
	{
	case GTO_ACT_365F:
		dccFactor = 1e0;
		break;
	case GTO_ACT_360:
		dccFactor = 365e0/360e0;
		break;
	default:
		throw KFailure("%s: Day count conv for zero curve %d must be"
			"\n either Act/360 or Act/365F, not %s.\n", 
			routine, curveIdx,
			GtoFormatDayCountConv(tCurve->fDayCountConv));
	}

	todayDT = (double) this->mTodayDate;

	// Interpolate zero rates
	for (idxT=0; idxT<=NbTP; idxT++) {

	    // Compute fractional day
	    offsetTime = TPTimes[idxT];
	    offsetDT = offsetTime * 365e0;
	    exactDT  = todayDT + offsetDT;

	    if (TPDates[idxT] > 0L) {
		// Exact timepoint
		// Check we didn't f.. up
		double	res = exactDT - (double) TPDates[idxT];
		if (!IS_ALMOST_ZERO(res))
		    throw KFailure("%s: date %d (%s, %ld) incompatible with"
			" time %12.8f\n", routine,
			idxT, printDate(TPDates[idxT]), TPDates[idxT],
			exactDT);

	    }

	    if (IS_ALMOST_ZERO(offsetDT)) {
		// case of zero time
		zeroPrice = 1e0;
	    } else {
		// Interpolate zero rate on exact time
		IF_FAILED_THROW( GtoInterpRateExact(
			exactDT,
			tCurve,
			GetKZCurve(curveIdx).ZeroInterpType(),
			&zeroRate));

		// Convert the rate to a discount factor.
		IF_FAILED_THROW( GtoRateToDiscountYearFrac(
			zeroRate,
			offsetTime * dccFactor,
			tCurve->fBasis,
			&zeroPrice));
	    }

	    GetZeroPrice(curveIdx)[idxT] = zeroPrice;
	}


	// Compute forward rate 1-period rate (not annualized)
	for (idxT=0; idxT<NbTP; idxT++) {
		GetForwardRate(curveIdx)[idxT] = 
	    	   GetZeroPrice(curveIdx)[idxT] 
		   / GetZeroPrice(curveIdx)[idxT+1] - 1e0;
	}


    }   // *iterCV

   }    // try block
   catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
   }

}





 
//--------------------------------------------------------------
// Setup tree timeline.
// 1. Run the zero bank date optimation and insert the "optimal"
//    dates in the critical date list.
// 2. Set up tree timeline after all product dates are inserted.
//
void
KPirTree::SetUpTimeline()
{
static  char    routine[] = "KPirTree::SetUpTimeline";
 
  try {

	//
	// 1. Given all the zero maturity and usage dates,
	// find an "optimal" set of zero-bank maturity dates
	// and add them to the critcal date list.
	//
	for(KVector(int)::iterator iterCV = mCVTypes.begin();
                                          iterCV != mCVTypes.end(); ++iterCV)
		GetZeroBank(*iterCV).InitializeZeroDates();
	
	//
	// 2. Set up the tree timeline
	//
	KMrNTree::SetUpTimeline();

   }  
   catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
   }
}





//--------------------------------------------------------------
// Calibrate performs the tree set up and drift calibration.
//
void
KPirTree::Calibrate()
{
static	char	routine[] = "KPirTree::Calibrate";

try{
	// Set up the tree after all dates are inserted.
	TreeSetUp();

	// Check validity of tree parameters
	CheckTreeValid();
	
	// Calibrate center offset
	CalibrateDrift();

   }  
   catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
   }
}





//--------------------------------------------------------------
// TreeSetUp is called after tree timeline is computeds,
// and does the following:
// 1. Call KMrNTree::TreeSetUp to set up the tree parameters that only 
//    depend on vols(jump size, orthogonal factors, and tree limits, etc.)
// 2. Initialize temp discount slices for each curve
// 3. Compute the zero prices and forward rates at each time step.
// 4. Sort and merge mDevDates

void
KPirTree::TreeSetUp()
{
static	char	routine[] = "KPirTree::TreeSetUp";

try{
	// 1. Set up the tree 
	KMrNTree::TreeSetUp();

	// 2. Initialize temp discount slices for each curve
	//    after tree limit is set up.
	for(KVector(int)::iterator iter  = mCVTypes.begin();
                                   iter != mCVTypes.end(); ++iter) 
	{
		mDiscountIR.insert(KMap(int, double*)::value_type(
					*iter,
					sliceNew(this->mIRDim)));
	}

	// 3. Compute the zero prices and compute the forward rates
	// at each time step
	ComputeZeroAndForward();	


	// 4. sort and merge mDevDates in ascending order
	sort(mDevDates.begin(), mDevDates.end());
	KVector(TDate)::iterator iterBadDateSt
			= unique(mDevDates.begin(), mDevDates.end());
	mDevDates.erase(iterBadDateSt, mDevDates.end());

	if (debugLevel > DEBUG_LEVEL_TIMELINE) {
		dppLog << format("%s: Express dev dates (state prices)\n",
			routine) << endl;
		for (int idx=0; idx<mDevDates.size(); idx++) {
			dppLog << format(" %3d/%3d %10s", idx, mDevDates.size(),
				DrlTDatePrint(NULL, mDevDates[idx])) << endl;
		}
	}


   }    // try block
   catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
   }
}




//--------------------------------------------------------------
// Update all tree parameters that depend on factor vols, such as
// tree limit, discount slice, etc.  Used by CET only.
//
void
KPirTree::UpdateFactorVol(
        KVector(TDate) &volDates,       // (I) volatility dates
        KVector(KVector(double)) &factVol)// (I) spot vols
 
{
 
	KMrNTree::UpdateFactorVol(volDates,
				  factVol);

	//
	// Initialize temp discount slices for each curve
	// after tree limit is set up.
	//
	for(KVector(int)::iterator iter  = mCVTypes.begin();
                                   iter != mCVTypes.end(); ++iter) 
	{
		mDiscountIR.insert(KMap(int, double*)::value_type(
					*iter,
					sliceNew(this->mIRDim)));
	}
 
}



//--------------------------------------------------------------
// Drift calibration.  Calibrate the drift in mIRDim IR tree,
// and save the state price in nDim tree slice.

void
KPirTree::CalibrateDrift()
{
static	char	routine[] = "KPirTree::CalibrateDrift";

	int	t;			// Time step index
	
	int	nDim, nIRDim;

	double	FwdShift;
	double	*StatePr  = NULL,
		*StatePrIR = NULL,
		*DiscountIR = NULL,
		*Discount = NULL;

	double 	*StatePrL = NULL,
		*StatePrIRL = NULL,
		*devStatePr = NULL,
		*devStatePrL = NULL;

	double	zeroPrice,
		p;

	double	Zt,
		DelZt;

	KVector(TDate)::iterator iterDev;

    try {
	
	nDim   = mNbFactor;	// Total dimension
	nIRDim = mIRDim;	// IR dimension

	// Check dimension valid
        if (nDim   < 1 || nDim   > 3 ||
            nIRDim < 1 || nIRDim > 3)
                throw KFailure("%s: invalid factor numbers.\n"
                "Total dimension is %d, IR dimension is %d.\n",
                routine, nDim, nIRDim);

	FwdShift = mFSh;

	//
	// Allocate tmp working space
	//
	StatePrIR  = sliceNew(nIRDim);
	DiscountIR = GetDiscount(KV_DIFF);

	//
	// Set initial state price to 1 at TPIdx = 0.
	//
	StatePrIRL = StatePrIR + NodeOffset(nIRDim, 0, 0, 0);
    	StatePrIRL[0] = 1e0;

	//
	// Set up ExpressDEV (r) tool.
	//
	mDevOn = false;

	if( !mDevDates.empty()){
	    	mDevOn = true;
	    	iterDev = mDevDates.begin();

		Discount = sliceNew(nDim);
		StatePr  = sliceNew(nDim);
		StatePrL = StatePr + NodeOffset(nDim, 0, 0, 0);
    		StatePrL[0] = 1e0;

		// If the start date is a express dev date
		if ( *iterDev == TPToday())
		{
	     		devStatePr = sliceNew(nDim);

			devStatePrL = devStatePr + NodeOffset(nDim, 0, 0, 0);
			devStatePrL[0] = 1e0;

			mDevStatePr.insert(KMap(TDate, double*)::value_type(
					   TPToday(), devStatePr));

			devStatePr = NULL;

			iterDev++;
		}
	}

 
	// initial guess of forward shift
	// It will be overwritten in CalcDiscDiff to
	// ensure that the rate distribution at X=0 
	// corresponds to the forward rate
	Zt = FwdShift;

	//
	// Calibrate IR drift at each point in the forward loop
	//
	for (t = 0; t < NbTP; t++)               
	{                                   
		tpIdxCurrent = t;	// set the current time point

		//
		// Calculate the discount between t and t+1
		// using the previous offset as a first guess
		//
		CalcDiscount(t, nIRDim, Zt, mDiscountIR);

		if (debugLevel >= DEBUG_LEVEL_GEOMETRY) {
			dppLog << format("%s: discount at tpIdx %3d "
				"before solve.\n", routine, t);
			slicePrint(
        			mDiscountIR[KV_DIFF],
        			nIRDim, 
        			t,   
        			FALSE, 
        			dppLog);
		}


		//
		// Solve for the offset to calibrate zero at t+1
		// using previous offset as initial guess
		//
		zeroPrice = GetZeroPrice(KV_DIFF)[t+1];

		SolveOffset(	
			t,
			nIRDim,
			DiscountIR,
			StatePrIR,
			zeroPrice,
        		Zt,
			&DelZt);

		//
		// Store new drift offset
		//
        	Zt += DelZt;
        	mTpZCenter[t] = Zt;

		//
		// Recompute the discount factor slice between (t and t+1)
		// with offset solved.
		//
		CalcDiscount(t, nIRDim, Zt, mDiscountIR);


		//
		// Compute transition probabilities between t and t+1
		//
		KMrNTree::Update(t);


		//
		// Compute state prices at t+1 by
		// forwarding the slice to t+1 (adjoint of discounting)
		//


		sliceFw(StatePrIR, nIRDim, DiscountIR, t);


		//
		// Set the current time point to t+1 after the forward
		//
		tpIdxCurrent = t+1;

		//
		// Test for correct zero coupon price:
		// sum of all state prices should equal zero bond
		//
		sliceSpecialOper(StatePrIR, nIRDim, "sum", &p);

#ifndef	__NO_CALIB__
		if (fabs (p - zeroPrice) > STATE_PRICE_TOL) {
			throw KFailure("%s: at time point %d, %d-D IR "
				"state prices don't add up to Z (res=%lf).\n",
				routine, t, nIRDim, fabs (p - zeroPrice));
		}
#endif

		if (debugLevel > DEBUG_LEVEL_TIMELINE) {
		    dppLog << format("%s: AFTER TPIDX=%4d Z_ERROR=%14.10f\n",
				routine, t, p - zeroPrice) << endl;
		}

		//
		// Express DEV tool: store state prices in full nDim dimension.
		//
		if (mDevOn && (iterDev != mDevDates.end()))
		{
		    if (nDim > nIRDim)
		    {
		    	//
		    	// Expand DiscountIR in nDim dimensions, 
		    	// which only varies in the first nIRDim 
		    	// dimensions and is constant in the 
		    	// higher nDim-nIRDim dimensions
		    	//
		    	sliceExpand(t, nDim, nIRDim, DiscountIR, Discount);

		    	//
		    	// Forward state price in nDim.
		    	//
		    	sliceFw(StatePr, nDim, Discount, t);

			//
			// sum of all state prices should equal zero bond
			//
			sliceSpecialOper(StatePr, nDim, "sum", &p);

#ifndef	__NO_CALIB__
			if (fabs (p - zeroPrice) > STATE_PRICE_TOL) {
			throw KFailure("%s: at time point %d, EDev %d-D state "
				"prices don't add up to Z (res=%lf).\n",
				routine, t, nDim, fabs (p - zeroPrice));
			}
#endif
		    }
	    
		    //
		    // Store state price at special dev dates
		    //
		    if ( *iterDev == TPDates[t+1])
		    {
	     		devStatePr = sliceNew(nDim);

			// Store corresponding state prices
			if(nDim > nIRDim)
			    sliceUnaryOper(devStatePr,
				           nDim,
				           StatePr,
				           COPY);
			else	// nDim=nIRDim
			    sliceUnaryOper(devStatePr,
				           nIRDim,
				           StatePrIR,
				           COPY);

			mDevStatePr.insert(KMap(TDate, double*)::value_type(
					  *iterDev, devStatePr));

			devStatePr = NULL;

			iterDev++;
		    }	// if iterDev
		}   // mDevOn

	}  // for t                                         

	// Debug print for drift calibration
	if (debugLevel > DEBUG_LEVEL_TIMELINE) {
		dppLog << format("%s: Drift calibration center offset "
				"(QHi=%lf  QLo=%lf FSh=%lf):\n",
				routine, mQHi, mQLo, mFSh);
		for (t = 0; t < NbTP; t++)               
		{
			dppLog << format(" %4d   %14.10f",
	    				t, mTpZCenter[t]) << endl;
		}
	}

	// Free tmp memory
	sliceDelete(StatePr);
	sliceDelete(StatePrIR);
	sliceDelete(Discount);

    }
    catch (KFailure) {
	// Free tmp memory
	sliceDelete(StatePr);
	sliceDelete(StatePrIR);
	sliceDelete(Discount);

	throw KFailure("%s: failed.\n", routine);
    }
}


//--------------------------------------------------------------
// Calculates the discount factor applicable between
// tpIdx and tpIdx+1.

void
KPirTree::CalcDiscount(
	int tpIdx,			// (I) time point index
	int nIRDim,			// (I) IR dimensions
	double Zt,			// (I) center offset 
	KMap(int, double*) &discount)	// (O) discount slice
{
static	char	routine[] = "KPirTree::CalcDiscount";

	int	t = tpIdx;		// more convenient ...
	double	ZRatio1,		// Zero coupon ratios for idx curve 1
		ZRatio2,		// Zero coupon ratios for idx curve 2 
		discDetZero;		// zero price of deterministic curve

	double	*discDetS,		// deterministic discount slice
		*discDiffS,		// diffuse discount slice
		*discIndxS;		// index discount slice

    try {

	//
	// Calculate the discount
	//
	CalcDiscDiff(
		tpIdx,			// (I) time point index
		nIRDim,			// (I) IR dimensions
		Zt,			// (I) center offset 
		discount[KV_DIFF]);


	//
	// Now apply the zero ratio zero method
	// to compute the index curves discount factors
	//

	KMap(int, double*)::iterator itIdx1 = mTpFRates.find(KV_IDX1);
	if (itIdx1 != mTpFRates.end()) {
		ZRatio1 = (1.+ GetForwardRate(KV_DIFF)[t])
		  	  /(1.+ GetForwardRate(KV_IDX1)[t]);

		discDiffS = discount[KV_DIFF];
		discIndxS = discount[KV_IDX1];

		sliceSpecialOper(
			discIndxS,
			nIRDim,
			"c*s",
			ZRatio1,
			discDiffS);
	}

	KMap(int, double*)::iterator itIdx2 = mTpFRates.find(KV_IDX2);
	if (itIdx2 != mTpFRates.end()) {
		ZRatio2 = (1.+ GetForwardRate(KV_DIFF)[t])
		  	  /(1.+ GetForwardRate(KV_IDX2)[t]);
		discDiffS = discount[KV_DIFF];
		discIndxS = discount[KV_IDX2];

		sliceSpecialOper(
			discIndxS,
			nIRDim,
			"c*s",
			ZRatio2,
			discDiffS);
	}



	//
	// Deterministic curve
	//
	KMap(int, double*)::iterator itIdx0 = mTpFRates.find(KV_DET);
	if (itIdx0 != mTpFRates.end()) {
		discDetS    = discount[KV_DET];
		discDetZero = GetZeroPrice(KV_DET)[t+1]/GetZeroPrice(KV_DET)[t];

		sliceScalarOper(
			discDetS,
			nIRDim,
			discDetZero,
			COPY);
	}


	

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}



//--------------------------------------------------------------
// Solves for the offset to calibrate the zero price at tpIdx+1.

void
KPirTree::SolveOffset(
	int tpIdx,		// (I) time point index
	int nIRDim,		// (I) IR dimensions
	double *discount,	// (I) slice discount bet t and t+1
	double *statePr,	// (I) slice state price at t
	double zeroPrice,	// (I) zero price at t+1
	double Zt,		// (I) intial offset
	double *DelZt)		// (O) incremental offset 

{
static	char	routine[] = "KPirTree::SolveOffset";

    try {
#ifndef	__NO_CALIB__
	SolveOffsetDiff(
		tpIdx,		// (I) time point index
		nIRDim,		// (I) IR dimensions
		discount,	// (I) slice discount bet t and t+1
		statePr,	// (I) slice state price at t
		zeroPrice,	// (I) zero price at t+1
		Zt,		// (I) intial offset
		DelZt);		// (O) incremental offset 
#else
	*DelZt = 0e0;
#endif




    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }

}





//--------------------------------------------------------------
// Calculates the discount factor applicable between
// tpIdx and tpIdx+1.
//

void
KPirTree::CalcDiscDiff(
	int tpIdx,			// (I) time point index
	int nIRDim,			// (I) IR dimensions
	double Zt,			// (I) center offset 
	double *Discount)		// (O) discount slice
{
static	char	routine[] = "KPirTree::CalcDiscDiff";


	int	t = tpIdx,		// more convenient ...
		i, j, k,		// dimension 1, 2, 3 indices
		Mid;

	double	du;
        double	pJump00,
        	pJump10, pJump11,
        	pJump20, pJump21, pJump22;

	double  QLeft   = this->mQLo;
	double  QRight  = this->mQHi;
	double  FwdShift = this->mFSh;

	double  FwdRate;	// Fwd rate 
	double  FwdRateA;	// Fwd rate adjusted 
	double  MLeft, SLeft;	// Multiple coeff and shift for grid pt
	double  MRight, SRight;
	double  VolBbq;		// Sigma used in bone mapping 
	
	double  QSh;		// q parameter for shift adjustment 
	
	double  Zidx;		// Zt index i,j,k adjusted 
	double  Grid;		// Grid points  

	double	RateJumpLeft,	// Rate space jump size for last dimension
		RateJumpRight,
		RateJumpLeft2,
		RateJumpRight2,
		RateJumpLeft5,
		RateJumpRight5;

	int     *Top1         = mTop1;
	int     *Bottom1      = mBottom1;
	int     **Top2        = mTop2;
	int     **Bottom2     = mBottom2;
	int     ***Top3       = mTop3;
	int     ***Bottom3    = mBottom3;
	int     *OutTop1      = mOutTop1;
	int     *OutBottom1   = mOutBottom1;
	int     **OutTop2     = mOutTop2;
	int     **OutBottom2  = mOutBottom2;
	int     ***OutTop3    = mOutTop3;
	int     ***OutBottom3 = mOutBottom3;


	double	*DiscountL;	// local slice pointer

    try {

	//
        // Precompute jumps 
	//

        du = sqrt (JUMPCOEFF * LengthJ[t-1]);
        
	switch (nIRDim) {
	case 3:
        	pJump20 = mAweight[t-1][3] * du;
        	pJump21 = mAweight[t-1][4] * du;
        	pJump22 = mAweight[t-1][5] * du;
	case 2:
        	pJump10 = mAweight[t-1][1] * du;
        	pJump11 = mAweight[t-1][2] * du;
	case 1:
        	pJump00 = mAweight[t-1][0] * du;
	}


	/*dppLog << format("%14.8f\n%14.8f\t%14.8f\n%14.8f\t%14.8f\t%14.8f\n",
		pJump00,
		pJump10, pJump11,
		pJump20, pJump21, pJump22);*/


	//
	// Calc needed constants
	//
	FwdRate = GetForwardRate(KV_DIFF)[t];

	// Avoid zero forward rate, which is allowed in normal case but
        // canceled out after the mapping.  Set the minimum as 0.1bp
	SHIFT_ZERO(FwdRate);


        FwdRateA  = FwdRate / (1. + FwdShift);


	// Vol backbone factor
	//
        VolBbq  = GetVolBbq(t);


        // Set initial Zt
        // Make sure that the rate distribution at X=0 
        // corresponds to the forward rate
        if (t == 0)
        {
            QSh = (FwdShift > 0) ? QRight : QLeft;
            if (IS_Q(QSh))
            {
                Zt = log(1. + QSh * FwdShift) / (QSh * VolBbq);
            }
            else
            {
                Zt = FwdShift / VolBbq;
            }
        }

	//
	// Dimension specific 
	//

	//----------------------------------------------
	// 1 F
	//----------------------------------------------

	if (nIRDim == 1) {

        DiscountL = Discount + NodeOffset(1, 0, 0, t);


        if (IS_Q(QLeft))
        {
            MLeft        = FwdRateA / QLeft;
            SLeft        = 1 + FwdRateA - FwdRateA / QLeft;     // 1 + r 
            RateJumpLeft = exp(QLeft * VolBbq * pJump00);
        }
        else
        {
            MLeft        = FwdRateA;
            SLeft        = 1 + FwdRateA;
            RateJumpLeft = FwdRateA * VolBbq * pJump00;
        }
        if (IS_Q(QRight))
        {
            MRight        = FwdRateA / QRight;
            SRight        = 1 + FwdRateA - FwdRateA / QRight;  // 1 + r      
            RateJumpRight = exp(QRight * VolBbq * pJump00);
        }
        else
        {
            MRight        = FwdRateA;
            SRight        = 1 + FwdRateA;
            RateJumpRight = FwdRateA * VolBbq * pJump00;
        }


        /* 
        *   Set up the grid points. 
        */

        /* LEFT part of distribution */

        Zidx  = Zt; 
        Mid   = (int) ceil(-Zidx / pJump00) - 1;
        Mid   = MIN ( MAX (Mid, Bottom1[t] - 1), Top1[t]);
        Zidx += (pJump00) * Bottom1[t];



        if (IS_Q(QLeft))
        {
            Grid = MLeft * exp (QLeft * VolBbq * Zidx);
        
            for (i = Bottom1[t]; i <= Mid; i++)
            {
                DiscountL[i] = 1. / (SLeft + Grid);
                
                Grid *= RateJumpLeft;
            }           
        }
        else
        {
            Grid = MLeft * VolBbq * Zidx;

            for (i = Bottom1[t]; i <= Mid; i++)
            {
                DiscountL[i] = 1. / (SLeft + Grid);
                
                Grid += RateJumpLeft;
            }
        }

        /* Right part of distribution */

        Zidx = Zt 
             + (pJump00) * (Mid + 1);

        if (IS_Q(QRight))
        {
            Grid = MRight * exp (QRight * VolBbq * Zidx);
        
            for (i = Mid + 1; i <= Top1[t]; i++)
            {
                DiscountL[i] = 1. / (SRight + Grid);
                
                Grid *= RateJumpRight;
            }           
        }
        else
        {
            Grid = MRight * VolBbq * Zidx;

            for (i = Mid +1; i <= Top1[t]; i++)
            {
                DiscountL[i] = 1. / (SRight + Grid);
                
                Grid += RateJumpRight;
            }
        }


	//----------------------------------------------
	// 2 F
	//----------------------------------------------
	} else if (nIRDim == 2) {


        if (IS_Q(QLeft))
        {
            MLeft         = FwdRateA / QLeft;
            SLeft         = 1 + FwdRateA - FwdRateA / QLeft;      
            RateJumpLeft2 = exp(QLeft * VolBbq * pJump11);
        }
        else
        {
            MLeft         = FwdRateA;
            SLeft         = 1 + FwdRateA;
            RateJumpLeft2 = FwdRateA * VolBbq * pJump11;
        }
        if (IS_Q(QRight))
        {
            MRight         = FwdRateA / QRight;
            SRight         = 1 + FwdRateA - FwdRateA / QRight;      
            RateJumpRight2 = exp(QRight * VolBbq * pJump11);
        }
        else
        {
            MRight         = FwdRateA;
            SRight         = 1 + FwdRateA;
            RateJumpRight2 = FwdRateA * VolBbq * pJump11;
        }



        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {
            /* LEFT part of distribution */

            Zidx  = Zt 
                  + (pJump00 + pJump10) * i;
            Mid   = (int) ceil(-Zidx / pJump11) - 1;
            Mid   = MIN ( MAX (Mid, Bottom2[t][i] - 1), Top2[t][i]);
            Zidx += (pJump11) * Bottom2[t][i];

            DiscountL = Discount + NodeOffset (2, i, 0, t);
    
            if (IS_Q(QLeft))
            {
                Grid = MLeft * exp (QLeft * VolBbq * Zidx);
            
                for (j = Bottom2[t][i]; j <= Mid; j++)
                {
                    DiscountL[j] = 1. / (SLeft + Grid);
                
                    Grid *= RateJumpLeft2;
                }
            }
            else
            {
                Grid = MLeft * VolBbq * Zidx;

                for (j = Bottom2[t][i]; j <= Mid; j++)
                {
                    DiscountL[j] = 1. / (SLeft + Grid);
                
                    Grid += RateJumpLeft2;
                }
            }

            /* Right part of distribution */

            Zidx = Zt 
                 + (pJump00 + pJump10) * i
                 + (pJump11) * (Mid + 1);

            if (IS_Q(QRight))
            {
                Grid = MRight * exp (QRight * VolBbq * Zidx);
            
                for (j = Mid + 1; j <= Top2[t][i]; j++)
                {
                    DiscountL[j] = 1. / (SRight + Grid);
                
                    Grid *= RateJumpRight2;
                }
            }
            else
            {
                Grid = MRight * VolBbq * Zidx;

                for (j = Mid + 1; j <= Top2[t][i]; j++)
                {
                    DiscountL[j] = 1. / (SRight + Grid);
                
                    Grid += RateJumpRight2;
                }
            }

        }  /* for i */



	//----------------------------------------------
	// 3 F
	//----------------------------------------------

	} else if (nIRDim == 3) {


        if (IS_Q(QLeft))
        {
            MLeft         = FwdRateA / QLeft;
            SLeft         = 1 + FwdRateA - FwdRateA / QLeft;      
            RateJumpLeft5 = exp(QLeft * VolBbq * pJump22);
        }
        else
        {
            MLeft         = FwdRateA;
            SLeft         = 1 + FwdRateA;
            RateJumpLeft5 = FwdRateA * VolBbq * pJump22;
        }
        if (IS_Q(QRight))
        {
            MRight         = FwdRateA / QRight;
            SRight         = 1 + FwdRateA - FwdRateA / QRight;      
            RateJumpRight5 = exp(QRight * VolBbq * pJump22);
        }
        else
        {
            MRight         = FwdRateA;
            SRight         = 1 + FwdRateA;
            RateJumpRight5 = FwdRateA * VolBbq * pJump22;
        }



        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {
            for (j = Bottom2[t][i]; j <= Top2[t][i]; j++)
            {
                /* LEFT part of distribution */

                Zidx  = Zt 
                      + (pJump00 + pJump10 + pJump20) * i
                      + (pJump11 + pJump21) * j;
                Mid   = (int) ceil(-Zidx / pJump22) - 1;
                Mid   = MIN ( MAX (Mid, Bottom3[t][i][j] - 1), Top3[t][i][j]);
                Zidx += (pJump22) * Bottom3[t][i][j];
            
                DiscountL = Discount + NodeOffset (3, i, j, t);
                
                if (IS_Q(QLeft))
                {
                    Grid = MLeft * exp (QLeft * VolBbq * Zidx);
    
                    for (k = Bottom3[t][i][j]; k <= Mid; k++)
                    {
                        DiscountL[k] = 1. / (SLeft + Grid);
                
                        Grid *= RateJumpLeft5;
                    }
                }
                else
                {
                    Grid = MLeft * VolBbq * Zidx;

                    for (k = Bottom3[t][i][j]; k <= Mid; k++)
                    {
                        DiscountL[k] = 1. / (SLeft + Grid);
                
                        Grid += RateJumpLeft5;
                    }
                }

                /* Right part of distribution */

                Zidx = Zt 
                     + (pJump00 + pJump10 + pJump20) * i
                     + (pJump11 + pJump21) * j
                     + (pJump22) * (Mid + 1);

                if (IS_Q(QRight))
                {
                    Grid = MRight * exp (QRight * VolBbq * Zidx);
    
                    for (k = Mid + 1; k <= Top3[t][i][j]; k++)
                    {
                        DiscountL[k] = 1. / (SRight + Grid);
                
                        Grid *= RateJumpRight5;
                    }
                }
                else
                {
                    Grid = MRight * VolBbq * Zidx;

                    for (k = Mid + 1; k <= Top3[t][i][j]; k++)
                    {
                        DiscountL[k] = 1. / (SRight + Grid);
                
                        Grid += RateJumpRight5;
                    }
                }

            }  /* for j */     
        
        }  /* for i */
       
 
	} // if (nIRDim == 
                



    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}



//--------------------------------------------------------------
// Solves for the offset to calibrate the zero price at tpIdx+1.
//

void
KPirTree::SolveOffsetDiff(
	int tpIdx,		// (I) time point index
	int nIRDim,		// (I) IR dimensions
	double *Discount,	// (I) slice discount bet t and t+1
	double *StatePr,	// (I) slice state price at t
	double zeroPrice,	// (I) zero price at t+1
	double Zt,		// (I) intial offset
	double *DelZt)		// (O) incremental offset 

{
static	char	routine[] = "KPirTree::SolveOffsetDiff";

	int	t = tpIdx,		// more convenient ...
		i, j, k,		// dimension 1, 2, 3 indices
		Mid;

	double	du;
        double	pJump00,
        	pJump10, pJump11,
        	pJump20, pJump21, pJump22;

	double  QLeft   = this->mQLo;
	double  QRight  = this->mQHi;
	double  FwdShift = this->mFSh;

	double  FwdRate;	// Fwd rate 
	double  FwdRateA;       // Fwd rate adjusted 
	double  VolBbq;         // Sigma used in bone mapping

	double  Zidx;		// Zt index i,j,k adjusted


	double	P[4],		// expansion
		Q, U,		// Q and U weight coeff for Z-polynomial
		D0,		// Consecutive powers of disc fact
		x,
		aL00,		// Coeff of Z-polynomial
		aL10, aL11,
		aL20, aL21, aL22,
		aR00,
		aR10, aR11,
		aR20, aR21, aR22;


	int     *Top1         = mTop1;
	int     *Bottom1      = mBottom1;
	int     **Top2        = mTop2;
	int     **Bottom2     = mBottom2;
	int     ***Top3       = mTop3;
	int     ***Bottom3    = mBottom3;
	int     *OutTop1      = mOutTop1;
	int     *OutBottom1   = mOutBottom1;
	int     **OutTop2     = mOutTop2;
	int     **OutBottom2  = mOutBottom2;
	int     ***OutTop3    = mOutTop3;
	int     ***OutBottom3 = mOutBottom3;


	double	*DiscountL,	// local slice pointer
		*StatePrL;

    try {

	//
        // Precompute jumps 
	//

        du = sqrt (JUMPCOEFF * LengthJ[t-1]);
        
	switch (nIRDim) {
	case 3:
        	pJump20 = mAweight[t-1][3] * du;
        	pJump21 = mAweight[t-1][4] * du;
        	pJump22 = mAweight[t-1][5] * du;
	case 2:
        	pJump10 = mAweight[t-1][1] * du;
        	pJump11 = mAweight[t-1][2] * du;
	case 1:
        	pJump00 = mAweight[t-1][0] * du;
	}

	/*dppLog << format("%14.8f\n%14.8f\t%14.8f\n%14.8f\t%14.8f\t%14.8f\n",
		pJump00,
		pJump10, pJump11,
		pJump20, pJump21, pJump22);*/



	//
	// Calc needed constants
	//
	FwdRate = GetForwardRate(KV_DIFF)[t];

	// Avoid zero forward rate, which is allowed in normal case but
        // canceled out after the mapping.  Set the minimum as 0.1bp
	SHIFT_ZERO(FwdRate);

        FwdRateA  = FwdRate / (1. + FwdShift);


	// Vol backbone factor
	//
        VolBbq  = GetVolBbq(t);

	//
	//
	//


        for (i = 0; i < 3; i++)                     
        {
            P[i] = 0e0; 
            
        }


        Q = QLeft * VolBbq;
        U = (FwdRateA * (1. - QLeft) - QLeft) * VolBbq; 
    
        aL00 =  1.;
        aL10 = -Q;          aL11 = -U;
        aL20 =  0.5*Q*Q;    aL21 =  1.5*U*Q;    aL22 =   U*U;

        Q = QRight * VolBbq;
        U = (FwdRateA * (1. - QRight) - QRight) * VolBbq; 
    
        aR00 =  1.;
        aR10 = -Q;          aR11 = -U;
        aR20 =  0.5*Q*Q;    aR21 =  1.5*U*Q;    aR22 =   U*U;

        Zidx = Zt; 


	//----------------------------------------------
	// 1 F
	//----------------------------------------------
	if (nIRDim == 1) {


	StatePrL  = StatePr  + NodeOffset(1, 0, 0, t);
	DiscountL = Discount + NodeOffset(1, 0, 0, t);


        Mid  = (int) ceil(-Zidx / pJump00) - 1;
        Mid  = MIN ( MAX (Mid, Bottom1[t] - 1), Top1[t]);

        for (i = Bottom1[t]; i <= Mid; i++)
        {
            x   = DiscountL[i];
            D0  = StatePrL[i] * x;
            P[0] += D0;

            P[1] += (aL10 + aL11 * x) * D0;

            P[2] += (aL20 + (aL21 + aL22 * x) * x) * D0;                  
        }

        for (i = Mid + 1; i <= Top1[t]; i++)
        {
            x   = DiscountL[i];
            D0  = StatePrL[i] * x;
            P[0] += D0;

            P[1] += (aR10 + aR11 * x) * D0;

            P[2] += (aR20 + (aR21 + aR22 * x) * x) * D0;                  
        }


	//----------------------------------------------
	// 2 F
	//----------------------------------------------
	} else if (nIRDim == 2) {


        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {    
            StatePrL  = StatePr  + NodeOffset(2, i, 0, t);
            DiscountL = Discount + NodeOffset(2, i, 0, t);
                                                
            Zidx = Zt 
                 + (pJump00 + pJump10) * i;
            Mid  = (int) ceil(-Zidx / pJump11) - 1;
            Mid  = MIN ( MAX (Mid, Bottom2[t][i] - 1), Top2[t][i]);

            for (j = Bottom2[t][i]; j <= Mid ; j++)
            {
                
                x   = DiscountL[j];
                D0  = StatePrL[j] * x;
                P[0] += D0;

                P[1] += (aL10 + aL11 * x) * D0;
            
                P[2] += (aL20 + (aL21 + aL22 * x) * x) * D0;                  

            }

            for (j = Mid + 1; j <= Top2[t][i] ; j++)
            {
                
                x   = DiscountL[j];
                D0  = StatePrL[j] * x;
                P[0] += D0;

                P[1] += (aR10 + aR11 * x) * D0;
            
                P[2] += (aR20 + (aR21 + aR22 * x) * x) * D0;                  

            }
        }  /* for i */


	//----------------------------------------------
	// 3 F
	//----------------------------------------------
	} else if (nIRDim == 3) {

        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {
            for (j = Bottom2[t][i]; j <= Top2[t][i]; j++)
            {
                StatePrL  = StatePr  + NodeOffset(3, i, j, t);
                DiscountL = Discount + NodeOffset(3, i, j, t);
                    
                Zidx = Zt 
                     + (pJump00 + pJump10 + pJump20) * i
                     + (pJump11 + pJump21) * j;
                Mid  = (int) ceil(-Zidx / pJump22) - 1;
                Mid  = MIN ( MAX (Mid, Bottom3[t][i][j] - 1), Top3[t][i][j]);

                for (k = Bottom3[t][i][j]; k <= Mid; k++)
                {

                    x   = DiscountL[k];
                    D0  = StatePrL[k] * x;
                    P[0] += D0;

                    P[1] += (aL10 + aL11 * x) * D0;
                              
                    P[2] += (aL20 + (aL21 + aL22 * x) * x) * D0;                  

                }

                for (k = Mid + 1; k <= Top3[t][i][j]; k++)
                {

                    x   = DiscountL[k];
                    D0  = StatePrL[k] * x;
                    P[0] += D0;

                    P[1] += (aR10 + aR11 * x) * D0;
                              
                    P[2] += (aR20 + (aR21 + aR22 * x) * x) * D0;                  

                }

            }  /* for j */
        }  /* for i */
        

	}



	//----------------------------------------------
	// Solve for zero price at t+1
	//----------------------------------------------

	P[0] -= zeroPrice;

	if (debugLevel >= DEBUG_LEVEL_DRIFT) {
		dppLog << format("%s: TPIDX=%4d "
			"Z=%14.10f P[0-2] %14.6e %14.6e %14.6e",
			routine, tpIdx, zeroPrice, P[0], P[1], P[2]) << endl;
	}

	if ( DrlNRPoly(
		0e0,
		P, 
		2,
		DelZt) != SUCCESS) {
		throw KFailure("%s: Failed at TPIDX=%4d P[0-2] "
			"Z=%14.10f P[]=%14.6e %14.6e %14.6e DelZT=%14.10f\n",
			routine, tpIdx, zeroPrice, P[0], P[1], P[2],
			*DelZt);
	}

	if (debugLevel >= DEBUG_LEVEL_DRIFT) {
		dppLog << format("%s: TPIDX=%4d "
			"Z=%14.10f P[0-2] %14.6e %14.6e %14.6e DelZT=%14.10f",
			routine, tpIdx, zeroPrice, P[0], P[1], P[2],
			*DelZt) << endl;
	}

	if (debugLevel >= DEBUG_LEVEL_GEOMETRY) {

		dppLog << "Input State Price:" << endl;
		slicePrint(
        		StatePr,
        		nIRDim, 
        		t,   
        		FALSE, 
        		dppLog);

		dppLog << "Discount Factor:" << endl;
		slicePrint(
        		Discount,
        		nIRDim, 
        		t,   
        		FALSE, 
        		dppLog);

		dppLog << endl;
		dppLog << format("\tCENTER_OFFSET: %15.10f\n", *DelZt);
	}


	if (fabs(*DelZt*VolBbq) > 1e0) {
		throw KFailure("%s: problem in calculation "
			"of drift at timepoint %4d (DelZt = %lf).\n",
			routine, tpIdx, *DelZt);
        }

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}



