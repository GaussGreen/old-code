/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      C. Daher, David Liu
 ************************************************************************/
#include "krate.h"		// Class definition
#include "kutilios.h"

extern	"C" {
#include "fltrate.h"
#include "duration.h"         /* GtoMoneyMarketModDuration */
#include "convex.h"           /* GtoBondConvexity */
#include "ldate.h"            /* GtoDtFwdAny GtoDayCountFraction GTO_ACT_365 */
#include "interp.h"           /* GTO_LINEAR_INTERP */
#include "stub.h"             /* GTO_STUB_SIMPLE */
#include "swapadj.h"          /* GtoCUPSCorrAdj */
#include "yearfrac.h"         /* GtoStringToDayCount */
#include "zr2coup.h"          /* GtoZerosToCouponsPoint */
#include "zr2fwd.h"           /* GtoForwardFromZcurve */

#include "drltime.h"

#include "crcrv.h"            /* CrxFwdParCDSSpread   */
};


//***************************************************************
//
//	KRate functions.
//
//***************************************************************
#define	MAX_HOL	256


static	void _TDateAdjIntvlCopy(
	TDateAdjIntvl *target,		// (I) 
	const TDateAdjIntvl* source)	// (O) 
{
	target->interval = source->interval;
	target->isBusDays = source->isBusDays;
	if (source->holidayFile != NULL) {
	    if (target->holidayFile == NULL) {
		target->holidayFile = NEW_ARRAY(char,MAX_HOL);
		ASSERT_OR_THROW(target->holidayFile != NULL);
	    }
	    strncpy(target->holidayFile, source->holidayFile, MAX_HOL-1);
	} else {
	    target->holidayFile = NULL;
	}
	target->badDayConv = source->badDayConv;
}

static	bool	_TDateAdjIntvlIsLess(
	const TDateInterval *intvl1,	// (I) 
	const TDateInterval *intvl2)	// (I) 
{
        double  x1, x2;
        TDateInterval   iv1, iv2;
 
        iv1 = *intvl1;
        IF_FAILED_THROW( GtoDateIntervalToYears(
                &iv1, &x1));
 
        iv2 = *intvl2;
        IF_FAILED_THROW( GtoDateIntervalToYears(
                &iv2, &x2));
 
        return(x1 < x2);
}


// Convert all the TDateInterval into economical form of 'M' or 'D'
// for comparison of equals.
// For TDateInterval in "W" and "D", convert to form of "D" and then
// compare.  
// This conversion assumes that there is no equivalence between "W", "D"
// and "A", "S", "Q", and "M". 
//
static	TDateInterval	_TDateIntvlEquivlant(
	TDateInterval intvl)	// (I) 
{
		
	char	prd_type;
	TDateInterval intvlEqv;


	// Default, not used in comparison
	intvlEqv.flag = 0;

	prd_type = intvl.prd_typ;
	switch (prd_type)
	{
	    case 'A':
		intvlEqv.prd_typ = 'M';
		intvlEqv.prd = (intvl.prd)*12;
		break;
	
	    case 'S':
		intvlEqv.prd_typ = 'M';
		intvlEqv.prd = (intvl.prd)*6;
		break;
	
	    case 'Q':
		intvlEqv.prd_typ = 'M';
		intvlEqv.prd = (intvl.prd)*3;
		break;
	
	    case 'M':
		intvlEqv.prd_typ = 'M';
		intvlEqv.prd = intvl.prd;
		break;
	    case 'W':
		intvlEqv.prd_typ = 'D';
		intvlEqv.prd = (intvl.prd)*7;
		break;
	    case 'D':
		intvlEqv.prd_typ = 'D';
		intvlEqv.prd = intvl.prd;
		break;
	    default:
		throw KFailure("%s: invalid TDateInverval for comparison. "
			"Only 'A', 'S', 'Q', 'M', 'W' and 'D' are allowd.\n",
			DrlTDateIntervalPrint(NULL, intvl));
	}

	return intvlEqv;
}




static	bool	_TDateAdjIntvlIsEqual(
	const TDateInterval *intvl1,	// (I) 
	const TDateInterval *intvl2)	// (I) 
{
static  char    routine[] = "_TDateAdjIntvlIsEqual";

        TDateInterval   iv1, iv2;
 
	// Convert to economical form
	//
        iv1 = _TDateIntvlEquivlant(*intvl1);
        iv2 = _TDateIntvlEquivlant(*intvl2);
 
	// GtoIntervalEqualsInterval compares both the number 
	// of periods and period type without conversion.  
	// So 0A != 0D.
	//
	if (iv1.prd == 0 && iv2.prd == 0)
		return(true);
	else 
        	return(GtoIntervalEqualsInterval(&iv1, &iv2));
}



//---------------------------------------------------------------


KRate::KRate(
	const KDateInterval &mat,	// (I) rate maturity
	const KDateInterval &freq,	// (I) rate payment freq
	KDayCc dayCc,			// (I) rate day count conv
	const KDateInterval &spotOffset,// (I) spot offset
	double spread,			// (I) spread
	double weight)			// (I) weight
{
	mFloatRate.matInterval = (TDateInterval) mat;
	mFloatRate.payInterval = (TDateInterval) freq;
	mFloatRate.dayCountConv = (long) dayCc;

	mFloatRate.spotOffset.holidayFile = NULL;
	_TDateAdjIntvlCopy(
		&mFloatRate.spotOffset, 
		(const TDateAdjIntvl*) spotOffset);

	mFloatRate.spread = spread;
	mFloatRate.weight = weight;

	mFloatRate.rateType = GTO_SIMPLE_BASIS;

	mCurveName = K_DEFAULT_NAME;
}


KRate::KRate(
	const String &curveName,	// (I) curve name
	const KDateInterval &mat,	// (I) rate maturity
	const KDateInterval &freq,	// (I) rate payment freq
	KDayCc dayCc,			// (I) rate day count conv
	const KDateInterval &spotOffset,// (I) spot offset
	double spread,			// (I) spread
	double weight)			// (I) weight
{
	mFloatRate.matInterval = (TDateInterval) mat;
	mFloatRate.payInterval = (TDateInterval) freq;
	mFloatRate.dayCountConv = (long) dayCc;

	mFloatRate.spotOffset.holidayFile = NULL;
	_TDateAdjIntvlCopy(
		&mFloatRate.spotOffset, 
		(const TDateAdjIntvl*) spotOffset);

	mFloatRate.spread = spread;
	mFloatRate.weight = weight;

	mFloatRate.rateType = GTO_SIMPLE_BASIS;

	mCurveName = curveName;
}



KRate::KRate(
	const String &name,		// (I) name
	const String &curveName,	// (I) curve name
	const KDateInterval &mat,	// (I) rate maturity
	const KDateInterval &freq,	// (I) rate payment freq
	KDayCc dayCc,			// (I) rate day count conv
	const KDateInterval &spotOffset,// (I) spot offset
	double spread,			// (I) spread
	double weight)			// (I) weight
{
	mFloatRate.matInterval = (TDateInterval) mat;
	mFloatRate.payInterval = (TDateInterval) freq;
	mFloatRate.dayCountConv = (long) dayCc;

	mFloatRate.spotOffset.holidayFile = NULL;
	_TDateAdjIntvlCopy(
		&mFloatRate.spotOffset, 
		(const TDateAdjIntvl*) spotOffset);

	mFloatRate.spread = spread;
	mFloatRate.weight = weight;

	mFloatRate.rateType = GTO_SIMPLE_BASIS;

	mCurveName = curveName;

	SetName(name);
}


KRate::KRate(
	const String 	&curveName,	// (I) curve name
	const TDate	startDate,	// (I) start date
	const TDate	endDate,	// (I) end date
	KDayCc 		dayCc,		// (I) rate day count conv
	const KDateInterval &spotOffset,// (I) spot offset
	double spread,			// (I) spread
	double weight)			// (I) weight
{
	KDateInterval	mat(endDate, startDate);

	mFloatRate.matInterval = (TDateInterval) mat;
	mFloatRate.payInterval = (TDateInterval) mat;
	mFloatRate.dayCountConv = (long) dayCc;

	mFloatRate.spotOffset.holidayFile = NULL;
	_TDateAdjIntvlCopy(
		&mFloatRate.spotOffset, 
		(const TDateAdjIntvl*) spotOffset);

	mFloatRate.spread = spread;
	mFloatRate.weight = weight;

	mFloatRate.rateType = GTO_SIMPLE_BASIS;

	mCurveName = curveName;
}




//---------------------------------------------------------------


KRate::KRate()
{
	TDateInterval	matInterval;
	TDateInterval	payInterval;
	long		dayCountConv = -1L;
	long		spotOffsetDays = 0L;
	double		spread = 0e0;
	double		weight = 0e0;

	mFloatRate.spotOffset.holidayFile = NULL;

	IF_FAILED_THROW( GtoYearsToDateInterval(
		0e0,
		&matInterval));

	IF_FAILED_THROW( GtoYearsToDateInterval(
		0e0,
		&payInterval));

	mFloatRate.spotOffset.holidayFile = NULL;
	IF_FAILED_THROW( GtoFloatRateSet(
		&matInterval,
		&payInterval,
		dayCountConv,
		spotOffsetDays,
		spread,
		weight,
		&mFloatRate));

	mFloatRate.rateType = GTO_SIMPLE_BASIS;

	mCurveName = K_DEFAULT_NAME;

}



//---------------------------------------------------------------


KRate::KRate(const KRate& rate)
	: KVPAtom(rate)
{
	mFloatRate.matInterval = rate.mFloatRate.matInterval;
	mFloatRate.payInterval = rate.mFloatRate.payInterval;
	mFloatRate.dayCountConv = rate.mFloatRate.dayCountConv;

	mFloatRate.spotOffset.holidayFile = NULL;
	_TDateAdjIntvlCopy(
		&mFloatRate.spotOffset, 
		&rate.mFloatRate.spotOffset);

	mFloatRate.spread = rate.mFloatRate.spread;
	mFloatRate.rateType = rate.mFloatRate.rateType;
	mFloatRate.weight = rate.mFloatRate.weight;

	mCurveName = rate.mCurveName;
}

//---------------------------------------------------------------



KRate&
KRate::operator=(const KRate& rate)
{

    if (this != &rate )
    {
	mFloatRate.matInterval = rate.mFloatRate.matInterval;
	mFloatRate.payInterval = rate.mFloatRate.payInterval;
	mFloatRate.dayCountConv = rate.mFloatRate.dayCountConv;

	_TDateAdjIntvlCopy(
		&mFloatRate.spotOffset, 
		&rate.mFloatRate.spotOffset);

	mFloatRate.spread = rate.mFloatRate.spread;
	mFloatRate.rateType = rate.mFloatRate.rateType;
	mFloatRate.weight = rate.mFloatRate.weight;

	mCurveName = rate.mCurveName;

    }

    return(*this);

}


//---------------------------------------------------------------

bool
KRate::operator<(const KRate& rate) const
{
	if (!(mCurveName < rate.mCurveName))
			return (FALSE);
	if (!_TDateAdjIntvlIsLess(
		&mFloatRate.matInterval, &rate.mFloatRate.matInterval))
			return (FALSE);
	if (!_TDateAdjIntvlIsLess(
		&mFloatRate.payInterval, &rate.mFloatRate.payInterval))
			return (FALSE);
	if (!(mFloatRate.dayCountConv < rate.mFloatRate.dayCountConv))
			return(FALSE);

	return (TRUE);
}


//---------------------------------------------------------------

bool
KRate::operator==(const KRate &rate) const
{
	if (mCurveName != rate.mCurveName)
			return (FALSE);
	if (!_TDateAdjIntvlIsEqual(
		&mFloatRate.matInterval, &rate.mFloatRate.matInterval))
			return (FALSE);
	if (!_TDateAdjIntvlIsEqual(
		&mFloatRate.payInterval, &rate.mFloatRate.payInterval))
			return (FALSE);

	if (mFloatRate.spotOffset.isBusDays != 
	    rate.mFloatRate.spotOffset.isBusDays)
			return (FALSE);
	if (strcmp(mFloatRate.spotOffset.holidayFile, 
	    	   rate.mFloatRate.spotOffset.holidayFile))
			return (FALSE);
	if (mFloatRate.spotOffset.badDayConv != 
	    rate.mFloatRate.spotOffset.badDayConv)
			return (FALSE);

	if (mFloatRate.dayCountConv != rate.mFloatRate.dayCountConv)
			return(FALSE);
	if (mFloatRate.spread != rate.mFloatRate.spread)
			return(FALSE);
	if (mFloatRate.weight != rate.mFloatRate.weight)
			return(FALSE);
	if (mFloatRate.rateType != rate.mFloatRate.rateType)
			return(FALSE);
	

	return (TRUE);
}

//---------------------------------------------------------------

int
KRate::IsSimple() const
{
	return _TDateAdjIntvlIsEqual(
		&mFloatRate.matInterval, &mFloatRate.payInterval);
}


//---------------------------------------------------------------
void 
KRate::SetFixedRate(double fixedRate)
{

	TDateInterval	matInterval;
	TDateInterval	payInterval;
	long		dayCountConv = -1L;
	long		spotOffsetDays = 0L;
	double		weight = 0e0;

	mFloatRate.spotOffset.holidayFile = NULL;

	IF_FAILED_THROW( GtoYearsToDateInterval(
		0e0,
		&matInterval));

	IF_FAILED_THROW( GtoYearsToDateInterval(
		0e0,
		&payInterval));

	mFloatRate.spotOffset.holidayFile = NULL;

	mFloatRate.rateType = GTO_SIMPLE_BASIS;

	IF_FAILED_THROW( GtoFloatRateSet(
		&matInterval,
		&payInterval,
		dayCountConv,
		spotOffsetDays,
		fixedRate,
		weight,
		&mFloatRate));

}


istream&
KRate::Get(istream& is, int drwFmt)
{

	TDateInterval	matInterval;
	TDateInterval	payInterval;
	long		dayCountConv = -1L;
	long		spotOffsetDays = 0L;
	double		spread = 0e0;
	double		weight = 0e0;

	if (drwFmt == FALSE) {
		THROW_NA
	} else {
		matInterval = getTDateInterval(is,
			"KRate::Get: matInterval.");
		payInterval = getTDateInterval(is,
			"KRate::Get: payInterval.");
		dayCountConv = getTDayCount(is,
			"KRate::Get: day count conv.");

		spotOffsetDays = 0L;

		weight = getDouble(is,
			"KRate::Get: weight.");
		spread = getDouble(is,
			"KRate::Get: spread.");

		if (mFloatRate.spotOffset.holidayFile != NULL) {
			FREE(mFloatRate.spotOffset.holidayFile);
			mFloatRate.spotOffset.holidayFile = NULL;
		}
		IF_FAILED_THROW( GtoFloatRateSet(
			&matInterval,
			&payInterval,
			dayCountConv,
			spotOffsetDays,
			spread,
			weight,
			&mFloatRate));
		mFloatRate.rateType = GTO_SIMPLE_BASIS;
	}

	return(is);
}





//---------------------------------------------------------------

ostream&
KRate::Put(ostream& os, int indent) const
{
	if (IsFloating()) {
		if (mFloatRate.dayCountConv > 0)
        {
	        os << format("RATE: %s %s %s x %8.4f + %8.4f Offset %s %-10s",
		    DrlTDateIntervalPrint(NULL, mFloatRate.matInterval),
		    DrlTDateIntervalPrint(NULL, mFloatRate.payInterval),
		    DrlTDayCountPrint(NULL, mFloatRate.dayCountConv),
		    mFloatRate.weight,
		    mFloatRate.spread,
		    DrlTDateIntervalPrint(NULL, mFloatRate.spotOffset.interval),
		    mCurveName.c_str());
        }
        else    // risky dcc
        {
	        os << format("RATE: %s %s %sD x %8.4f + %8.4f Offset %s %-10s",
		    DrlTDateIntervalPrint(NULL, mFloatRate.matInterval),
		    DrlTDateIntervalPrint(NULL, mFloatRate.payInterval),
		    DrlTDayCountPrint(NULL, labs(mFloatRate.dayCountConv)),
		    mFloatRate.weight,
		    mFloatRate.spread,
		    DrlTDateIntervalPrint(NULL, mFloatRate.spotOffset.interval),
		    mCurveName.c_str());
        }
	} else {
	    os << format("RATE: %8.4f FIXED",
		mFloatRate.spread);
	}

	return(os);
}



//---------------------------------------------------------------

ostream&
KRate::WriteSimple(ostream& os)
{
	if (IsFloating()) {
			if (mFloatRate.dayCountConv > 0)
			{

				os  << "FLOATRATE("
	    			<< format("%s, %s, %s, \"%s\"",
							DrlTDateIntervalPrint(NULL, mFloatRate.matInterval),
							DrlTDateIntervalPrint(NULL, mFloatRate.payInterval),
							DrlTDayCountPrint(NULL, mFloatRate.dayCountConv),
							mCurveName.c_str())
					<< ")";
			}
			else  // risky dcc
			{
				os  << "FLOATRATE("
	    			<< format("%s, %s, %sD, \"%s\"",
							DrlTDateIntervalPrint(NULL, mFloatRate.matInterval),
							DrlTDateIntervalPrint(NULL, mFloatRate.payInterval),
							DrlTDayCountPrint(NULL, labs(mFloatRate.dayCountConv)),
							mCurveName.c_str())
					<< ")";
		    
			}
	} else {
	        os  << "FIXEDRATE("
	    	    << format("%f", mFloatRate.spread)
		    << ")";
	}

	WriteDone();

	return(os);
}




//---------------------------------------------------------------

ostream&
KRate::YacctionWrite(ostream& os, int indent)
{

	if (GetWriteFlag())
	{
	    os  << GetName() << "=";

	    WriteSimple(os); 

	    os << ";" << endl << endl;

	    WriteDone();
	}

	return(os);
}




//---------------------------------------------------------------

void
KRate::SetSpotOffset(int numDays, int isBusDays)
{
	KDateInterval iv = KDateInterval(numDays, isBusDays);

	_TDateAdjIntvlCopy(
		&mFloatRate.spotOffset, 
		(const TDateAdjIntvl*) iv);
}


//---------------------------------------------------------------
// Convenience routine to compute a spot value of rate.
//

double
KRate::Spot(const KZCurve& zcCurve) const	// (I) zero curve
{
	double	retVal;
	TCurve	*zcCurveAlib = (TCurve*)(zcCurve);

	KRate	rateCopy = *this;	//$$$ Not very efficient, 

        // IR Rate only, no default accrual 
        rateCopy.mFloatRate.dayCountConv = labs(rateCopy.mFloatRate.dayCountConv);
	if (IsFloating()) {
		IF_FAILED_THROW( GtoForwardRate(
			zcCurveAlib,
			zcCurve.ZeroInterpType(),
			&rateCopy.mFloatRate,
			zcCurveAlib->fBaseDate,
			&retVal));
	} else {
		retVal = mFloatRate.spread;
	}

	return(retVal);

}


//---------------------------------------------------------------
// Convenience routine to compute a forward rate.
//

double
KRate::Forward(
	const KZCurve& zcCurve,	// (I) zero curve
	TDate rateEffDate)	// (I) rate effective date
		const
{
	double	retVal;
	TCurve	*zcCurveAlib = (TCurve*)(zcCurve);
	KRate	rateCopy = *this;	//$$$ Not very efficient, 

        // IR Rate only, no default accrual 
        rateCopy.mFloatRate.dayCountConv = labs(rateCopy.mFloatRate.dayCountConv);
	if (IsFloating()) {
		IF_FAILED_THROW( GtoForwardRate(
			zcCurveAlib,
			zcCurve.ZeroInterpType(),
			&(rateCopy.mFloatRate),
			rateEffDate,
			&retVal));
	} else {
		retVal = mFloatRate.spread;
	}

	return(retVal);

}



//---------------------------------------------------------------
//
//

double
KRate::Forward(
	TCurve *zcCurve,	// (I) zero curve
	TDate resetDate) const	// (I) reset date
{
	double	retVal;
	TDate	effDate;
	TCurve	*zcCurveAlib = (TCurve*)(zcCurve);
	KRate	rateCopy = *this;	//$$$ Not very efficient, 

        // IR Rate only, no default accrual 
        rateCopy.mFloatRate.dayCountConv = labs(rateCopy.mFloatRate.dayCountConv);

	if (IsFloating()) {
		// Compute effective date
		effDate = resetDate + SpotOffset();

		IF_FAILED_THROW( GtoForwardRate(
			zcCurveAlib,
			GTO_LINEAR_INTERP,
			&(rateCopy.mFloatRate),
			effDate,
			&retVal));
	} else {
		retVal = mFloatRate.spread;
	}

	return(retVal);

}



//---------------------------------------------------------------
// Convenience routine to compute a forward CDS spread.
//

double
KRate::ForwardCDS(
	const KZCurve& irZCurve,	// (I) IR  zero curve
	const KZCurve& crZCurve,	// (I) CDS zero curve
	double         recovery,	// (I) Recovery
	TDate rateEffDate)	        // (I) rate effective date
		const
{
	double	fwdSpread, annuity;
	TCurve	*irZCurveAlib = (TCurve*)(irZCurve);
	TCurve	*crZCurveAlib = (TCurve*)(crZCurve);
	KRate	rateCopy = *this;	//$$$ Not very efficient, 

        TDateInterval  intgrlFreq,
                       delayItvl;

        // Default payment delay
        IF_FAILED_THROW( GtoMakeDateInterval(
                               0,
                               'D',
                               &delayItvl));

        // Protection leg integration freq
        IF_FAILED_THROW( GtoMakeDateInterval(
                               1,
                               'M',
                               &intgrlFreq));

	if (IsFloating()) {
	    IF_FAILED_THROW( CrxFwdParCDSSpread(
                        &fwdSpread,
                        &annuity,
			MIN(irZCurve.BaseDate(), crZCurve.BaseDate()), // ?today
			MIN(irZCurve.BaseDate(), crZCurve.BaseDate()), // ?val
			rateEffDate,
			rateEffDate + Maturity(),
			(TDateInterval)(rateCopy.Frequency()),
			labs((long)(rateCopy.DayCc())),
			SHORT_FRONT,
			(rateCopy.DayCc().isNegative() ? ACCRUAL_PAY_ALL : ACCRUAL_PAY_NONE),
                        PAY_DEF,
                        delayItvl,
                        recovery,
                        intgrlFreq,	// "1M" interval
			irZCurveAlib,
			crZCurveAlib));
	} else {
	    fwdSpread = mFloatRate.spread;
	}

	return(fwdSpread);

}

//---------------------------------------------------------------


double
KRate::ForwardAdjusted(
	const KZCurve& zcDiscount,// (I) Discount zero curve
	const KZCurve& zcIndex,	// (I) Zero curve for fwdrate
	TDate rateResetDate,	// (I) Reset date for forward rate
	TDate rateEffDate,	// (I) rate effective date
	TDate ratePayDate,	// (I) Payment date (for delay adj only)
	TDate todayDate,	// (I) Where volatity starts
	double volFwdRate,	// (I) Vol of rate def. by fwdRateInfo
	double volResetToPay)	// (I) Vol of rate from reset to pay
{
static	char	routine[] = "KRate::ForwardAdjusted";
	double	retVal;

	double	duration,
		convexity,
		convexAdj,
		delayDuration,
		delayAdj,
		adjustment,

		forwardRate,
		fwdDisRate,
		payFreqD,		
		maturityYears,
		yearsToReset,

		correlation = 1e0;

	long	payFreq;
	long	daysResetToPay;


	TCurve	*zcDiscountAlib = (TCurve*)(zcDiscount);
	TCurve	*zcIndexAlib = (TCurve*)(zcIndex);



   try {
	// Compute fwd
	forwardRate = Forward(zcIndex, rateEffDate);

	//
	IF_FAILED_THROW( GtoDateIntervalToFreq(
		&mFloatRate.payInterval,
		&payFreqD));
	payFreq = (long) payFreqD;

	IF_FAILED_THROW( GtoDateIntervalToYears(
		&mFloatRate.matInterval,
		&maturityYears));

	IF_FAILED_THROW( GtoDayCountFraction(
		todayDate,
		rateResetDate,
		GTO_ACT_365F,
		&yearsToReset));


	daysResetToPay = ratePayDate - rateResetDate;





	//
	// Convexity adjustment for forwardRate
	//
	if (payFreq ISNT 0)
	{
		IF_FAILED_THROW( GtoBondModDuration(
			forwardRate,
			forwardRate,
			payFreq,
			maturityYears,
			GTO_STUB_SIMPLE,
			&duration));

		IF_FAILED_THROW( GtoBondConvexity(
			forwardRate,
			forwardRate,
			payFreq,
			maturityYears,
			GTO_STUB_SIMPLE,
			&convexity));

	} else { // Frequency = 0 implies zero-coupon rate */
		long dccDenom;          /* Either 365 or 360 */

		IF_FAILED_THROW( GtoDayCountDenom(
			mFloatRate.dayCountConv,
			&dccDenom));

		IF_FAILED_THROW( GtoMoneyMarketModDuration(
			forwardRate, 
			(long)(maturityYears*dccDenom),
			dccDenom,
			&duration));

		IF_FAILED_THROW( GtoMoneyMarketConvexity(
			forwardRate, 
			(long)(maturityYears*dccDenom),
			dccDenom,
			&convexity));
	}


	IF_FAILED_THROW( GtoSwapConvexityAdj(
		forwardRate,
		duration,
		convexity,
		volFwdRate,
		yearsToReset,
		&convexAdj));

	//
	// delay adjustment for forwardRate
	//

	// Get Act/365F discount rate from reset date to pay date.
	if (rateResetDate < ratePayDate) {
		IF_FAILED_THROW( GtoForwardFromZCurve(
			zcDiscountAlib,
			zcDiscount.ZeroInterpType(),
			rateResetDate,
			ratePayDate,
			GTO_ACT_365F,
			GTO_SIMPLE_BASIS,
			&fwdDisRate));
	} else {
		fwdDisRate = 0e0;	// resetDate = paymentDate
	}



	IF_FAILED_THROW( GtoMoneyMarketModDuration(
		fwdDisRate,
		daysResetToPay,
		365, /* Since rate is Act/365F */
		&delayDuration));

	IF_FAILED_THROW( GtoSwapDelayAdj(
		forwardRate,
		volFwdRate,
		yearsToReset,
		fwdDisRate,
		volResetToPay,
		delayDuration,
		correlation,
		&delayAdj));


	//
	// Total adjustment for forwardRate
	//
	adjustment = convexAdj + delayAdj;

	retVal = forwardRate + adjustment;

	return(retVal);

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}


//***************************************************************
//
//	KRateReset functions.
//
//***************************************************************


//---------------------------------------------------------------

KRateReset::KRateReset(
	TDate resetDate,		// (I) observation date
	const KRate& rate)		// (I) rate description
{
	mRate = rate;
	mResetDate = resetDate;
	mEffDate = mResetDate + rate.SpotOffset();
}



//---------------------------------------------------------------

KRateReset::KRateReset(
	TDate resetDate,		// (I) observation date
	TDate effDate,			// (I) effective date
	const KRate& rate)		// (I) rate description
{
	mRate = rate;
	mResetDate = resetDate;
	mEffDate = effDate;

	// If default of spotOffset is zero, set it to be 
	// the number of days between effDate and resetDate.
	if (IS_ALMOST_ZERO(rate.SpotOffset().Years()))
	{
		mRate.SetSpotOffset(effDate - resetDate, 
			      GTO_DATE_ADJ_TYPE_CALENDAR);
	}
}


//---------------------------------------------------------------

KRateReset::KRateReset(
	const String &curveName,	// (I) curve name
	TDate resetDate,		// (I) observation date
	TDate effDate,			// (I) effective date
	const KRate& rate)		// (I) rate description
{
	mRate = rate;

	mRate.SetCurveName(curveName);
	mResetDate = resetDate;
	mEffDate = effDate;

	// If default of spotOffset is zero, set it to be 
	// the number of days between effDate and resetDate.
	if (IS_ALMOST_ZERO(rate.SpotOffset().Years()))
	{
		mRate.SetSpotOffset(effDate - resetDate, 
			      GTO_DATE_ADJ_TYPE_CALENDAR);
	}
}


//---------------------------------------------------------------
//


KRateReset::KRateReset(
	TDate resetDate,		// (I) observation date
	TDate effDate,			// (I) effective date
	TDate endDate,			// (I) end date
	const KRate& rate)		// (I) rate description
{
	mRate = KRate(rate.CurveName(),
		effDate,
		endDate,
		rate.DayCc(),
		rate.SpotOffset(),
		rate.Spread(),
		1e0);

	mResetDate = resetDate;
	mEffDate = effDate;

	// If default of spotOffset is zero, set it to be 
	// the number of days between effDate and resetDate.
	if (IS_ALMOST_ZERO(rate.SpotOffset().Years()))
	{
		mRate.SetSpotOffset(effDate - resetDate, 
			      GTO_DATE_ADJ_TYPE_CALENDAR);
	}
}



//---------------------------------------------------------------
// Copy constructor

KRateReset::KRateReset(const KRateReset& rateReset)
{
	mRate = rateReset.mRate;
	mResetDate = rateReset.mResetDate;
	mEffDate = rateReset.mEffDate;
}



//---------------------------------------------------------------
bool
KRateReset::operator<(const KRateReset& rate) const
{
	if (!(mResetDate < rate.mResetDate))
		return (FALSE);
	if (!(mRate< rate.mRate))
		return (FALSE);
	return (TRUE);
}


//---------------------------------------------------------------
// Generate a series of zero reset dates from KRateReset
void
KRateReset::GetZeroResetList(KVector(KZeroReset) &zeroResetList) const
{
static  char    routine[] = "KRateReset::GetZeroResetListFromRateReset";

	TDate	resetDate = mResetDate;
	TDate   matDate, cfDate;
	String	curveName = mRate.CurveName();
	KDateInterval	freqIntvl = mRate.Frequency();
	KDateInterval	matIntvl  = mRate.Maturity();
	int	numPrd;

	matDate = mEffDate + matIntvl;

	// First cash flow date is the effective date
	numPrd = 0;
	cfDate = mEffDate;
	
	while(cfDate < matDate) {
		zeroResetList.push_back(KZeroReset(curveName,		
						   resetDate,
					 	   cfDate));
		numPrd++;
		cfDate = mEffDate + numPrd*freqIntvl;	
	}	

	// Check the last cash flow date is the same as 
	// the maturity date, and add to the list
	if(cfDate == matDate)
		zeroResetList.push_back(KZeroReset(curveName,		
						   resetDate,
					 	   cfDate));
	else
		throw KFailure("%s: the KRate maturity (%s) is NOT "
				"multiples of frequency interval (%s).\n",
				routine,
				matIntvl.Str(),
				freqIntvl.Str());
				
}


//---------------------------------------------------------------

ostream&
KRateReset::Put(ostream& os, int indent) const
{
	os << format("RATERESET: %10s %10s %s ",
		DrlTDatePrint(NULL, mResetDate),
		DrlTDatePrint(NULL, mEffDate),
		mRate.CurveName().c_str());
	os << this->Rate();
	return (os);
}



//***************************************************************
//
//	KZeroReset functions.
//
//***************************************************************


//---------------------------------------------------------------
// Copy constructor

KZeroReset::KZeroReset(const KZeroReset& zeroReset)
{
	mCurveName = zeroReset.mCurveName;
	mMaturityDate = zeroReset.mMaturityDate;
	mEarliestDate = zeroReset.mEarliestDate;
}



//--------------------------------------------------------------
// Check date validity
//
void
KZeroReset::IsValid() const
{
static  char    routine[] = "KZeroReset::IsValid";
 
 try{
	// Check the validity of curve name
	// to be done
 
	if ((mMaturityDate <= 1L) || (mMaturityDate >= 500000L))
		throw KFailure("%s: invalid maturity date %ld.\n",
				routine, mMaturityDate);
 
	if ((mEarliestDate <= 1L) || (mEarliestDate >= 500000L))
		throw KFailure("%s: invalid early usage date %ld.\n",
				routine, mEarliestDate);
 
	if (mMaturityDate < mEarliestDate)
		throw KFailure("%s: mMaturityDate (%s) < mEarliestDate (%s).\n",
				routine, GtoFormatDate(mMaturityDate),
				GtoFormatDate(mEarliestDate));
    }
    catch(KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
 
}



//---------------------------------------------------------------

ostream&
KZeroReset::Put(ostream& os, int indent) const
{
	os << format("ZERORESET: %10s %10s %s ",
		DrlTDatePrint(NULL, mEarliestDate),
		DrlTDatePrint(NULL, mMaturityDate),
		mCurveName.c_str());
	return (os);
}

