/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      David Liu
 ************************************************************************/
#include "kswap.h"

#include "kutilios.h"

extern "C" {
#include "datelist.h"           // GtoNewDateList
#include "zr2coup.h"		// GtoZerosToCouponsPoint
#include "zr2simp.h"		// GtoZerosToSimplePoint
#include "interp.h"	
#include "swopblk.h"
#include "drltime.h"
};

#define DEBUG_CURVE
#undef	DEBUG_CURVE

//===============================================================
//
//===============================================================


//---------------------------------------------------------------
//
KBMSwap::KBMSwap()
{
	mZC = NULL;
}


KBMSwap::~KBMSwap()
{
	mZC = NULL;
}



//---------------------------------------------------------------
//
KBMSwap::KBMSwap(
	const char *name,		// (I) name
	const TDate today,		// (I) today's date
	const TDate startDate,          // (I) start dates
	const TDate matDate,		// (I) maturity dates
	const KDateInterval &freq,	// (I) frequency as interval
	const double    vol,		// (I) volatility
	KDayCc &dayCc,			// (I) pay day count
	TBoolean stubAtEnd,             // (I) T=stub at end, F=at beg
	const String &curveName,	// (I) curve name
	const KZCurve	&zc)		// (I) Zero curve
{
static	char	routine[] = "KBMSwap::KBMSwap";

	TDateList	*datesDL = (TDateList *)NULL;
	int		idx;

    try {

	mName	  = name;
	mTodayDate = today;
	mStDate = startDate;
	mFreq	  = freq;
	mVol      = vol;
	mDayCc    = dayCc;
	mStubAtEnd = stubAtEnd;
	mCurveName = curveName;
	mZC	  = &zc;

	// Set up payment dates. The first reset index is the next
	// to last payment date (since we're moving backwards.)
	datesDL = GtoNewDateList(
			startDate,
			matDate,
			&((TDateInterval) freq),
			mStubAtEnd);
	ASSERT_OR_THROW(datesDL != NULL);


	// Assume payment date = next reset date
	for (idx=0; idx < datesDL->fNumItems-1; idx++)
		mPayDates.push_back(datesDL->fArray[idx+1]);


	GtoFreeDateList(datesDL);

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}

	

//--------------------------------------------------------------
//
void
KBMSwap::Initialize()
{
static	char	routine[] = "KBMSwap::Initialize";

	double	parYield;
	double	dcfrn;

    try {
	

	// Compute the coupon payments.
	//

	// Par yield
	parYield = ParYield();

	// The first payment
	KVector(TDate)::iterator it = mPayDates.begin();

	// Day count fraction
	ASSERT_OR_THROW(GtoDayCountFraction(
				StartDate(),
				*it,
				GTO_ACT_365,
                                &dcfrn) == SUCCESS);


	mCoupons.push_back(parYield*dcfrn);

	
	// Move on to next payment date
	++it;
	for (; it != mPayDates.end(); ++it)
	{
		ASSERT_OR_THROW(GtoDayCountFraction(
					*(it-1),
					*it,
					GTO_ACT_365,
                                        &dcfrn) == SUCCESS);

		mCoupons.push_back(parYield*dcfrn);
	}

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}




//--------------------------------------------------------------
//
double
KBMSwap::ParYield()
{
static	char	routine[] = "KBMSwap::ParYield";

	double	parYield;

    try {

	if(mPayDates.size() > 1)
	{
		IF_FAILED_THROW(GtoZerosToCouponsPoint(
				(TCurve*)(*mZC),
				(*mZC).ZeroInterpType(),
				mStDate,
				&(TDateInterval)mFreq,
				mPayDates.back(),
				(long)mDayCc,
				GTO_STUB_NONE,
				mStubAtEnd,
				&parYield));
	}
	else
	{
		IF_FAILED_THROW(GtoZerosToSimplePoint(
				(TCurve*)(*mZC),
				(*mZC).ZeroInterpType(),
				mStDate,
				mPayDates.back(),
				(long)mDayCc,
				&parYield));
	}

	return (parYield);

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}



	
//--------------------------------------------------------------
//
double
KBMSwap::BSPrice()
{
static	char	routine[] = "KBMSwap::BSPrice";

	double	parYield;

	double	price;
	TOptionPropertiesIR result;

    try {
	// Compute the par yield
	//
	IF_FAILED_THROW(GtoZerosToCouponsPoint(
				(TCurve*)(*mZC),
				(*mZC).ZeroInterpType(),
				mStDate,
				&(TDateInterval)mFreq,
				mPayDates.back(),
				(long)mDayCc,
				GTO_STUB_NONE,
				mStubAtEnd,
				&parYield));

	// Compute the B-S price
	//
	IF_FAILED_THROW(GtoSwaptionPVBlackVanilla(
				(TCurve*)(*mZC),
				(*mZC).ZeroInterpType(),
				mStDate,
				&(TDateInterval)mFreq,
				mPayDates.back(),
				(long)mDayCc,
				GTO_STUB_NONE,
				mStubAtEnd,
				parYield,
				mStDate,
				mTodayDate,
				mVol,
				TRUE,	// true=call, false=put
				TRUE,	// true=ytm, false=zcurve
				GTO_OPTION_PRICE,
				&result));
	
	price = result.fPrice;

	return (price);
				
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}



//--------------------------------------------------------------
//
double
KBMSwap::BSVega()
{
static	char	routine[] = "KBMSwap::BSVEga";

	double	parYield;
	double	vega;

	TOptionPropertiesIR result;

    try {
	// Compute the par yield
	//
	IF_FAILED_THROW(GtoZerosToCouponsPoint(
				(TCurve*)(*mZC),
				(*mZC).ZeroInterpType(),
				mStDate,
				&(TDateInterval)mFreq,
				mPayDates.back(),
				(long)mDayCc,
				GTO_STUB_NONE,
				mStubAtEnd,
				&parYield));

	// Compute the B-S price
	//
	IF_FAILED_THROW(GtoSwaptionPVBlackVanilla(
				(TCurve*)(*mZC),
				(*mZC).ZeroInterpType(),
				mStDate,
				&(TDateInterval)mFreq,
				mPayDates.back(),
				(long)mDayCc,
				GTO_STUB_NONE,
				mStubAtEnd,
				parYield,
				mStDate,
				mTodayDate,
				mVol,
				TRUE,	// true=call, false=put
				TRUE,	// true=ytm, false=zcurve
				GTO_OPTION_VEGA,
				&result));
	
	vega = result.fVega;

	return (vega);
				
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}




//---------------------------------------------------------------
//
void
KBMSwap::GetZeroResetList(KVector(KZeroReset) &zeroResetList)
{
static	char	routine[] = "KBMSwap::GetZeroResetList";


    try {

	for (KVector(TDate)::iterator it = mPayDates.begin();
		it != mPayDates.end(); ++it)
	{
		zeroResetList.push_back(KZeroReset(mCurveName,
                                                   mStDate,
                                                   *it));
	}

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}



//---------------------------------------------------------------
// 
ostream&
KBMSwap::Put(ostream& os) const
{

	int	idx;

	os << "NAME: `" << GetName() << "'" << endl;
	os << "Swap Curve: " << mCurveName << endl;
	os << "Swap Volatility: " << mVol << endl;
	os << "Swap Frequency: " << mFreq << endl;
	os << "Start Date: " << GtoFormatDate(mStDate) << endl;
	os << "Pay Dates: " << endl;

	for (idx=0; idx <= mPayDates.size()-1; idx++)
		os << "\t" << GtoFormatDate(mPayDates[idx]) << endl;
		
	return(os);

}
