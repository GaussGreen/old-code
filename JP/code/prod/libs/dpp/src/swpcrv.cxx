/****************************************************************
 * Module:	PenGuin
 * Submodule:	
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#define	_dpp_SRC
#include <ctype.h>
#include "kswpcrv.h"		// class header file

#include "kutilios.h"

extern	"C" {
#include "tcurve.h"
#include "gtonpi.h"
#include "duration.h"
#include "date_sup.h"
#include "ldate.h"
#include "cerror.h"
#include "convert.h"
#include "yearfrac.h"
#include "zr2coup.h"		/* Analytics C Library */
#include "zr2simp.h"		/* Analytics C Library */


#include "drlvtype.h"
#include "drlstr.h"
#include "drlio.h"
#include "drlinter.h"
#include "drlmem.h"

#include "drltime.h"
#include "drllineq.h"		/* DrlRealLinSysSvd */
};


//---------------------------------------------------------------
// Default constructor.

KPCurve::KPCurve()
	: KZCurve()
{
	mFormat = KPCurve::STD;

	mZcCurve = NULL;

	mRateType = NULL;
	mEffDur = NULL;
	mMmDen = 0;

}

//---------------------------------------------------------------
// Copy constructor.


KPCurve::KPCurve(const KPCurve& swpcrv)
	: KZCurve()
{
	int	i;

	create(
		swpcrv.BaseDate(),
		swpcrv.NumItems(),
		swpcrv.Basis(),
		swpcrv.DayCc());

	mMmDen = swpcrv.mMmDen;

	for (i=0; i<=NumItems()-1; i++) {
	    Date(i) = swpcrv.Date(i);
	    Rate(i) = swpcrv.Rate(i);
	    mRateType[i] = swpcrv.mRateType[i];
	    mEffDur[i] = swpcrv.mEffDur[i];
	}
}



//---------------------------------------------------------------
// Copy constructor.


KPCurve::KPCurve(const TCurve* zcCurve)
	: KZCurve()
{
	int	i;

	create(
		zcCurve->fBaseDate,
		zcCurve->fNumItems,
		zcCurve->fBasis,
		zcCurve->fDayCountConv);

	for (i=0; i<=NumItems()-1; i++) {
	    Date(i) = zcCurve->fArray[i].fDate;
	    Rate(i) = zcCurve->fArray[i].fRate;
	    mRateType[i] = 'Z';
	    mEffDur[i] = 0e0;
	}
}





//---------------------------------------------------------------
// Destructor.

KPCurve::~KPCurve()
{
	destroy();
}

//---------------------------------------------------------------
// private memory management

void
KPCurve::create(
	TDate baseDate,
	int size,
	double basis,
	TDayCount dayCountConv)
{
	mZcCurve = GtoNewTCurve(baseDate, size, basis, dayCountConv);
        if (!mZcCurve) throw KFailure("GtoNewTCurve failed.\n");

	mRateType = new char[size];
        if (!mRateType) throw KFailure("Memory failure\n");
	mEffDur   = new double[size];
        if (!mEffDur) throw KFailure("Memory failure\n");

}



void
KPCurve::destroy()
{
	GtoFreeTCurve(mZcCurve);
	delete mRateType;
	delete mEffDur;

	mFormat = 0;
	mZcCurve = NULL;
	mRateType = NULL;
	mEffDur = NULL;
	mMmDen = 0;
}



//---------------------------------------------------------------
// Assignment operator.

KPCurve&
KPCurve::operator=(const KPCurve& swpcrv)
{
	int	i;

	destroy();
	// empty
	if (swpcrv.isempty()) return(*this);

	create(
		swpcrv.BaseDate(),
		swpcrv.NumItems(),
		swpcrv.Basis(),
		swpcrv.DayCc());

	mMmDen = swpcrv.mMmDen;

	for (i=0; i<=NumItems()-1; i++) {
	    Date(i) = swpcrv.Date(i);
	    Rate(i) = swpcrv.Rate(i);
	    mRateType[i] = swpcrv.mRateType[i];
	    mEffDur[i] = swpcrv.mEffDur[i];
	}

	return(*this);
}




//---------------------------------------------------------------
// Assignment operator.

KPCurve&
KPCurve::operator=(const TCurve* zcCurve)
{
	int	i;

	destroy();
	if (!zcCurve) return(*this);
	create(
		zcCurve->fBaseDate,
		zcCurve->fNumItems,
		zcCurve->fBasis,
		zcCurve->fDayCountConv);

	mMmDen = 0;

	for (i=0; i<=NumItems()-1; i++) {
	    Date(i) = zcCurve->fArray[i].fDate;
	    Rate(i) = zcCurve->fArray[i].fRate;
	    mRateType[i] = 'Z';
	    mEffDur[i] = 0e0;
	}

	return(*this);
}


//---------------------------------------------------------------
// Read class from a stream.

istream&
KPCurve::Get(istream& is)
{
static	char	routine[] = "KPCurve::Get";

	int		i;
	TDate		xBaseDate;		// tmp storage
	double		xBasis;
	TDayCount	xDayCountConv;
	int		xMmDen;
	int		xNumItems;


    try {

	xBaseDate = getTDate(is, "BaseDate");
	xBasis = getDouble(is, "Basis");
	xDayCountConv = getTDayCount(is, "DayCountConv");
	xMmDen = getInt(is, "mMmDen" );
	xNumItems = getInt(is, "NumItems");

	destroy();
	create(
		xBaseDate,
		xNumItems,
		xBasis,
		xDayCountConv);
	mMmDen = xMmDen;

	for (i=0; i<=NumItems()-1; i++) {
		Date(i) = getTDate(is, "mDate");
		getString(is, "mIntervalStr");
		mRateType[i] = getChar(is, "mRateType");
		Rate(i) = getDouble(is, "mRate");
		mEffDur[i] = getDouble(is, "mEffDur");
	}

	return(is);
    }
    catch (KFailure) {
	    throw KFailure("%s: failed.\n", routine);
    }
}



//---------------------------------------------------------------
// Write class to a stream.

ostream&
KPCurve::Put(ostream& os, int oneLine) const
{
	int	i;

	if (isempty()) {
		os << "EMPTY\n";
		return(os);
	}

	switch (mFormat) {
	case KPCurve::BVDRWRAP:
		os << "# Base Volatility Frequency (A,S,Q,M)\n";
		switch (Freq()) {
			case 1:  os << "A\n"; break;
			case 2:  os << "S\n"; break;
			case 4:  os << "Q\n"; break;
			case 12: os << "M\n"; break;
			default: os << "ERROR\n"; break;
		}

		os << "# No of Base Volatility Points\n %d\n";
		os << NumItems() << '\n';

		os << "# Base Volatility Dates and Volatilities in %%\n";
		for (i=0; i<=NumItems()-1; i++) {
		    os << format(" %s    %15.10f\n",
			DrlTDatePrintYMD(NULL,
			    mZcCurve->fArray[i].fDate),
			mZcCurve->fArray[i].fRate*1e2);
		}

	    break;
	case KPCurve::LIVERATE:
	    	for (i=0; i<=NumItems()-1; i++) {
			TDateInterval	intvl;

		    IF_FAILED_THROW( GtoDateSubtract(
				Date(i),
	    			BaseDate(),
				&intvl));

		    os << format("\"%s%c\"\t%12.8f%%", 
				DrlTDateIntervalPrint(NULL, intvl),
				mRateType[i],
				Rate(i)*1e2);
		    os << "\t0.000000";
		    os << "\t0.000000";
		    os << "\t0.000000";
		    os << "\n";
	    	}
		break;



	case KPCurve::DELTA:
	    os << "KPCurve::DELTA:" << endl;
	    os << "BASE_DATE " << DrlTDatePrint(NULL, BaseDate()) << endl;

	    os << "NUM_ITEMS: " << NumItems() << endl;

	    os << " DATE             MAT   T       FACE           "
			"  DUR          1BPTWK" << endl;

	    for (i=0; i<=NumItems()-1; i++) {
		TDateInterval	intvl;

		IF_FAILED_THROW( GtoDateSubtract(
			Date(i),
			BaseDate(),
			&intvl));

		os << format(" %10s  %8s  %c  %20s  %8.4f  %15s\n", 
			DrlTDatePrint(NULL, Date(i)),
			DrlTDateIntervalPrint(NULL, intvl),
			mRateType[i],
			DrlCurPrint(NULL, Rate(i), 2),
			mEffDur[i],
			DrlCurPrint(NULL, Rate(i)*mEffDur[i] *1e-4, 2));
	    }
	    break;


	case KPCurve::BP:
	/*
	    os << "#BaseDate\n" << GtoFormatDate(BaseDate()) << '\n';
	    os << "#Basis\n" << Basis() << '\n';
	    os << "#DayCc\n"
		<< GtoFormatDayCountConv(DayCc()) << '\n';
	    os << "#mMmDen\n" << mMmDen << '\n';

	    os << "#NumItems\n" << NumItems() << '\n';
	*/

	    for (i=0; i<=NumItems()-1; i++) {
		TDateInterval	intvl;

		IF_FAILED_THROW( GtoDateSubtract(
			Date(i),
			BaseDate(),
			&intvl));

		os << format(" %10s  %8s  %c  %10.6fe-4   0\n", 
			DrlTDatePrint(NULL, Date(i)),
			DrlTDateIntervalPrint(NULL, intvl),
			mRateType[i],
			Rate(i)*1e4);
	    }
	    break;


	case KPCurve::STD:
	default:
	    os << "#BaseDate\n" << GtoFormatDate(BaseDate()) << '\n';
	    os << "#Basis\n" << Basis() << '\n';
	    os << "#DayCc\n"
		<< GtoFormatDayCountConv(DayCc()) << '\n';
	    os << "#mMmDen\n" << mMmDen << '\n';

	    os << "#NumItems\n" << NumItems() << '\n';

	    for (i=0; i<=NumItems()-1; i++) {
		TDateInterval	intvl;

		IF_FAILED_THROW( GtoDateSubtract(
			Date(i),
			BaseDate(),
			&intvl));

		os << format(" %10s  %8s  %c  %12.10f  %10.6f\n", 
			DrlTDatePrint(NULL, Date(i)),
			DrlTDateIntervalPrint(NULL, intvl),
			mRateType[i],
			Rate(i),
			mEffDur[i]);
	    }
	    break;
	}


	return(os.flush());
}





//---------------------------------------------------------------


void
KZCurveExportToWrapper(
	const char *fnam,
	const KPCurve &parCurve,
	const KZCurve &zcCurve)
{
static	char	routine[] = "KZCurveExportToWrapper";
	FILE    *fp = NULL;
	int	idx;

    try {

	ASSERT_OR_THROW((fp = fopen(fnam, "w")) != NULL);

	DrlFPrintf(fp, "# Start date\n%s\n", DrlTDatePrintYMD(NULL,
		zcCurve.BaseDate()));


	if (!parCurve.IsEmpty()) {

		DrlFPrintf(fp, "# Money Market basis (360 or 365)\n %d\n",
			parCurve.MmDen());

		DrlFPrintf(fp, "# Annual or semi-annual curve (A, S)\n");
		switch (parCurve.Freq()) {
		case 12: DrlFPrintf(fp, "M\n"); break;
		case  4: DrlFPrintf(fp, "Q\n"); break;
		case  2: DrlFPrintf(fp, "S\n"); break;
		case  1: DrlFPrintf(fp, "A\n"); break;
		default: DrlFPrintf(fp, "ERROR(%d)\n", parCurve.Freq()); break;
		}

//$$$WARNING
		DrlFPrintf(fp, "# Year basis for swaps (ACT, 365, 360)\n %s\n",
			"ACT");
	} else {
		//$$$$
		//$$$$ WARNING !!! USD SPECIFIC
		//$$$$

		DrlFPrintf(fp, "# Money Market basis (360 or 365)\n360\n");
		DrlFPrintf(fp, "# Annual or semi-annual curve (A, S)\nS\n");
		DrlFPrintf(fp, "# Year basis for swaps (ACT, 365, 360)\nACT\n");
	}




	DrlFPrintf(fp, "# No of entries\n %d\n",
                        zcCurve.NumItems());

	DrlFPrintf(fp, "# Zero Dates & Rates in %% (ACT/365F A)\n");
	for (idx=0; idx<zcCurve.NumItems(); idx++) {
	DrlFPrintf(fp, " %s    %15.10f\n",
		DrlTDatePrintYMD(NULL, zcCurve.Date(idx)),
		zcCurve.Rate(idx)*1e2);
	}
	fclose(fp);

    }
    catch (KFailure) {
	if (fp) fclose(fp);
	throw KFailure("%s: failed\n", routine);
    }
}



//---------------------------------------------------------------
//

int
CheckSameType(const KPCurve& crv1, const KPCurve& crv2)
{
	int	i;
	TDateInterval	intvl1, intvl2;
    try {
	ASSERT_OR_THROW(
		( crv1.isempty() &&  crv2.isempty()) ||
		(!crv1.isempty() && !crv2.isempty()));

	ASSERT_OR_THROW(crv1.NumItems() == crv2.NumItems());

	/*
	ASSERT_OR_THROW(crv1.BaseDate() == crv2.BaseDate());
	for (i=0; i<=crv1.NumItems()-1; i++) {
		ASSERT_OR_THROW(crv1.Date(i) == crv2.Date(i));
	}
	*/
	for (i=0; i<=crv1.NumItems()-1; i++) {
		IF_FAILED_THROW( GtoDateSubtract(
			crv1.Date(i), crv1.BaseDate(), &intvl1));
		IF_FAILED_THROW( GtoDateSubtract(
			crv2.Date(i), crv2.BaseDate(), &intvl2));
		if (!GtoIntervalEqualsInterval(&intvl1, &intvl2))
			throw KFailure();
	}




	ASSERT_OR_THROW(crv1.Freq() == crv2.Freq());

	return(TRUE);
    }
    catch (KFailure) {
	throw KFailure("CheckSameType(KPCurve&,KPCurve&): inconsistent.\n");
    }
}



//---------------------------------------------------------------
// Generates a zero curve from a swap curve using NPI methodology.
// 

TCurve*
ZCGenNPI(
	KPCurve& parCrv)	// (I) par swap curve 
{
static	char	routine[] = "KPCurveNPIZeroGen";
	int	status = FAILURE;
	int	i;

	double	rateArr[DEF_SWAPRATE_MAX];
	TDate	dateArr[DEF_SWAPRATE_MAX];
	double	priceArr[DEF_SWAPRATE_MAX];
	char	nameArr[DEF_SWAPRATE_MAX];
	int	useDefault = GtoSETINSTR;
	int	numRate;
	double	*futArr = NULL;
	TDate 	*futDateArr = NULL;
	int	numFuture = 0;
	double	yieldVolatility = 0;
	int	stubMethod = GtoSTUB3M;
	double	*stubData = NULL;
	int	stubDataSize = 0;
	int	couponInterpMethod = GtoLINEARINTERP;
	int	zeroInterpMethod = parCrv.ZeroInterpType();
	char  	fwdLength = 'Q';

	TCurve	*zcCurve;	// zero curve 


	numRate = parCrv.NumItems();
	for (i=0; i<= parCrv.NumItems()-1; i++) {
		rateArr[i] = parCrv.Rate(i);
		dateArr[i] = parCrv.Date(i);
		priceArr[i] = 1.0;
		nameArr[i] = (toupper(parCrv.mRateType[i]) == 'M' ?
				GtoMONEYNAME :
				GtoSWAPNAME);
	}


	// generate zero curve
	zcCurve = GtoNPiZC(
		parCrv.BaseDate(),
		parCrv.mMmDen,
		(int) parCrv.Basis(),	// swap freq
		parCrv.DayCc(),
		rateArr,
		dateArr,
		priceArr,
		nameArr,
		useDefault,
		numRate,
		futArr,
		futDateArr,
		numFuture,
		0, NULL,
		yieldVolatility,
		stubMethod,
		stubData,
		stubDataSize,
		couponInterpMethod,
		zeroInterpMethod,
		fwdLength,
		0
		);

	return (zcCurve);
}


//---------------------------------------------------------------
// Generates a zero curve from a swap curve using NPI methodology.
// 

void
ZCGenNPI(
	KZCurve& zcCrv,		// (O) zero curve 
	KPCurve& parCrv)	// (I) par swap curve 
{
	TCurve	*zcCurve;

	if (parCrv.isempty())
		throw KFailure("ZCGenNPI: empty input curve.\n");

	zcCurve = ZCGenNPI(parCrv);
	if (zcCurve == NULL) throw KFailure();

	zcCrv = zcCurve;
	GtoFreeTCurve(zcCurve);
	return;
}


//---------------------------------------------------------------
//

KPCurve&
KPCurve::operator=(double argument)
{
	int	i;
	if (isempty()) return (*this);
	for (i=0; i<= NumItems()-1; i++)
		Rate(i) = argument;
	return(*this);
}

//---------------------------------------------------------------
//

KPCurve&
KPCurve::operator+=(double argument)
{
	int	i;
	if (isempty()) return (*this);
	for (i=0; i<= NumItems()-1; i++)
		Rate(i) += argument;
	return(*this);
}

//---------------------------------------------------------------
//

KPCurve&
KPCurve::operator*=(double argument)
{
	int	i;
	if (isempty()) return (*this);
	for (i=0; i<= NumItems()-1; i++)
		Rate(i) *= argument;
	return(*this);
}


//---------------------------------------------------------------
//

KPCurve&
KPCurve::operator+=(const KPCurve& crv)
{
	int	i;

	if (isempty()) {
		*this = crv;
	} else {
		CheckSameType(*this, crv);
		for (i=0; i<= NumItems()-1; i++) 
			Rate(i) += crv.Rate(i);
	}
	return(*this);
}


//---------------------------------------------------------------
// Recalculates all effecive durations using a zero curve {\tt zcCurve}.
//

KPCurve&
KPCurve::RecalcDuration(int upDown)
{
static	char	routine[] = "KPCurve::RecalcDuration";
	int	status = FAILURE;
	int	idx;

	TDate	maturityDate;
	double	yieldSave, yield, effDur, zeroRate,
		rateTwkSize = 1e-4;
	TCurve	*zcCurve = NULL;

    try {

	for (idx=0; idx<= NumItems()-1; idx++) {
	    //
	    // Up Tweak
	    //

	    yieldSave = Rate(idx);
	    Rate(idx) = yieldSave + rateTwkSize;

	    // Generate Zero Curve
	    zcCurve = ZCGenNPI(*this);

	    yield = Rate(idx);
	    maturityDate = Date(idx);

	    switch (toupper(mRateType[idx])) {
	    case 'M':
		IF_FAILED_THROW(GtoMoneyMarketModDuration(
			yield,
			maturityDate-BaseDate(),
			mMmDen,
			&mEffDur[idx]));
		break;
	    case 'S':
		IF_FAILED_THROW(GtoInterpRate(
			maturityDate,
			zcCurve,
			ZeroInterpType(),
			&zeroRate));

		IF_FAILED_THROW(GtoBondEffDuration(
			yield,
			zcCurve->fBasis,
			(maturityDate-BaseDate())/365e0,
			zeroRate,
			0e0,
			&mEffDur[idx]));
			
		break;
	    default:
		throw KFailure("%s: bad rate type.\n", routine);
	    }
	    GtoFreeTCurve(zcCurve); zcCurve = NULL;

	    //
	    // Down Tweak
	    //
	    if (upDown != FALSE) {
	    	Rate(idx) = yieldSave - rateTwkSize;

	    	// Generate Zero Curve 
	        zcCurve = ZCGenNPI(*this);

	    	yield = Rate(idx);
	    	maturityDate = Date(idx);

	    	switch (toupper(mRateType[idx])) {
	    	case 'M':
			IF_FAILED_THROW(GtoMoneyMarketModDuration(
				yield,
				maturityDate-BaseDate(),
				mMmDen,
				&effDur));
			break;
	    	case 'S':
			IF_FAILED_THROW(GtoInterpRate(
				maturityDate,
				zcCurve,
				ZeroInterpType(),
				&zeroRate));
	
			IF_FAILED_THROW(GtoBondEffDuration(
				yield,
				zcCurve->fBasis,
				(maturityDate-BaseDate())/365e0,
				zeroRate,
				0e0,
				&effDur));
			
			break;
	    	default:
			throw KFailure("%s: bad rate type.\n", routine);
	    	}
	    	GtoFreeTCurve(zcCurve); zcCurve = NULL;

		mEffDur[idx] = 0.5e0 * (mEffDur[idx] + effDur);

	    }


	    /* Restore old yield */
	    Rate(idx) = yieldSave;
	}
	return (*this);
    }
    catch (KFailure) {
	GtoFreeTCurve(zcCurve);
	throw KFailure("%s: failed\n", routine);
    }
}

//---------------------------------------------------------------

KPCurve&
KPCurve::ShiftDate(TDate baseDate, int forward)
{

	if (forward) {
		//
		// Forward the environment
		//
		KZCurve	zcCurve;

		// Compute zc curve
		ZCGenNPI(zcCurve, *this);

		// fwd zc curve
		zcCurve.ShiftDate(baseDate, TRUE);

		// Shift par curve dates
		long	diff = baseDate - BaseDate();
		int	idx;

		BaseDate() += diff;
		for (idx=0; idx<= NumItems()-1; idx++) {
		    	Date(idx) += diff;
		}

		// Recalc yields
		RecalcYields((TCurve*) zcCurve);


	} else {

		long	diff = baseDate - BaseDate();
		int	idx;

		BaseDate() += diff;
		for (idx=0; idx<= NumItems()-1; idx++) {
		    	Date(idx) += diff;
		}
	}

	return(*this);


}


//---------------------------------------------------------------
// Shifts all rates by {\tt value}.

void
KPCurve::ShiftRates(double value)
{
static	char	routine[] = "KPCurve::ShiftRates";
	int	idx;

	for (idx=0; idx<= NumItems()-1; idx++) {
	    Rate(idx) += value;
	}

}



//---------------------------------------------------------------
// Recalculates all dates and par yields
// from a zero curve (the value date is changed).

void
KPCurve::RecalcYields(TCurve *zcCurve)
{
static	char	routine[] = "KPCurve::RecalYield";

	int		idx;
	TDate		oldValueDate = BaseDate();
	TDateInterval	matInterval,
			payInterval;
	TDayCount	mmDayCount;


    try {

	oldValueDate = BaseDate();
	BaseDate() = zcCurve->fBaseDate;

	// convert frequency to interval
	IF_FAILED_THROW( GtoFreq2TDateInterval(
		(long) Basis(),
		&payInterval));

	for (idx=0; idx<= NumItems()-1; idx++) {

	    // Recalculate  maturity dates 

	    IF_FAILED_THROW( GtoDateSubtract(
		Date(idx),
		oldValueDate,
		&matInterval));

	    IF_FAILED_THROW( GtoDtFwdAny(
		BaseDate(),
		&matInterval,
		&Date(idx)));


	    // Compute par yields 
	    switch (toupper(mRateType[idx])) {
	    case 'M':
		IF_FAILED_THROW( GtoGetDayCountConv(
			'A',
			mMmDen,
			&mmDayCount));

                IF_FAILED_THROW( GtoZerosToSimplePoint(
			zcCurve,
                        ZeroInterpType(),
			BaseDate(),
                        Date(idx),
                        mmDayCount,
                        &Rate(idx)));

		break;
	    case 'S':
                IF_FAILED_THROW( GtoZerosToCouponsPoint(
			zcCurve,
                        ZeroInterpType(),
                        BaseDate(),
                        &payInterval,
                        Date(idx),
			DayCc(),
                        GTO_STUB_BOND,
                        FALSE,
                        &Rate(idx)));
                break;
	    default:
		throw KFailure("%s: bad rate type.\n", routine);
	    }
	}

    }
    catch (...) {
	throw KFailure("%s: failed.\n", routine);
    }
}


//---------------------------------------------------------------

double
KPCurve::Get10YEq()
{
	return(0e0);
}

//---------------------------------------------------------------


double
KPCurve::EffDur(const KDateInterval& intvl, int upDown) 
{
static	char	routine[] = "KPCurve::EffDur";
	int	idx;

	TDate	maturityDate;
	double	yieldSave, yield, effDur, effDur2, zeroRate,
		rateTwkSize = 1e-4;
	TCurve	*zcCurve = NULL;

    try {
	// FInd maturity date
	maturityDate = BaseDate() + intvl;

	for (idx=0; idx<= NumItems()-1; idx++)
		if (Date(idx) == maturityDate) break;
	if (idx >= NumItems()) 
		throw KFailure("%s: can't find intvl %s.\n",
				routine,
				intvl.Str()); 


	    //
	    // Up Tweak
	    //

	    yieldSave = Rate(idx);
	    Rate(idx) = yieldSave + rateTwkSize;

	    // Generate Zero Curve
	    zcCurve = ZCGenNPI(*this);

	    yield = Rate(idx);
	    maturityDate = Date(idx);

	    switch (toupper(mRateType[idx])) {
	    case 'M':
		IF_FAILED_THROW(GtoMoneyMarketModDuration(
			yield,
			maturityDate-BaseDate(),
			mMmDen,
			&effDur));
		break;
	    case 'S':
		IF_FAILED_THROW(GtoInterpRate(
			maturityDate,
			zcCurve,
			ZeroInterpType(),
			&zeroRate));

		IF_FAILED_THROW(GtoBondEffDuration(
			yield,
			zcCurve->fBasis,
			(maturityDate-BaseDate())/365e0,
			zeroRate,
			0e0,
			&effDur));
			
		break;
	    default:
		throw KFailure("%s: bad rate type.\n", routine);
	    }
	    GtoFreeTCurve(zcCurve); zcCurve = NULL;

	    //
	    // Down Tweak
	    //
	    if (upDown != FALSE) {
	    	Rate(idx) = yieldSave - rateTwkSize;

	    	// Generate Zero Curve 
	        zcCurve = ZCGenNPI(*this);

	    	yield = Rate(idx);
	    	maturityDate = Date(idx);

	    	switch (toupper(mRateType[idx])) {
	    	case 'M':
			IF_FAILED_THROW(GtoMoneyMarketModDuration(
				yield,
				maturityDate-BaseDate(),
				mMmDen,
				&effDur2));
			break;
	    	case 'S':
			IF_FAILED_THROW(GtoInterpRate(
				maturityDate,
				zcCurve,
				ZeroInterpType(),
				&zeroRate));
	
			IF_FAILED_THROW(GtoBondEffDuration(
				yield,
				zcCurve->fBasis,
				(maturityDate-BaseDate())/365e0,
				zeroRate,
				0e0,
				&effDur2));
			
			break;
	    	default:
			throw KFailure("%s: bad rate type.\n", routine);
	    	}
	    	GtoFreeTCurve(zcCurve); zcCurve = NULL;

		effDur = 0.5e0 * (effDur + effDur2);

	    }


	    /* Restore old yield */
	    Rate(idx) = yieldSave;

	return (effDur);
    }
    catch (KFailure) {
	GtoFreeTCurve(zcCurve);
	throw KFailure("%s: failed\n", routine);
    }
}







//---------------------------------------------------------------
//

TDateInterval
KPCurve::Intvl(int i) const
{
	TDateInterval	intvl;
	IF_FAILED_THROW( GtoDateSubtract(
		Date(i),
   		BaseDate(),
		&intvl));
	return(intvl);
}

//---------------------------------------------------------------
//


KMap(KDateInterval,double)
KPCurve::IntvlMap() const
{
	int	idx;
	KMap(KDateInterval,double)	matValMap;
	for (idx=0; idx<= NumItems()-1; idx++) {
		matValMap[Intvl(idx)] = Rate(idx);
	}	
	return matValMap;
}


//---------------------------------------------------------------


double
KPCurve::InterpYears(double yrs) const
{
static	char	routine[] = "KPCurve::InterpYears(double)";
	double	value;
	TDate	date;
    try {
	IF_FAILED_THROW( GtoTDateAdvanceYears(
		BaseDate(),
		yrs,
		&date));

	IF_FAILED_THROW( GtoInterpRate(
		date,
		mZcCurve,
		ZeroInterpType(),
		&value));
	return(value);
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}



//---------------------------------------------------------------
// 

void
KPCurveRegress(
	int numFact,			// (I) number of factors
	KPCurve *fact,			// (I) array of factors
	KPCurve &weiMat,		// (I) weights (or NULL) 
	KPCurve &beforeMat,		// (I) start (or NULL) 
	KPCurve &afterMat,		// (I) end 
	double *amplitudes,		// (O) amplitudes [0..numFact-1] 
	double *residual)		// (O) residual 
{
static	char	routine[] = "KPCurveRegress";

	int	i, j;
        int	m;		// (I) number of equations (# points)
        int	n;		// (I) number of unknowns (# fact)
        double	**a = NULL;	// (I) input matrix [0..m-1][0..n-1] 
        double	*x = NULL;	// (O) solution [0..n-1] 
        double	*b = NULL;	// (I) inhomogeneous term [0..m-1] 
        int	regType = 0;	// (I) see DrlMatrixSvdRegularizeEV 
        double	alpha = 1e-9;	// (I) see DrlMatrixSvdRegularizeEV 
        double	param = 0e0;	// (I) see DrlMatrixSvdRegularizeEV 
	double	res0, rfct;

    try {
	// 
	n = numFact;
	m = beforeMat.NumItems();

	// Check consistency
	if (afterMat)
		CheckSameType(beforeMat, afterMat);
	for (j=0; j<=n-1; j++) 
		CheckSameType(beforeMat, fact[j]);

	// 
	ASSERT_OR_THROW((a = DrlDoubleMatrAlloc(0, m-1, 0, n-1)) != NULL);
	ASSERT_OR_THROW((x = DrlDoubleVectAlloc(0, n-1)) != NULL);
	ASSERT_OR_THROW((b = DrlDoubleVectAlloc(0, m-1)) != NULL);

	for (i=0; i<=m-1; i++) {
	    if (weiMat) {
		for (j=0; j<=n-1; j++) {
			a[i][j] = fact[j].Rate(i) *
				sqrt(weiMat.Rate(i));
		}
		b[i] = afterMat.Rate(i) *
				sqrt(weiMat.Rate(i));
		if (beforeMat) {
			b[i] = beforeMat.Rate(i) *
				sqrt(weiMat.Rate(i));
		}
	    } else {
		for (j=0; j<=n-1; j++) {
			a[i][j] = fact[j].Rate(i);
		}
		b[i] = afterMat.Rate(i);
		if (beforeMat) {
			b[i] -= beforeMat.Rate(i);
		}
	    }
	}


	IF_FAILED_THROW (DrlRealLinSysSvd(
		m,
		n,
		a,
		x,
		b,
		regType,
		alpha,
		param));

	for (j=0; j<=n-1; j++) {
		amplitudes[j] = x[j];
	}

	// Compute residual : absolute everage deviation 
	*residual = 0e0;
	res0 = 0e0;
	for (i=0; i<=m-1; i++) {
	    rfct = 0e0;
	    for (j=0; j<=n-1; j++) {
		rfct += a[i][j]*x[j];
	    }
	    *residual += fabs(b[i]-rfct);
	    res0 += fabs(b[i]);
	}
	*residual /= m;
	res0 /= m;

	dppLog << format("%s: RES0=%8.4f %%   RES1=%8.4f %%\n",
			routine, res0*1e2, *residual*1e2);

	DrlDoubleMatrFree(a, 0, n-1, 0, m-1);
	DrlDoubleVectFree(x, 0, n-1);
	DrlDoubleVectFree(b, 0, m-1);
	return;

    }
    catch (KFailure) {
	DrlDoubleMatrFree(a, 0, n-1, 0, m-1);
	DrlDoubleVectFree(x, 0, n-1);
	DrlDoubleVectFree(b, 0, m-1);
	throw KFailure("%s: failed.\n");
    }
}















