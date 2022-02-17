/****************************************************************
 * Module:	PenGuin
 * Submodule:	
 * File:	
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#include <ctype.h>
#include "kzcurve.h"

#include "kutilios.h"		// 

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
#include "drlts.h"

#include "drltime.h"
#include "drllineq.h"		/* DrlRealLinSysSvd */
};


//---------------------------------------------------------------
// Default constructor.

KZCurve::KZCurve()
{
	mFormat = KZCurve::STD;
	mZcCurve = NULL;
	mZeroShift = 1e0;
	mZeroInterpType = GTO_LINEAR_INTERP;
}

//---------------------------------------------------------------
// Copy constructor.


KZCurve::KZCurve(const KZCurve& swpcrv)
{

	mFormat = KZCurve::STD;
	mZcCurve = NULL;
	mZeroShift = 1e0;
	if (swpcrv.mZcCurve == NULL) return;


	create(
		swpcrv.BaseDate(),
		swpcrv.NumItems(),
		swpcrv.Basis(),
		swpcrv.DayCc());

	for (int i=0; i<=NumItems()-1; i++) {
	    Date(i) = swpcrv.Date(i);
	    Rate(i) = swpcrv.Rate(i);
	}

	mFormat = swpcrv.mFormat;
	mZeroShift = swpcrv.ZeroShift();

	mZeroInterpType = swpcrv.ZeroInterpType();
}

//---------------------------------------------------------------


//
//---------------------------------------------------------------

KZCurve::KZCurve(
	TDate baseDate,			// (I) 
	const KVector(TDate)& dates,	// (I) 
	const KVector(double)& rates,	// (I) 
	double basis,			// (I) 
	KDayCc dayCountConv,		// (I) 
	int    interpType)	        // (I)
{
static	char	routine[] = "KZCurve::KZCurve";

	if (dates.size() != rates.size()) {
		throw KFailure("%s: dates .size (%d) != rates size(%d.\n",
			routine, dates.size(), rates.size());
	}

	create(baseDate, dates.size(), basis, dayCountConv);

	for (int i=0; i<=NumItems()-1; i++) {
	    Date(i) = dates[i];
	    Rate(i) = rates[i];
	}

	mFormat = KZCurve::STD;
	mZeroShift = 1e0;

	mZeroInterpType = interpType;
}



//---------------------------------------------------------------
// Copy constructor.

/*

KZCurve::KZCurve(const TCurve* zcCurve)
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
	}

	mFormat = KZCurve::STD;
	mZeroShift = 1e0;

	mZeroInterpType = GTO_LINEAR_INTERP;
}

*/

//---------------------------------------------------------------
// Copy constructor with interp type.


KZCurve::KZCurve(const TCurve* zcCurve, int interpType)
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
	}

	mFormat = KZCurve::STD;
	mZeroShift = 1e0;

	mZeroInterpType = interpType;
}




//---------------------------------------------------------------
// Copy constructor. Shift base date from value date to todayDt.
//

KZCurve::KZCurve(TCurve* zcCurve, TDate todayDt, int interpType)
{
static	char	routine[] = "KZCurve::KZCurve";

	int	i, j, numItems;
	double	discount;

    ASSERT_OR_THROW (zcCurve != NULL);
    ASSERT_OR_THROW (zcCurve->fNumItems > 0);
    ASSERT_OR_THROW (zcCurve->fArray != NULL);
    
    mZeroInterpType = interpType;

    // initialize;
	mZeroShift = 1e0;

	if (todayDt == zcCurve->fBaseDate)	// just copy
	{
		create(
			zcCurve->fBaseDate,
			zcCurve->fNumItems,
			zcCurve->fBasis,
			zcCurve->fDayCountConv);

		for (i=0; i<=NumItems()-1; i++) {
	    		Date(i) = zcCurve->fArray[i].fDate;
	    		Rate(i) = zcCurve->fArray[i].fRate;
		}

		mFormat = KZCurve::STD;
		mZeroShift = 1e0;
	}
	else if (todayDt < zcCurve->fBaseDate)	// insert value date in zc
	{
		mFormat = KZCurve::STD;

		// Insert value date in zc if the first date > value date
		//
		if (zcCurve->fBaseDate < zcCurve->fArray[0].fDate)
		{
			numItems = zcCurve->fNumItems + 1;

			create(
				todayDt,
				numItems,
				zcCurve->fBasis,
				zcCurve->fDayCountConv);

			// The first date is value date with 
			// flat rate.
			// 
			Date(0) = zcCurve->fBaseDate;
			Rate(0) = zcCurve->fArray[0].fRate;

		}
		else {
			// Only shift the base date to today,
			// but the first date and rate are unchanged

			numItems = zcCurve->fNumItems;

			create(
				todayDt,
				numItems,
				zcCurve->fBasis,
				zcCurve->fDayCountConv);

			// 
			Date(0) = zcCurve->fArray[0].fDate;
			Rate(0) = zcCurve->fArray[0].fRate;

		}

		// Zero from today to first date
		//
		IF_FAILED_THROW(GtoRateToDiscount(
				Rate(0),
				todayDt,
				Date(0),
				zcCurve->fDayCountConv,
				(long) zcCurve->fBasis,
				&mZeroShift));

		//
		//
		for (i=1; i<=NumItems()-1; i++) {
			if (zcCurve->fBaseDate < zcCurve->fArray[0].fDate)
				j = i - 1;
			else
				j = i;

	    		Date(i) = zcCurve->fArray[j].fDate;
			
			// discount relative to fBaseDate	
			//
			IF_FAILED_THROW(GtoRateToDiscount(
					zcCurve->fArray[j].fRate,
					zcCurve->fBaseDate,
					zcCurve->fArray[j].fDate,
					zcCurve->fDayCountConv,
					(long) zcCurve->fBasis,
					&discount));

			// discount relative to todayDt
			//
			discount *= mZeroShift;
	
			// convert to rate
			//
			IF_FAILED_THROW(GtoDiscountToRate(
					discount,
					todayDt,
					Date(i),
					zcCurve->fDayCountConv,
					(long) zcCurve->fBasis,
					&Rate(i)));
		}
	}
	else
		throw KFailure("%s: today's date (%s) > "
			      "zero curve base date (%s).\n",
			      routine, GtoFormatDate(todayDt),
			      GtoFormatDate(zcCurve->fBaseDate));

}



//---------------------------------------------------------------
// private memory management

void
KZCurve::create(TDate baseDate, int size, double basis,
	TDayCount dayCountConv)
{
    mZcCurve = NULL;
	mZcCurve = GtoNewTCurve(baseDate, size, basis, dayCountConv);
        if (!mZcCurve) throw KFailure("GtoNewTCurve failed.\n");
}



void
KZCurve::destroy()
{
	GtoFreeTCurve(mZcCurve);
	mFormat = 0;
	mZcCurve = NULL;
}



//---------------------------------------------------------------
// Destructor.

KZCurve::~KZCurve()
{
	destroy();
}


//---------------------------------------------------------------
// Assignment operator.

KZCurve&
KZCurve::operator=(const KZCurve& swpcrv)
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


	for (i=0; i<=NumItems()-1; i++) {
	    Date(i) = swpcrv.Date(i);
	    Rate(i) = swpcrv.Rate(i);
	}

	mZeroShift = swpcrv.ZeroShift();
	mZeroInterpType = swpcrv.ZeroInterpType();

	return(*this);
}




//---------------------------------------------------------------
// Assignment operator.

KZCurve&
KZCurve::operator=(const TCurve* zcCurve)
{
	int	i;

	destroy();
	if (!zcCurve) return(*this);
	create(
		zcCurve->fBaseDate,
		zcCurve->fNumItems,
		zcCurve->fBasis,
		zcCurve->fDayCountConv);

	for (i=0; i<=NumItems()-1; i++) {
	    Date(i) = zcCurve->fArray[i].fDate;
	    Rate(i) = zcCurve->fArray[i].fRate;
	}

	mZeroShift = 1.0;
	mZeroInterpType = GTO_LINEAR_INTERP;

	return(*this);
}


//---------------------------------------------------------------
// Read class from a stream.

istream&
KZCurve::Get(istream& is)
{
static	char	routine[] = "KZCurve::Get";

	int		i;
	TDate		xBaseDate;		// tmp storage
	double		xBasis;
	TDayCount	xDayCountConv;
	int		xNumItems;


    try {

	xBaseDate = getTDate(is);
	xBasis = getDouble(is);
	xDayCountConv = getTDayCount(is);
	xNumItems = getInt(is);

	destroy();
	create(
		xBaseDate,
		xNumItems,
		xBasis,
		xDayCountConv);

	for (i=0; i<=NumItems()-1; i++) {
		Date(i) = getTDate(is);
		getString(is);
		Rate(i) = getDouble(is);
	}

	mZeroInterpType = GTO_LINEAR_INTERP;

	return(is);
    }
    catch(...) {
	    throw KFailure("%s: failed.\n", routine);
    }
}



//---------------------------------------------------------------
// Write class to a stream.

ostream&
KZCurve::Put(ostream& os, int indent) const
{
	int	i;

	if (this->isempty()) {
		os << "EMPTY\n";
		return(os);
	}

	switch (this->mFormat) {
	case KZCurve::BVDRWRAP:
		os << "# Base Volatility Frequency (A,S,Q,M)\n";
		switch (this->Freq()) {
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
			    this->mZcCurve->fArray[i].fDate),
			this->mZcCurve->fArray[i].fRate*1e2);
		}

	    break;
	case KZCurve::STD:
	default:
	    os << "#mBaseDate\n" << GtoFormatDate(this->BaseDate()) << '\n';
	    os << "#mBasis\n" << Basis() << '\n';
	    os << "#mDayCountConv\n"
		<< GtoFormatDayCountConv(DayCc()) << '\n';
	    os << "#mZeroInterpType\n"
		<< GtoFormatInterpType(this->ZeroInterpType()) << '\n';
		    	

	    os << "#mNumItems\n" << NumItems() << '\n';

	    for (i=0; i<=NumItems()-1; i++) {
		TDateInterval	intvl;

		IF_FAILED_THROW( GtoDateSubtract(
			this->Date(i),
			this->BaseDate(),
			&intvl));

		os << format(" %10s  %8s  %12.10f\n", 
			DrlTDatePrint(NULL, this->Date(i)),
			DrlTDateIntervalPrint(NULL, intvl),
			Rate(i));
	    }
	    break;
	}


	return(os.flush());
}


//---------------------------------------------------------------
//

int
CheckSameType(const KZCurve& crv1, const KZCurve& crv2)
{
	int	i;
	TDateInterval	intvl1, intvl2;
    try {
	ASSERT_OR_THROW(
		( crv1.isempty() &&  crv2.isempty()) ||
		(!crv1.isempty() && !crv2.isempty()));

	ASSERT_OR_THROW(crv1.ZeroInterpType() == crv2.ZeroInterpType());

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
	throw KFailure("CheckSameType(KZCurve&,KZCurve&): inconsistent.\n");
    }
}



//---------------------------------------------------------------
//

KZCurve&
KZCurve::operator=(double argument)
{
	int	i;
	if (isempty()) return (*this);
	for (i=0; i<= NumItems()-1; i++)
		Rate(i) = argument;
	return(*this);
}

//---------------------------------------------------------------
//

KZCurve&
KZCurve::operator+=(double argument)
{
	int	i;
	if (isempty()) return (*this);
	for (i=0; i<= NumItems()-1; i++)
		Rate(i) += argument;
	return(*this);
}

//---------------------------------------------------------------
//

KZCurve&
KZCurve::operator*=(double argument)
{
	int	i;
	if (isempty()) return (*this);
	for (i=0; i<= NumItems()-1; i++)
		Rate(i) *= argument;
	return(*this);
}


//---------------------------------------------------------------
//

KZCurve&
KZCurve::operator+=(const KZCurve& crv)
{
	int	i;

	CheckSameType(*this, crv);
	for (i=0; i<= NumItems()-1; i++) 
		Rate(i) += crv.Rate(i);
	return(*this);
}


//---------------------------------------------------------------

KZCurve&
KZCurve::ShiftDate(TDate baseDate, int forward)
{
static	char	routine[] = "KZCurve::ShiftDate";

    try {
	if (forward)  {
#ifdef	_SKIP

		long	diffDays;
		int	idx;
		double	z0, z1;
		TCurve	*newCurve = NULL;
	

		if ((diffDays = baseDate - BaseDate()) < 0)
			throw KFailure(
				"%s: new base date (%s) < base Date (%s).\n",
				GtoFormatDate(baseDate),
				GtoFormatDate(BaseDate()));

		ASSERT_OR_THROW((newCurve = GtoCopyCurve(mZcCurve)) != NULL);
		IF_FAILED_THROW( GtoDiscountDate(
			baseDate,
			mZcCurve,
			mZeroInterpType,
			&z0));

		newCurve->fBaseDate = baseDate;

		for (idx=0; idx<=mZcCurve->fNumItems-1; idx++) {
		    newCurve->fArray[idx].fDate += diffDays;
		    IF_FAILED_THROW( GtoDiscountDate(
			newCurve->fArray[idx].fDate,
			mZcCurve,
			mZeroInterpType,
			&z1));

		    IF_FAILED_THROW( GtoDiscountToRate(
			z1/z0,
			newCurve->fBaseDate,
			newCurve->fArray[idx].fDate,
			newCurve->fDayCountConv,
			(long)newCurve->fBasis,
			&newCurve->fArray[idx].fRate));
		}
		GtoFreeTCurve(mZcCurve);
		mZcCurve = newCurve;
#else


		long	diffDays;
		int	idx, idxOffset;
		double	z0, z1;
		KZCurve	newCurve;
	

		if ((diffDays = baseDate - BaseDate()) < 0)
			throw KFailure(
				"%s: new base date (%s) < base Date (%s).\n",
				GtoFormatDate(baseDate),
				GtoFormatDate(BaseDate()));

		for (idx=0; idx<=mZcCurve->fNumItems-1; idx++) {
			if (Date(idx) > baseDate)
				break;
		}
		idxOffset = idx;


		if (idx == mZcCurve->fNumItems)
			throw KFailure();


		newCurve.create(
			baseDate,
			mZcCurve->fNumItems-idxOffset,
			mZcCurve->fBasis,
			mZcCurve->fDayCountConv);



		IF_FAILED_THROW( GtoDiscountDate(
			baseDate,
			mZcCurve,
			mZeroInterpType,
			&z0));



		for (idx=0; idx<newCurve.NumItems(); idx++) {
			newCurve.Date(idx) = Date(idx+idxOffset);
			newCurve.Rate(idx) = Rate(idx+idxOffset);

		    IF_FAILED_THROW( GtoDiscountDate(
			Date(idx+idxOffset),
			mZcCurve,
			mZeroInterpType,
			&z1));

		    IF_FAILED_THROW( GtoDiscountToRate(
			z1/z0,
			newCurve.BaseDate(),
			newCurve.Date(idx),
			newCurve.DayCc(),
			(long)newCurve.Basis(),
			&newCurve.Rate(idx)));
		}



		*this = newCurve;


#endif

/*
		// Forward zc
		IF_FAILED_THROW(DrlTCurveForwardDateTCurve(
			mZcCurve,
			baseDate,
			&newZcCurve));

		GtoFreeTCurve(mZcCurve);
		mZcCurve = newZcCurve;
*/

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
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}


//---------------------------------------------------------------
// Shifts all rates by {\tt value}.

void
KZCurve::ShiftRates(double value)
{
static	char	routine[] = "KZCurve::ShiftRates";
	int	idx;

	for (idx=0; idx<= NumItems()-1; idx++) {
	    Rate(idx) += value;
	}

}


//---------------------------------------------------------------
//

KDateInterval
KZCurve::Intvl(int i) const
{
	TDateInterval	intvl;
	IF_FAILED_THROW( GtoDateSubtract(
		this->Date(i),
   		this->BaseDate(),
		&intvl));
	return(KDateInterval(intvl));
}


//---------------------------------------------------------------

double
KZCurve::InterpYears(double yrs) const
{
static	char	routine[] = "KZCurve::InterpYears(double)";
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
		mZeroInterpType,
		&value));
	return(value);
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}




//---------------------------------------------------------------

double
KZCurve::DiscFact(TDate discDate) const
{
static	char	routine[] = "KZCurve::DiscFact";
	double	value;

    try {
	IF_FAILED_THROW( GtoDiscountDate(
		discDate,
		mZcCurve,
		mZeroInterpType,
		&value));

	return(value);
    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}














