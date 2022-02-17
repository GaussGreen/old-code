/***************************************************************
 * Module:	PenGuin
 * Submodule:	
 * File:	ktypes.h
 * Function:	Standard Definition and Include files.
 * Author:	Christian Daher
 * Revision:	$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/crxdpp/include/ktypes.h,v 1.1.1.1 2005/06/27 19:16:10 dliu Exp $
 ***************************************************************/
#ifndef	_ktypes_H
#define	_ktypes_H

#include "kstdinc.h"		/* Exceptions, io, etc. */
#include "kstub.h"


extern	"C" {
#include "drlstd.h"

#include "ldate.h"
#include "convert.h"
#include "cmemory.h"
#include "macros.h"
#include "bastypes.h"

#include "date_sup.h"
};

// For compatibility with ALIB 9.2
#if !defined(GTO_DATE_ADJ_TYPE_CALENDAR)
# define GTO_DATE_ADJ_TYPE_CALENDAR 0
# define GTO_DATE_ADJ_TYPE_BUSINESS 1
#endif



//--------------------------------------------------------------
/**
 * A class that wraps the ALIB types
 * TDateInterval (for calendar intervals) and
 * TDateAdjIntvl (for business days intervals)
 * so that it can be used with the standard containers.
 */

class	KDateInterval {
public:
	/** Default constructor.  Creates a null interval. */
	KDateInterval()
	{
		IF_FAILED_THROW( GtoYearsToDateInterval(
			0e0,
			&mDateAdjIntvl.interval));
		mDateAdjIntvl.isBusDays = FALSE;
		mDateAdjIntvl.holidayFile = NULL;
		SetHoliday("NONE");
		mDateAdjIntvl.badDayConv = GTO_BAD_DAY_NONE;
	}


	/**
	 * Copy constructor.
	 */
	KDateInterval(const KDateInterval &intvl)
	{
		mDateAdjIntvl.interval = intvl.mDateAdjIntvl.interval;
		mDateAdjIntvl.isBusDays = intvl.mDateAdjIntvl.isBusDays;
		mDateAdjIntvl.holidayFile = NULL;
		SetHoliday(intvl.mDateAdjIntvl.holidayFile);
		mDateAdjIntvl.badDayConv = intvl.mDateAdjIntvl.badDayConv;
	}



	/** Creates a KDateInterval from a TDateInterval. */
	KDateInterval(const TDateInterval &intvl)
	{
		mDateAdjIntvl.interval = intvl;
		mDateAdjIntvl.isBusDays = GTO_DATE_ADJ_TYPE_CALENDAR;
		mDateAdjIntvl.holidayFile = NULL; SetHoliday("NONE");
		mDateAdjIntvl.badDayConv = GTO_BAD_DAY_NONE;
	}

	/** Creates a KDateInterval from a TDateAdjIntvl. */
	KDateInterval(const TDateAdjIntvl &adjIntvl)
	{
		mDateAdjIntvl.interval = adjIntvl.interval;
		mDateAdjIntvl.isBusDays = adjIntvl.isBusDays;
		mDateAdjIntvl.holidayFile = NULL;
				SetHoliday(adjIntvl.holidayFile);
		mDateAdjIntvl.badDayConv = adjIntvl.badDayConv;
	}

	/** Creates a KDateInterval from a double. */
	KDateInterval(double yrs)
	{
		TDateInterval	intvl;
		IF_FAILED_THROW( GtoYearsToDateInterval(
			yrs,
			&intvl));
		mDateAdjIntvl.interval = intvl;
		mDateAdjIntvl.isBusDays = GTO_DATE_ADJ_TYPE_CALENDAR;
		mDateAdjIntvl.holidayFile = NULL; SetHoliday("NONE");
		mDateAdjIntvl.badDayConv = GTO_BAD_DAY_NONE;
	}


	/** Creates a KDateInterval from a number of days. */
	KDateInterval(int numDays, int isBusDays)
	{
		TDateInterval	intvl;
		IF_FAILED_THROW( GtoMakeDateInterval(
			numDays,
			'D',
			&intvl));

		mDateAdjIntvl.interval = intvl;
		mDateAdjIntvl.isBusDays = isBusDays;
		mDateAdjIntvl.holidayFile = NULL; SetHoliday("NONE");
		mDateAdjIntvl.badDayConv = GTO_BAD_DAY_NONE;
	}

	/**
	 * Creates a KDateInterval from two dates.
	 */
	KDateInterval(TDate date1, TDate date2)
	{
		TDateInterval	intvl;
		IF_FAILED_THROW( GtoDateSubtract(
			date1,
			date2,
			&intvl));

		mDateAdjIntvl.interval = intvl;
		mDateAdjIntvl.isBusDays = GTO_DATE_ADJ_TYPE_CALENDAR;
		mDateAdjIntvl.holidayFile = NULL; SetHoliday("NONE");
		mDateAdjIntvl.badDayConv = GTO_BAD_DAY_NONE;
	}


	/**
	 * Creates a KDateInterval from a string.
	 * Recognized format are 1D, 1M, etc..
	 */
	KDateInterval(const char *str);


	/**
	 * Destructor.
	 */
	~KDateInterval()
	{
		if (mDateAdjIntvl.holidayFile)
			FREE(mDateAdjIntvl.holidayFile);
	}



	/**
	 * Copy operator.
	 */
	KDateInterval& operator=(const KDateInterval& intvl)
	{
		mDateAdjIntvl.interval = intvl.mDateAdjIntvl.interval;
		mDateAdjIntvl.isBusDays = intvl.mDateAdjIntvl.isBusDays;
		SetHoliday(intvl.mDateAdjIntvl.holidayFile);
		mDateAdjIntvl.badDayConv = intvl.mDateAdjIntvl.badDayConv;
		return(*this);
	}


	/**
	 * Conversion to years.
	 * Throws an exception if the intervals is a business days interval.
	 */
	double Years() const
	{
		double	yrs;
		TDateInterval	intvl = mDateAdjIntvl.interval;
		IF_FAILED_THROW( GtoDateIntervalToYears(
			&intvl,
			&yrs));
		return(yrs);
	}

	/**
	 * Conversion to frequency (1,2,4,12).
	 * Throws an exception if the intervals is a business days interval.
	 */
	int Freq() const
	{
		double	freq;
		TDateInterval	intvl = mDateAdjIntvl.interval;
		IF_FAILED_THROW( GtoDateIntervalToFreq(
			&intvl,
			&freq));
		return((int) freq);
	}

	/**
	 * Sets the holiday definition.
	 */
	KDateInterval& SetHoliday(const char *holidayFile)
	{
		if (mDateAdjIntvl.holidayFile == NULL) {
			mDateAdjIntvl.holidayFile = NEW_ARRAY(char, 256);
			ASSERT_OR_THROW(mDateAdjIntvl.holidayFile != NULL);
		}
		if (holidayFile)
			strncpy(mDateAdjIntvl.holidayFile, holidayFile, 256);
		else
			strcpy(mDateAdjIntvl.holidayFile, "NONE");
		return (*this);
	}


	/**
	 * Conversion operator to TDateInterval.
	 * Throws an exception if the intervals is a business days interval.
	 */
	operator TDateInterval() const
	{
		return mDateAdjIntvl.interval;
	}

	/**
	 * Conversion operator to TDateInterval.
	 */
	operator const TDateAdjIntvl*() const
	{
		return &mDateAdjIntvl;
	}


	/**
	 * Equality test operator.
	 */
	bool operator==(const KDateInterval& intvl) const;

	/**
	 * Non equality test operator.
	 */
	bool operator!=(const KDateInterval& intvl) const
		{ return !(this->operator==(intvl)); }


	/**
	 * Compares two KDateIntervals.
	 * Warning: IMM intervals are treated as number of
	 * quarters.
	 */
	bool operator< (const KDateInterval&) const;


	/** Read object from stream. */
friend istream& operator>> (istream& is, KDateInterval &intvl);

	/** Write object to stream. */
friend ostream& operator<< (ostream& os, const KDateInterval &intvl);

	/** Returns a (static) string. */
	const char* Str() const;



	/**
	 * Reads and returns a KDateInterval from a stream.
	 * Skips comments.
	 */
friend	KDateInterval	getKDateInterval(istream& is,
				const char *debugName = "");


	/**
	 * Adds an interval to a date
	 */
friend	TDate	operator+(TDate date, const KDateInterval &intvl)
	{
		TDate		newDate;
		TDateAdjIntvl iv = intvl.mDateAdjIntvl;
		IF_FAILED_THROW( GtoDtFwdAdj(date, &iv, &newDate));
		return(newDate);
	}


	/**
	 * Subtracts an interval to a date
	 */
friend	TDate	operator-(TDate date, const KDateInterval &intvl)
	{
		TDate		newDate;
		TDateAdjIntvl	iv = intvl.mDateAdjIntvl;
		IF_FAILED_THROW( GtoDtBackAdj(date, &iv, &newDate));
		return(newDate);
	}

	/**
	 * Multiply an interval 
	 */
friend 	KDateInterval	operator*(int numIntvl, const KDateInterval &intvl)
	{
		KDateInterval newIntvl = intvl;
		newIntvl.mDateAdjIntvl.interval.prd *= numIntvl;
		return(newIntvl);
	}


private:
	TDateAdjIntvl	mDateAdjIntvl;
};


/**
 * Loads calendar files from the disk.
 * The input string "currencies" contains a list 3-letter
 * currency codes separted by the column (":") character
 * and "dir" the directory containg the holiday files
 * with name of the form "holiday.cur".\\
 * Example: \\
 *     DppHolidayLoadFromDisk("usd:gbp", "/opt/calendar") \\
 * looks for the files holiday.usd and holiday.gbp in the
 * directory /opt/calendar and loads the calendars
 * as "USD" and "GBP" respectively.\\
 * 
 * If a NULL "currencies" string is passed, it looks for
 * an environment variable DR_HOLIDAY_DIR which 
 * consists of the concatenation of a directory
 * and a list of currencies separated by ":".
 * Example: \\
 *     DR_HOLIDAY_DIR="/opt/calendar:usd:gbp".\\
 */


void	DppHolidayLoadFromDisk(
	const char *currencies,	// (I) list of cur separated by ":"
	const char *dir);	// (I) directory (or NULL)


//--------------------------------------------------------------
/**
 * A class for knock in/out definition.
 */

enum KKnockIO { CRX_KNOCK_IN, CRX_KNOCK_OUT, CRX_NONE };

istream& operator>>(istream& os, KKnockIO &x);
ostream& operator<<(ostream& os, const KKnockIO &x);

	/**
	 * Read data from an input stream.
	 *
	 * Scans a KKnockIO value in an istream and returns it.
	 * Skips all commented lines (#) and ignores double quotes (")
	 * around strings.
	 * @exception KFailure
	 */
	KKnockIO	getKnockIO(istream& is, const char *debugName = "");



//--------------------------------------------------------------
/**
 * A class for smoothing type definition.
 */

enum KSmooth { SINGLE_SMOOTH, DOUBLE_SMOOTH, NO_SMOOTH };

istream& operator>>(istream& os, KSmooth &x);
ostream& operator<<(ostream& os, const KSmooth &x);

	/**
	 * Read data from an input stream.
	 *
	 * Scans a KSmooth value in an istream and returns it.
	 * Skips all commented lines (#) and ignores double quotes (")
	 * around strings.
	 * @exception KFailure
	 */
	KSmooth		getSmooth(istream& is, const char *debugName = "");



//--------------------------------------------------------------
/**
 * A class for observation freqency type definition.
 */

enum KObsType { DATES, FREQ, NODES, CONTINUOUS };

istream& operator>>(istream& os, KObsType &x);
ostream& operator<<(ostream& os, const KObsType &x);

	/**
	 * Read data from an input stream.
	 *
	 * Scans a KObsType value in an istream and returns it.
	 * Skips all commented lines (#) and ignores double quotes (")
	 * around strings.
	 * @exception KFailure
	 */
	KObsType	getKObsType(istream& is, const char *debugName = "");






//--------------------------------------------------------------
/**
 * A class for day count convention.
 * Wraps the ALIB long type.
 */

class KDayCc {
public:
	/**
	 * Default constructor.
	 */
	KDayCc()
	{}

	/**
	 * Constructor from ALIB long day count conv.
	 */
	KDayCc(long alibDcc)
	{
		mDayCc = alibDcc;
	}

	/**
	 * Constructor from a string.
	 */
	KDayCc(const char *);


	/**
	 * Cast to ALIB long day count conv.
	 */
	operator long()
	{
		return(mDayCc);
	}

	/**
	 * Default accrual in risky cashflow is specified as same integer
     * representation for dcc but negative.
	 */
	bool isNegative() const
	{
		return(mDayCc < 0 ? true : false);
	}



	/**
	 * Read data from a stream.
	 * Skips all commented lines (#) and ignores double quotes (")
	 * around strings.
	 * @exception KFailure
	 */
	istream& Get(istream& is, const char *debugName = "");


	/**
	 * Reads from a stream.
	 */
friend	istream& operator>>(istream& is, KDayCc &dcc)
	{
		return dcc.Get(is);
	}

	/**
	 * Writes to a stream.
	 */
friend	ostream& operator<<(ostream& os, const KDayCc &dcc);


	/**
	 * Calculate the day count fraction.
	 */
friend	double	DayCountFraction(TDate startDate, TDate endDate, KDayCc &dcc);

private:
			/** Wrapped long. */
	long	mDayCc;
};


//--------------------------------------------------------------
/**
 * A class for stub convention.
 */

class KStubConv {
public:
	/**
	 * Default constructor.
	 */
	KStubConv()
	{mStubConv = KV_STUB_PAR;}

	/**
	 * Constructor from ALIB long day count conv.
	 */
	KStubConv(long alibStubConv)
	{
		mStubConv = alibStubConv;
	}

	/**
	 * Constructor from a string.
	 */
	KStubConv(const char *s);

	/**
	 * Cast to ALIB long day count conv.
	 */
	operator long() const
	{
		return(mStubConv);
	}

	/**
	 * Read data from a stream.
	 * Skips all commented lines (#) and ignores double quotes (")
	 * around strings.
	 * @exception KFailure
	 */
	istream& Get(istream& is, const char *debugName = "");

	/**
	 * Reads from a stream.
	 */
friend	istream& operator>>(istream& is, KStubConv &stubConv)
	{
		return stubConv.Get(is);
	}

	/**
	 * Writes to a stream.
	 */
friend	ostream& operator<<(ostream& os, const KStubConv &dcc);


private:
	long	mStubConv;
};





#endif
