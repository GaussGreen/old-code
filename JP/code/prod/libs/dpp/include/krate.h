/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	krate.h
 * Function:	
 * Author:	David Liu, Christian Daher
 ***************************************************************/
#ifndef	_krate_H
#define	_krate_H
#include "kstdinc.h"
#include "kzcurve.h"
#include "ktypes.h"
#include "kvpatom.h"

extern "C" {
#include "bastypes.h"		// ALIB TFloatRate
};


const   int     K_DEFAULT_IDX = -1024;
const   String  K_DEFAULT_NAME("Default");

//
// Assumes a TFLoatRate ALIB implementation:
//
// typedef struct _TFloatRate
// {
//     TDateInterval   matInterval;   /* Time to maturity of rate */
//     TDateInterval   payInterval;   /* Time between payments for rate */
//     long            dayCountConv;  /* Day count convention of rate */
//     TDateAdjIntvl   spotOffset;    /* From reset to rate effective date */
//     double          spread;        /* Added to the rate  */
//     long            rateType;      /* GTO_SIMPLE_BASIS, GTO_ANNUAL_BASIS*/
//     double          weight;        /* Multiplied by rate */
// } TFloatRate;
//
//



//--------------------------------------------------------------
/**
 * Class for a floating rate definition.
 * @version 2.0
 */


class KRate : public KVPAtom {
public:
	/**
	 * Default constructor (constructs a fixed rate 0).
	 */
	KRate();

	/**
	 * General constructor.
	 */
	KRate(
		const KDateInterval &mat,	// (I) rate maturity
		const KDateInterval &freq,	// (I) rate payment freq
		KDayCc dayCc,			// (I) rate day count conv
		const KDateInterval &spotOffset,// (I) spot offset
		double spread,			// (I) spread
		double weight);			// (I) weight

	/**
	 * General constructor with curve name.
	 */
	KRate(
		const String  &curveName,	// (I) curve name
		const KDateInterval &mat,	// (I) rate maturity
		const KDateInterval &freq,	// (I) rate payment freq
		KDayCc dayCc,			// (I) rate day count conv
		const KDateInterval &spotOffset,// (I) spot offset
		double spread,			// (I) spread
		double weight);			// (I) weight

	/**
	 * General constructor for a simple rate between two
	 * specified dates.
	 */
	KRate(
		const String  &curveName,	// (I) curve name
		const TDate   startDate,	// (I) start date
		const TDate   endDate,		// (I) end date
		KDayCc dayCc,			// (I) rate day count conv
		const KDateInterval &spotOffset,// (I) spot offset
		double spread,			// (I) spread
		double weight);			// (I) weight


	/**
	 * General constructor with debugging and curve name.
	 */
	KRate(
		const String  &name,		// (I) debug name
		const String  &curveName,	// (I) curve name
		const KDateInterval &mat,	// (I) rate maturity
		const KDateInterval &freq,	// (I) rate payment freq
		KDayCc dayCc,			// (I) rate day count conv
		const KDateInterval &spotOffset,// (I) spot offset
		double spread,			// (I) spread
		double weight);			// (I) weight

	/**
	 * Copy constructor.
	 */
	KRate(const KRate& rate);

	/**
	 * Convenience constructor (creates a fixed rate).
	 */
	KRate(double fixedRate)
	{SetFixedRate(fixedRate);}

	/**
	 * Convenience constructor (creates a fixed rate with debugging name).
	 */
	KRate(String& name, double fixedRate)
	{mName = name; SetFixedRate(fixedRate);}


	/**
	 *
	 */
	~KRate()
	{
		FREE(mFloatRate.spotOffset.holidayFile);
	}

	/**
	 * Copy operator.
	 */
	KRate& operator=(const KRate& rate);


	/**
	 * Comparison operator (used by STL).
	 * Compares rate curve name, maturity, frequency and day count
	 * (the other fields are NOT compared).
	 */
	bool	operator<(const KRate& rate) const;
	
	/**
	 * Equality testing operator.
	 */
	bool	operator==(const KRate &rate) const;


	/**
	 * Type name.
	 */
virtual const char*	TypeName() const {return("KRate");}

	/**
	 * Returns TRUE/FALSE if the rate is floating/fixed.
	 */
	int	IsFloating() const
	{
		return !IS_ALMOST_ZERO(mFloatRate.weight);
	}

	/**
	 * Returns TRUE/FALSE if the rate is simple rate.
	 */
	int	IsSimple() const;

	/**
	 * Convenience fixed rate constructor 
	 */
	void SetFixedRate(double fixedRate);


	/**
	 *
	 */
	double	Spread() const
	{
		return mFloatRate.spread;
	}

	/**
	 *
	 */
	KRate& SetSpread(double value)
	{
		mFloatRate.spread = value;
		return *this;
	}


	/**
	 * Reads the object from a stream.
	 */
virtual	istream& Get(istream& is, int drwFmt = FALSE);

	/**
	 * Writes the object to a stream.
	 */
virtual	ostream& Put(ostream& os, int indent=FALSE) const;

	/**
	 * Writes in yacction format with the name.
	 */
virtual	ostream& YacctionWrite( ostream& os, int indent=FALSE);

	/**
	 * Writes in yacction format without the name.
	 */
virtual	ostream& WriteSimple( ostream& os);

	/**
	 * Reads from a stream.
	 */
friend	istream& operator>>(istream& is, KRate& rate)
		{ rate.Get(is); return (is);}


	/**
	 * Cast to TFloatRate.
	 */
	operator TFloatRate&()
	{
		return(mFloatRate);
	}


	/**
	 * Returns rate maturity as an interval.
	 */
	KDateInterval	Maturity() const
	{
		return KDateInterval(mFloatRate.matInterval);
	}


	/**
	 * Returns rate payment frequency as an interval.
	 */
	KDateInterval	Frequency() const
	{
		return KDateInterval(mFloatRate.payInterval);
	}

	/**
	 * Returns day count convention as a KDayCc.
	 */
	KDayCc		DayCc() const
	{
		return KDayCc(mFloatRate.dayCountConv);
	}

	/**
	 * Set the day count convention
	 */
	void    SetDayCc(KDayCc dcc)
	{
		mFloatRate.dayCountConv = (long)dcc;
	}

	/**
	 * Returns the days to spot interval.
	 */
	KDateInterval	SpotOffset() const
	{
		return KDateInterval(mFloatRate.spotOffset);
	}

	/**
	 * Get the curve name of the rate
	 */
	const String&   CurveName() const
		{ return mCurveName;}

	/**
	 * Set the days to spot interval in mFloatRate from date1 - date2.
	 */
	void	SetSpotOffset(int numDays, int isBusDays);

	/**
	 * Set the curve name.
	 */
	void	SetCurveName(const String curveName)
		{ mCurveName = curveName;}

	void	SetCurveName(const char *curveName)
		{ mCurveName = curveName;}


	/**
	 * Sets the name of KRate object.
	 */
	void	SetName(const char *name)
		{mName = name;}

	/**
	 * Sets the name of KRate object.
	 */
	void	SetName(const String name)
		{mName = name;}
 
 
	/**
	 * Gets the name of KRate object.
	 */
virtual const char*	GetName() const
	{
		return mName.c_str();
	}


	/**
	 * Convenience routine to compute a spot value of rate.
	 */
	double	Spot(const KZCurve& zcCurve) const;	// (I) zero curve

	/**
	 * Convenience routine to compute a forward.
	 */
	double	Forward(
		const KZCurve& zcCurve,		// (I) zero curve
		TDate rateEffDate) const;	// (I) rate effective date

	/**
	 * Convenience routine to compute a forward given a reset date
	 * (the spot offset is added to the reset date to obtain the effective
	 * date and and then calculated of the zero curve).
	 */
	double	Forward(
		TCurve *zcCurve,	// (I) zero curve
		TDate resetDate) const;	// (I) reset date

	/**
	 * Convenience routine to compute a forward CDS par spread.
	 */
	double	ForwardCDS(
		const KZCurve& irZCurve,	// (I) IR zero curve
		const KZCurve& crZCurve,	// (I) CDS zero curve
                double         recovery,        // (I) recovery   
		TDate rateEffDate) const;	// (I) rate effective date



	/**
	 * Convenience routine to compute a forward adjusted rate.
	 */
	double	ForwardAdjusted(
		const KZCurve& zcDiscount,// (I) Discount zero curve
		const KZCurve& zcIndex,	// (I) Zero curve for fwdrate
		TDate rateResetDate,	// (I) Reset date for forward rate
		TDate ratEffDate,	// (I) rate effective date
		TDate ratePayDate,	// (I) Payment date (for delay adj)
		TDate todayDate,	// (I) Where volatity starts
		double volFwdRate,	// (I) Vol of rate def. by fwdRateInfo
		double volResetToPay);	// (I) Vol of rate from reset to pay


private:
				/** Wrapped ALIB TFloatRate. */
	TFloatRate	mFloatRate;

protected:
				/** Curve name. */
	String		mCurveName;

};




//--------------------------------------------------------------
/**
 * Class for the definition of a zero coupon (i.e. discount factor)
 * reset. It contains the maturity of the zero coupon,
 * the date at which the zero coupon is to be observed (earliest date or
 * reset date) and the curve name over the discount factor is to be computed.
 */

class   KZeroReset : public KVPAtom
{
public:
	/** Constructor */
	KZeroReset(
		const String& curveName,// (I) Curve name.
		TDate earliestDate,	// (I) Maturity date.
		TDate maturityDate)	// (I) Earliest start date.
        {
                mCurveName = curveName;
                mMaturityDate = maturityDate;
                mEarliestDate = earliestDate;
        }

    /** Default constructor needed to initialize arrays */
    KZeroReset() {};

	/** Copy constructor */
	KZeroReset(const KZeroReset&);

	/**
	 * Type name.
	 */
virtual const char*	TypeName() const {return("KZeroReset");}

        /**
	 * Check the validity of dates and curve name.
	 */
virtual void    IsValid() const;


	/**
	 * Writes to a stream.
	 */
virtual ostream& Put(ostream& os, int indent=FALSE) const;

	/**
	 * Writes in yacction format.
	 */
virtual	ostream& YacctionWrite( ostream& os, int indent=FALSE)
		{Put(os, indent); return os;}

					/** Curve name. */
	String	mCurveName;
					/** Zero maturity date. */
	TDate	mMaturityDate;
					/** Earliest usage date. */
	TDate	mEarliestDate;
};


//--------------------------------------------------------------
/**
 * A class to represent a floating rate reset, i.e. a floating rate index
 * with specified reset and effetive dates.
 */

class	KRateReset : public KVPAtom {

public:	
	/**
	 * Constructor of KRateReset from a rate and
	 * a reset (i.e. observation) date.
	 * The effective date is computed using the spot offset
	 * convention contained in the KRate.
	 */
	KRateReset(
		TDate resetDate,		// (I) observation date
		const KRate& rate);		// (I) rate description


        /**
         * General constructor.
	 * Constructs a KRateReset from a KRate and reset information:
	 * curve name, reset date and effective date.
         */
        KRateReset(
		TDate resetDate,		// (I) observation date
		TDate effDate,			// (I) effective date
		const KRate& rate);		// (I) rate description

        KRateReset(
		const String& curveName,	// (I) curve name
		TDate resetDate,		// (I) observation date
		TDate effDate,			// (I) effective date
		const KRate& rate);		// (I) rate description

	/**
	 * Creates a simple stub rate reset between effDate and endDate,
	 * using the same day count convention, spread, etc. as "rate".
	 */
        KRateReset(
		TDate resetDate,		// (I) observation date
		TDate effDate,			// (I) effective date
		TDate endDate,			// (I) End date
		const KRate& rate);		// (I) rate description

	/** Copy constructor */
	KRateReset(const KRateReset&);

	/**
	 * Comparison operator (used by STL).
	 * Compares rate curve name, maturity, frequency, day count
	 * and reset date (the other fields are NOT compared).
	 */
	bool	operator<(const KRateReset& rate) const;

	/**
	 * Type name.
	 */
virtual const char*	TypeName() const {return("KRateReset");}

	/**
	 * Returns the underlying rate.
	 */
	const KRate& Rate() const
		{ return mRate; }

	/**
	 * Returns the underlying rate.
	 */
	KRate& Rate()
		{ return mRate; }

	/**
	 * Returns the reset date.
	 */
	TDate	ResetDate() const
		{ return mResetDate; }

	/**
	 * Returns the effective date.
	 */
	TDate	EffDate() const
		{ return mEffDate; }


	/**
	 * Generate a series of zero reset dates from rate reset.
	 */
	void	GetZeroResetList(KVector(KZeroReset)&) const; 

	/**
	 * Convenience routine to compute a forward.
	 */
	double	Forward(const KZCurve& zcCurve)	const // (I) zero curve
	{
		return mRate.Forward(zcCurve, mEffDate);
	}

	/**
	 * Writes to a stream.
	 */
virtual	ostream& Put(ostream& os, int indent=FALSE) const;

	/**
	 * Writes in yacction format.
	 */
virtual	ostream& YacctionWrite( ostream& os, int indent=FALSE)
		{Put(os, indent); return os;}

private:
					/** Floating rate definition. */
	KRate	mRate;
					/** Reset (observation) date. */
	TDate	mResetDate;
					/** Effective date. */
	TDate	mEffDate;

};


#endif  // _krate_H
