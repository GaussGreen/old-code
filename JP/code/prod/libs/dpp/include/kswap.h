/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	kswap.h
 * Function:	
 * Author:	David Liu
 *****************************************************************/
#ifndef	_kbmswap_H
#define	_kbmswap_H

#include "ktypes.h"
#include "krate.h"


//--------------------------------------------------------------
/**
 * This class describes the calibration of spot volatilities to
 * price a series of benchmark swaptionsa correctly.
 * It uses Newton-Raphson's method to solve iteratively.
 */


class KBMSwap {
public:

	/**
	 * Default constructor.
	 */
	KBMSwap();


	/**
	 * General constructor.
	 */
	KBMSwap(const char *name,		// (I) object name
		const TDate today,		// (I) today's date
		const TDate startDate,		// (I) start date
		const TDate matDate,		// (I) maturity date
		const KDateInterval &freq,	// (I) frequency as interval
		const double 	vol,		// (I) volatility
		KDayCc &dayCc,			// (I) pay day count
		TBoolean   stubAtEnd,		// (I) stub at end */
		const String  &curveName,	// (I) curve name
		const KZCurve &zc);             // (I) zero curve

	/**
	 * Destructor.
	 */
	~KBMSwap();

	/**
	 * Initialize zero resets and par coupons.
	 */
void	Initialize();
	
	/**
	 * Calculate the par yield of swap
	 */
double	ParYield();

	/**
	 * Calculate the B-S price 
	 */
double	BSPrice();

	/**
	 * Calculate the B-S vega 
	 */
double	BSVega();

	/**
	 * Get name
	 */
const String&   GetName() const
	{ return mName;}

	/**
	 * Get curve name
	 */
const String&   GetCurveName() const
	{ return mCurveName;}

	/**
	 * Calculate list of zero resets
	 */
void	GetZeroResetList(KVector(KZeroReset)&);

	/**
	 * Swap start date
	 */
TDate	StartDate()
	{ return mStDate; }


	/**
	 * Writes the object to a stream.
	 */
virtual ostream& Put(ostream& os) const;
 
        /**
         * Writes to a stream.
         */
friend  ostream& operator<<(ostream& os, const KBMSwap &swap)
	{ swap.Put(os); return (os);}



public:
				/** Array of pay dates. */
	KVector(TDate)	mPayDates;
				/** Par coupons */
	KVector(double)	mCoupons;



private:
				/** Name */
	String		mName;	
				/** Today's date */
	TDate		mTodayDate;
				/** Swap Start date */
	TDate		mStDate;
				/** Frequency */
	KDateInterval 	mFreq;
				/** volatility */
	double		mVol;
				/** Day count convention. */
	KDayCc		mDayCc;
				/** Stub at end. */
	TBoolean	mStubAtEnd;
				/** Curve name */
	String		mCurveName;
				/** Zero curve. */
	const KZCurve	*mZC;
};


#endif	/* _kbmswap_H */

