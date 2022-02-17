/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vpfleg.h
 * Function:	
 * Author:	David Liu
 *****************************************************************/
#ifndef	_vpfleg_H
#define	_vpfleg_H

#include "vpbase.h"
#include "krate.h"



//--------------------------------------------------------------
/**
 * This class describes a payment leg, i.e. a sequence of payments.
 * Currently it is modeled after the ALIB's TFloatLeg,
 * but this may change in the future
 * @author Christian Daher
 * @version 2.0
 */


class KVPFloatLeg : public KVPInstr {
public:

	/**
	 * General constructor with explicit effective reset dates.
	 */
	KVPFloatLeg(
		const char *name,			// (I) object name
		const KVector(TDate)& resetDates,	// (I) reset dates
		const KVector(TDate)& resetEffDates,	// (I) reset eff dates
		const KVector(TDate)& accStartDates,	// (I) acc start dates
		const KVector(TDate)& accEndDates,	// (I) acc end dates
		const KVector(TDate)& payDates,		// (I) pay dates array
		const KVector(double)& notionals,	// (I) notionals array
		KDayCc dayCc,				// (I) pay day count
		KStubConv stubConv,			// (I) stub conv
		const SharedPointer<KRate> floatRate,		// (I) paid rate
		const char *formula,			// (I) paymt formula 
		const char *discZcName);		// (I) disc curve name

	/**
	 * General constructor without explicit reset dates
	 * (embeded in floatRate.)
	 */
	KVPFloatLeg(
		const char *name,			// (I) object name
		const KVector(TDate)& resetEffDates,	// (I) eff reset dates
		const KVector(TDate)& accStartDates,	// (I) acc start dates
		const KVector(TDate)& accEndDates,	// (I) acc end dates
		const KVector(TDate)& payDates,		// (I) pay dates array
		const KVector(double)& notionals,	// (I) notionals array
		KDayCc dayCc,				// (I) pay day count
		KStubConv stubConv,			// (I) stub conv
		const SharedPointer<KRate> floatRate,		// (I) paid rate
		const char *formula,			// (I) paymt formula 
		const char *discZcName);		// (I) disc curve name

	/**
	 * General constructor with arbitrary formula schedule.
	 */
	KVPFloatLeg(
		const char *name,			// (I) object name
		const KVector(TDate)& resetDates,	// (I) reset dates
		const KVector(TDate)& resetEffDates,	// (I) reset eff dates
		const KVector(TDate)& accStartDates,	// (I) acc start dates
		const KVector(TDate)& accEndDates,	// (I) acc end dates
		const KVector(TDate)& payDates,		// (I) pay dates array
		const KVector(double)& notionals,	// (I) notionals array
		const KVector(String)& formulas,	// (I) pay formula array
		KDayCc dayCc,				// (I) pay day count
		KStubConv stubConv,			// (I) stub conv
		const SharedPointer<KRate> floatRate,		// (I) paid rate
		const char *discZcName);		// (I) disc curve name

	/**
	 * General constructor with arbitrary rate spread schedule.
	 * Since "floatRate" will be modified with spreads, it should
	 * not be declared as const structure.
	 */
	KVPFloatLeg(
		const char *name,			// (I) object name
		const KVector(TDate)& resetDates,	// (I) reset dates
		const KVector(TDate)& resetEffDates,	// (I) reset eff dates
		const KVector(TDate)& accStartDates,	// (I) acc start dates
		const KVector(TDate)& accEndDates,	// (I) acc end dates
		const KVector(TDate)& payDates,		// (I) pay dates array
		const KVector(double)& spreads,		// (I) spread array
		const KVector(double)& notionals,	// (I) notionals array
		KDayCc dayCc,				// (I) pay day count
		KStubConv stubConv,			// (I) stub conv
		const SharedPointer<KRate> floatRate,		// (I) paid rate
		const char *formula,			// (I) paymt formula 
		const char *discZcName);		// (I) disc curve name


	/**
	 * Convenience constructor.
	 * Creates a generic KVPFloatLeg from startDate, 
	 * payInterval, and maturityDate. If there
	 * is an implied stub, is assumed to be a front stub.
	 * The routine assumes  accrueStartDates=resetDates, 
	 * and accrueEndDates = payDates.
	 */
	KVPFloatLeg(
		const char *name,	// (I) name
		TDate startDate,	// (I) This date is not included
		TDate matDate,		// (I) maturity date
		const KDateInterval &freq,// (I) frequency as interval
		KDayCc dayCc,		// (I) 
		KStubConv stubConv,	// (I) 
		TBoolean stubAtEnd,	// (I) F=front, T=back
		const SharedPointer<KRate> floatRate,	// (I) 
		const char *formula,	// (I)
		const char *discZcName);// (I) 


	/**
	 * General constructor, arbitrary everything.
	 */
	KVPFloatLeg(
		const char *name,			// (I) object name
		const KVector(TDate)& resetDates,	// (I) reset dates
		const KVector(TDate)& resetEffDates,	// (I) reset eff dates
		const KVector(TDate)& accStartDates,	// (I) acc start dates
		const KVector(TDate)& accEndDates,	// (I) acc end dates
		const KVector(TDate)& payDates,		// (I) pay dates array
		const KVector(double)& notionals,	// (I) notionals array
		const KVector(String)& formulas,	// (I) pay formula array
		const KVector(SharedPointer<KRate>) &floatRates, // (I) flt rate array
		KDayCc dayCc,				// (I) pay day count
		KStubConv stubConv,			// (I) stub conv
		const char *discZcName);		// (I) disc curve name


	/**
	 * Destructor.
	 */
	~KVPFloatLeg();

	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPFloatLeg");};

    /**
     * Check validity.
     */
virtual void    CheckValid() const;

	/**
	 * Reads the object from a stream.
	 */
virtual	istream& Get(istream& is, int drw=FALSE);

	/**
	 * Writes the object to a stream.
	 */
virtual	ostream& Put(ostream& os, int indent = FALSE) const;

	/**
	 * Writes recursively the dependents. In this case,
	 * only single rate index plus constant spreads or step-up
	 * coupons are supported.
	 */
virtual ostream& YacctionWriteRecursive(ostream& os);

	/**
	 * Writes in yacction format.
	 */
virtual ostream& YacctionWrite( ostream& os, int indent=FALSE);

	/**
	 * Reads from a stream.
	 */
friend	istream& operator>>(istream& is, KVPFloatLeg& asset)
		{ asset.Get(is); return (is);}

	/**
	 * Writes to a stream.
	 */
friend	ostream& operator<<(ostream& os, const KVPFloatLeg& asset)
		{ asset.Put(os); return (os);}


	/**
	 * Default constructor.
	 */
	KVPFloatLeg(const char *name = NULL)
		: KVPInstr(name)
	{};

//protected:
public:

				/** Array of reset dates. */
	KVector(TDate)	mResetDates;
				/** Array of reset effective dates. */
	KVector(TDate)	mResetEffDates;
				/** Array of acc start dates. */
	KVector(TDate)	mAccStartDates;
				/** Array of acc end dates. */
	KVector(TDate)	mAccEndDates;
				/** Array of pay dates. */
	KVector(TDate)	mPayDates;
				/** Array of notionals. */
	KVector(double)	mNotionals;
				/** Array of payment formulas */
	KVector(String)	mFormulas;
				/** Day count convention. */
	KDayCc		mDayCc;
				/** Stub convention. */
	KStubConv	mStubConv;
				/** Rate reset definitions. */
	KVector(KRateReset*)	mRateResets;


friend	class	KVPToolFloatLeg;
friend	class	KVPToolFloatLeg2Idx;
friend	class	KVPToolFloatLeg3Idx;
};






#endif	/* _vpfleg_H */

