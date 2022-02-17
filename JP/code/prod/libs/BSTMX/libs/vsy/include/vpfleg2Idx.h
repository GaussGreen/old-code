/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vpfleg.h
 * Function:	
 * Author:	David Liu
 *****************************************************************/
#ifndef	_vpfleg2Idx_H
#define	_vpfleg2Idx_H

#include "vpbase.h"
#include "krate.h"
#include "vpfleg.h"



//--------------------------------------------------------------
/**
 * This class describes a payment leg, i.e. a sequence of payments
 * based on combination of two floating indices with arbitrary
 * payoff formula.
 */


class KVPFloatLeg2Idx : public KVPFloatLeg {
public:

	/**
	 * General constructor.
	 */
	KVPFloatLeg2Idx(
		const char *name,			// (I) object name
		const KVector(TDate)& resetDates,	// (I) reset dates
		const KVector(TDate)& resetEffDates,	// (I) reset eff dates
		const KVector(TDate)& accStartDates,	// (I) acc start dates
		const KVector(TDate)& accEndDates,	// (I) acc end dates
		const KVector(TDate)& payDates,		// (I) pay dates array
		const KVector(double)& notionals,	// (I) notionals array
		KDayCc dayCc,				// (I) pay day count
		KStubConv stubConv,			// (I) stub conv
		const SharedPointer<KRate> floatRate1,	// (I) 1st rate index
		const SharedPointer<KRate> floatRate2,	// (I) 2nd rate index
		const char *formula,			// (I) paymt formula 
		const char *discZcName);		// (I) disc curve name

	/**
	 * General constructor with schedule of payment formula.
	 */
	KVPFloatLeg2Idx(
		const char *name,			// (I) object name
		const KVector(TDate)& resetDates,	// (I) reset dates
		const KVector(TDate)& resetEffDates,	// (I) reset eff dates
		const KVector(TDate)& accStartDates,	// (I) acc start dates
		const KVector(TDate)& accEndDates,	// (I) acc end dates
		const KVector(TDate)& payDates,		// (I) pay dates array
		const KVector(double)& notionals,	// (I) notionals array
		const KVector(String)&formula,		// (I) paymt formula arr
		KDayCc dayCc,				// (I) pay day count
		KStubConv stubConv,			// (I) stub conv
		const SharedPointer<KRate> floatRate1,	// (I) 1st rate index
		const SharedPointer<KRate> floatRate2,	// (I) 2nd rate index
		const char *discZcName);		// (I) disc curve name

	/**
	 * General constructor with schedule of payment formula and two rates.
	 */
	KVPFloatLeg2Idx(
		const char *name,			// (I) object name
		const KVector(TDate)& resetDates,	// (I) reset dates
		const KVector(TDate)& resetEffDates,	// (I) reset eff dates
		const KVector(TDate)& accStartDates,	// (I) acc start dates
		const KVector(TDate)& accEndDates,	// (I) acc end dates
		const KVector(TDate)& payDates,		// (I) pay dates array
		const KVector(double)& notionals,	// (I) notionals array
		const KVector(String)& formula,		// (I) paymt formula arr
		const KVector(SharedPointer<KRate>) &floatRate1, // (I) 1st rate index
		const KVector(SharedPointer<KRate>) &floatRate2, // (I) 2nd rate index
		KDayCc dayCc,				// (I) pay day count
		KStubConv stubConv,			// (I) stub conv
		const char *discZcName);		// (I) disc curve name


	/**
	 * Convenience constructor.
	 * Creates a generic KVPFloatLeg from startDate, 
	 * payInterval, and maturityDate. If there
	 * is an implied stub, is assumed to be a front stub.
	 * The routine assumes  accrueStartDates=resetDates, 
	 * and accrueEndDates = payDates.
	 */
	KVPFloatLeg2Idx(
		const char *name,		// (I) name
		TDate startDate,		// (I) This date is not included
		TDate matDate,			// (I) maturity date
		const KDateInterval &freq,	// (I) frequency as interval
		KDayCc dayCc,			// (I) 
		KStubConv stubConv,		// (I) 
		TBoolean stubAtEnd,		// (I) F=front, T=back
		const SharedPointer<KRate> floatRate1,// (I) 1st rate index
		const SharedPointer<KRate> floatRate2,// (I) 2nd rate index
		const char *formula,		// (I)
		const char *discZcName);	// (I) 

	/**
	 * Destructor.
	 */
	~KVPFloatLeg2Idx();

	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPFloatLeg2Idx");};

    /**
     * Check validity.
     */
virtual void    CheckValid() const;

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
	 * Writes to a stream.
	 */
friend	ostream& operator<<(ostream& os, const KVPFloatLeg2Idx& asset)
		{ asset.Put(os); return (os);}


	/**
	 * Default constructor.
	 */
	KVPFloatLeg2Idx(const char *name = NULL)
		: KVPFloatLeg(name)
	{};

private:

				/** 2nd index rate reset definitions. */
	KVector(KRateReset*)	mRateResets2;

friend	class	KVPToolFloatLeg2Idx;
};































#endif	/* _vpfleg_H */

