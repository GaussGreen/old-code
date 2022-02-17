/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vpfleg3Idx.h
 * Function:	
 * Author:	RObert Mattingly
 *****************************************************************/
#ifndef	_vpfleg3Idx_H
#define	_vpfleg3Idx_H

#include "vpbase.h"
#include "krate.h"
#include "vpfleg.h"



//--------------------------------------------------------------
/**
 * This class describes a payment leg, i.e. a sequence of payments
 * based on combination of two floating indices with arbitrary
 * payoff formula.
 */


class KVPFloatLeg3Idx : public KVPFloatLeg {
public:

	/**
	 * General constructor with schedule of payment formula.
	 */
	KVPFloatLeg3Idx(
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
		const SharedPointer<KRate> floatRate3,	// (I) 3rd rate index
		const char *discZcName);		// (I) disc curve name

	/**
	 * Destructor.
	 */
	~KVPFloatLeg3Idx();

	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPFloatLeg3Idx");};

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
friend	ostream& operator<<(ostream& os, const KVPFloatLeg3Idx& asset)
		{ asset.Put(os); return (os);}


	/**
	 * Default constructor.
	 */
	KVPFloatLeg3Idx(const char *name = NULL)
		: KVPFloatLeg(name)
	{};

protected:

	KVector(KRateReset*)	mRateResets2;
	KVector(KRateReset*)	mRateResets3;

friend	class	KVPToolFloatLeg3Idx;
};


#endif	/* _vpfleg3Idx_H */

