/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vpprotleg.h
 * Function:	
 * Author:	David Liu
 *          Modified by Charles Morcom to handle varying 
 *          notional
 *****************************************************************/
#ifndef	_vpprotfleg_H
#define	_vpprotfleg_H

#include "vpbase.h"

/**
 * Pay upon default (PAY_DEF) or at maturity (PAY_MAT)
 */
enum KProtPayConv {
    PAY_DEF,         // pay immediate upon default 
    PAY_MAT          // pay at a fixed date, typically at maturity
};



//--------------------------------------------------------------
/**
*  This structure defines a CDS protection leg.
*/


class KVPProtLeg : public KVPInstr {
public:

	/**
	 * General constructor 
	 */
	KVPProtLeg(
		const char     *name,		// (I) object name
		const TDate     stDate,	    // (I) start date
		const TDate     endDate,	// (I) end date
		const double    notional,	// (I) notional
		const char     *recovery,	// (I) recovery rate
        KProtPayConv    payType,    // (I) default payment type
		const char     *discZcName);// (I) CDS curve name

    /**
     * Varying notionals
     */
    KVPProtLeg (
        const char*             name,		// (I) object name
	    const KVector(TDate)&   startDates,	// (I) Notional start dates
        const KVector(TDate)&   endDates,   // (I) Notional period end dates
	    const KVector(double)&  notionals,	// (I) notional amounts
	    const char*             recovery,	// (I) recovery rate
        KProtPayConv            payType,    // (I) default payment type
	    const char*             discZcName  // (I) CDS curve name
    );

	/**
	 * Destructor.
	 */
	~KVPProtLeg() {};

	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPProtLeg");};

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
	 * Writes recursively the dependents. 
	 */
virtual ostream& YacctionWriteRecursive(ostream& os);

	/**
	 * Writes in yacction format.
	 */
virtual ostream& YacctionWrite( ostream& os, int indent=FALSE);

	/**
	 * Reads from a stream.
	 */
friend	istream& operator>>(istream& is, KVPProtLeg& asset)
		{ asset.Get(is); return (is);}

	/**
	 * Writes to a stream.
	 */
friend	ostream& operator<<(ostream& os, const KVPProtLeg& asset)
		{ asset.Put(os); return (os);}


	/**
	 * Default constructor.
	 */
	KVPProtLeg(const char *name = NULL)
		: KVPInstr(name)
	{};

protected:
    /**Protection period start dates*/
    KVector(TDate)  mNtlStartDts;
    /**Protection period end dates*/
    KVector(TDate)  mNtlEndDts;
    /**Notional protected amount*/
    KVector(double) mNtlAmts;
                /** Recovery  */
    double          mRecovery;
                /** Default payment type */
    KProtPayConv    mPayType;
                /** Use default/binary recovery rate */
    bool            mIsDefRecovery;

friend	class	KVPToolProtLeg;
};


#endif	/* _vpprotleg_H */

