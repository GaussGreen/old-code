/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vpcashfl.h
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_vpcashfl_H
#define	_vpcashfl_H

#include "vpbase.h"


//--------------------------------------------------------------
/**
 * @author Christian Daher
 * @version 2.0
 */


class KVPCashFlows : public KVPInstr {
public:


	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPCashFlows");};

	/**
	 * Reads the object from a stream.
	 */
virtual	istream& Get(istream& is, int drw=FALSE);

	/**
	 * Writes the object to a stream.
	 */
virtual	ostream& Put(ostream& os, int indent = FALSE) const;

	/**
	 * Writes in yacction format.
	 */
virtual ostream& YacctionWrite( ostream& os, int indent=FALSE);

	/**
	 * Reads from a stream.
	 */
friend	istream& operator>>(istream& is, KVPCashFlows& asset)
		{ asset.Get(is); return (is);}

	/**
	 * Writes to a stream.
	 */
friend	ostream& operator<<(ostream& os, const KVPCashFlows& asset)
		{ asset.Put(os); return (os);}


	/**
	 * Default Constructor.
	 */
	KVPCashFlows(const char *name = NULL)
		: KVPInstr(name)
	{}


	/**
	 * Constructor.
	 */
	KVPCashFlows(
		const char* name,
		const KVector(TDate) payDates,
		const KVector(double) amounts,
		const char* discZcName);


				/** Array of reset dates. */
	KVector(TDate)	mPayDates;
				/** Array of acc start dates. */
	KVector(double)	mPayAmounts;
};































#endif	/* _vpcashfl_H */

