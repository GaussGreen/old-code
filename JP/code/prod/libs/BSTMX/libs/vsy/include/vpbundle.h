/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vpbundle.h
 * Function:	
 * Author:	Christian Daher, David Liu
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_vpbundle_H
#define	_vpbundle_H

#include "vpbase.h"

//--------------------------------------------------------------
/**
 * Class for a weighted vector of assets.
 */

class KVPWBundle : public KVPInstr {
public:
	/**
	 * Constructor.
	 */
	KVPWBundle(const char* name)
		: KVPInstr(name) {}

	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPWBundle");};


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
	 * Adds another instrument "instr" to the list
	 * of dependencies of the instrument.
	 */
virtual	KVPWBundle&	AddDepWeight(
						const SharedPointer<KVPInstr> instr, 
						double weight);

				/** weights */
	KVector(double)	mWeights;

    /**
     * Gets discount curve name
     * Use the first dependent discount name (arbitrary, might affact
     * the ZeroShift value from today to value date.)
     */
virtual const String& GetDiscName() const;

};



#endif


