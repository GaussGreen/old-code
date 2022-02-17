/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vpbase.h
 * Function:	
 * Author:	Christian Daher, David Liu
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_vpbase_H
#define	_vpbase_H

#include "kstdinc.h"
#include "ktypes.h"
#include "kvpatom.h"		// KVPAtom class


//--------------------------------------------------------------
/**
 * Base class for an instrument, i.e. a stream to be valued
 * in the tree. This class just contains a list pointing
 * to other instruments on which the current instrument
 * is assumed to depend (e.g. if the current instrument
 * is an swap option, it contains the underlying fixed and
 * floating legs as dependencies).
 */

class KVPInstr : public KVPAtom {
public:

	/**
	 * Default constructor
	 */
	KVPInstr() {}

	/**
	 * Default constructor
	 */
	KVPInstr(const char *name) : KVPAtom(name) {}

	/**
	 * Destructor.
	 */
virtual	~KVPInstr();
	


	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPInstr");};


	/**
	 * Reads the object from a stream.
	 */
virtual	istream& Get(istream& is, int drw=FALSE)
	{
		throw KFailure("KVPInstr::Get: NA.\n");
		return(is);
	}

	/**
	 * Reads from a stream.
	 */
friend	istream& operator>>(istream& is, KVPInstr& asset)
		{ asset.Get(is); return (is);}


	/**
	 * Sets the name of the discount curve.
	 */
	void    SetDiscName(const char *zcName);



        /**
	 * Gets discount curve name
	 */
virtual const String& GetDiscName() const
	{
		return mDiscZcName;
	}



protected:
				/** Discount zero curve name */
	String		mDiscZcName;

};




#endif


