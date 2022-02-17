/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtlkio.h
 * Function:	
 * Author:	David Liu
 *****************************************************************/
#ifndef	_vtlkio2Idx_H
#define	_vtlkio2Idx_H

#include "vtlbase.h"
#include "vpkio2idx.h"
#include "vtlkio.h"


//--------------------------------------------------------------
/**
 * Class for a dual knock-in/out component
 */

class KVPToolKnockIO2Idx : public KVPToolKnockIO {
public:
	/**
	 * Constructor.
	 */
	KVPToolKnockIO2Idx(SharedPointer<KVPKnockIO2Idx> ins, KVTree &vt);

	/**
	 * Destructor.
	 */
virtual	~KVPToolKnockIO2Idx();


	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPToolKnockIO2Idx");};

	/**
	 * Returns the KVPAtom to which the KVPToolAtom is associated.
	 */
virtual	SharedPointer<KVPAtom>	Atom()
	{
		SharedPointer<KVPAtom> vp;
		SharedPointerConvertTo(mKnockIO2Idx, vp);
		return (vp);
	}

	/**
	 * Update tool.
	 */
virtual	void	Update();



private:
				/** Pointer to instrument. */
	SharedPointer<KVPKnockIO2Idx>	mKnockIO2Idx;
				/** Curve name for slice geometry. */
	String		mCurve2Name;
				/** Curve name for trigger rate.   */
	String		mRateIdx2Name;
};


#endif




