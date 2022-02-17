/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtlfleg2Idx.h
 * Function:	
 * Author:	David Liu
 *****************************************************************/
#ifndef	_vtlfleg2Idx_H
#define	_vtlfleg2Idx_H

#include "vtlbase.h"
#include "vtlfleg.h"
#include "vpfleg2Idx.h"


//--------------------------------------------------------------
/**
 * Class for floating leg with 2 rate indices.
 */

class KVPToolFloatLeg2Idx : public KVPToolFloatLeg {
public:
	/**
	 * Constructor.
	 */
	KVPToolFloatLeg2Idx(
		SharedPointer<KVPFloatLeg2Idx> ins, // (I) inetrument to value
		KVTree &vt);			    // (I) virtual tree

	/**
	 * Destructor.
	 */
virtual	~KVPToolFloatLeg2Idx();


	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPToolFloatLeg2Idx");};

	/**
	 * Returns the name for the geometry of the slice
	 * (see description of method in KVPToolAtom).
	 */
virtual	const String&	GetCurveName() { return mCurveName;}


	/**
	 * Returns the KVPAtom to which the KVPToolAtom is associated.
	 */
virtual	SharedPointer<KVPAtom>	Atom()
	{
		SharedPointer<KVPAtom> vp;
		SharedPointerConvertTo(mFloatLeg2Idx, vp);
		return (vp);
	}

	/**
	 * Update tool.
	 */
virtual	void	Update();


	/**
	 * Returns the current value of the stream.
	 */
virtual	KTSlice&	GetValue();



private:
					/** Pointer to instrument. */
	SharedPointer<KVPFloatLeg2Idx>	mFloatLeg2Idx;
					/** Array of model reset dates 2. */
	KVector(TDate)		mModelResetDates2;
};




#endif




