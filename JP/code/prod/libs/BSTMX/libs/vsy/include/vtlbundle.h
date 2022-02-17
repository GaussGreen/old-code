/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtlbundle.h
 * Function:	
 * Author:	Christian Daher, David Liu
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_vtlbundle_H
#define	_vtlbundle_H

#include "vpbundle.h"
#include "vtlbase.h"

//--------------------------------------------------------------
/**
 * Class for a weighted vector of assets.
 */

class KVPToolWBundle : public KVPToolAtom {
public:
	/**
	 * Constructor.
	 */
	KVPToolWBundle(SharedPointer<KVPWBundle> ins, KVTree &vt);

	/**
	 * Destructor.
	 */
	~KVPToolWBundle();


	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPToolWBundle");};

	/**
	 * Returns the name for the geometry of the slice
	 * (see description of method in KVPToolAtom).
	 */
virtual	const String&	GetCurveName();

	/**
	 * Returns the KVPAtom to which the KVPToolAtom is associated.
	 */
virtual	SharedPointer<KVPAtom>	Atom()
	{
		SharedPointer<KVPAtom> vp;
		SharedPointerConvertTo(mWBundle, vp);
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



protected:
				/** Curve name used for the slice geometry. */
	String				mCurveName;

				/** WBundle associated. */
	SharedPointer<KVPWBundle>	mWBundle;

				/** Slice value. */
	KTSlice				*mValue;

};





#endif // _vtlbundle_H




