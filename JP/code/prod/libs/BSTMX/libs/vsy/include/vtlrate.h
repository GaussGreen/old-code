/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtlrate.h
 * Function:	
 * Author:	David Liu
 *****************************************************************/
#ifndef	_vtlrate_H
#define	_vtlrate_H

#include "krate.h"
#include "vtlbase.h"

//--------------------------------------------------------------
/**
 * Class for a weighted vector of assets.
 */

class KVPToolRate : public KVPToolAtom {
public:
	/**
	 * Constructor.
	 */
	KVPToolRate(SharedPointer<KRate> ins, KVTree &vt);

	/**
	 * Destructor.
	 */
	~KVPToolRate();


	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPToolRate");};

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
		SharedPointerConvertTo(mRate, vp);
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
	SharedPointer<KRate>		mRate;

				/** Slice value. */
	KTSlice				*mValue;

};





#endif // _vtlrate_H




