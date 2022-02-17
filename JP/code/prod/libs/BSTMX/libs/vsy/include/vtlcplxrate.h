/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtlcplxrate.h
 * Function:	
 * Author:	Changhong He
 *****************************************************************/
#ifndef	_vtlcplxrate_H
#define	_vtlcplxrate_H

#include "kcplxrate.h"
#include "vtlbase.h"

//--------------------------------------------------------------
/**
 * Class for a weighted vector of assets.
 */

class KVPToolCplxRate : public KVPToolAtom {
public:
	/**
	 * Constructor.
	 */
	KVPToolCplxRate(SharedPointer<KCplxRate> ins, KVTree &vt);

	/**
	 * Destructor.
	 */
	~KVPToolCplxRate();


	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPToolCplxRate");};

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
		SharedPointerConvertTo(mCplxRate, vp);
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
	SharedPointer<KCplxRate>		mCplxRate;

				/** Slice value. */
	KTSlice				*mValue;

};





#endif // _vtlcplxrate_H



