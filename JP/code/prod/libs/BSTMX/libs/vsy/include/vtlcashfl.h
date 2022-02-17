/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtlcashfl.h
 * Function:	
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_vtlcashfl_H
#define	_vtlcashfl_H

#include "vtlbase.h"
#include "vpcashfl.h"


//--------------------------------------------------------------
/**
 * Class for a weighted vector of assets.
 */

class KVPToolCashFlows : public KVPToolAtom {
public:
	/**
	 * Constructor.
	 */
	KVPToolCashFlows(
		SharedPointer<KVPCashFlows> ins, 
		KVTree &vt);

	/**
	 * Destructor.
	 */
virtual	~KVPToolCashFlows();


	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPToolCashFlows");};

	/**
	 * Returns the name for the geometry of the slice
	 * (see description of method in KVPToolAtom).
	 */
virtual	const String	&GetCurveName()
			{ return(mCashFlows->GetDiscName()); }


	/**
	 * Returns the KVPAtom to which the KVPToolAtom is associated.
	 */
virtual	SharedPointer<KVPAtom>	Atom()
	{
		SharedPointer<KVPAtom> vp;
		SharedPointerConvertTo(mCashFlows, vp);
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
				/** Pointer to inetrument. */
	SharedPointer<KVPCashFlows>	mCashFlows;
				/** Current Value */
	KTSlice		*mValue;
				/** Default Value */
	KTSlice		*mDefValue;
};




#endif





