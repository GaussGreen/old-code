/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtlcplxrate.h
 * Function:	
 * Author:	Changhong He
 *****************************************************************/
#ifndef	_vtlindicator_H
#define	_vtlindicator_H

#include "kindicator.h"
#include "vtlbase.h"

const String INRANGE  = "INSIDE_RANGE";
const String OUTRANGE = "OUTSIDE_RANGE";
const String SMOOTHNO = "SMOOTH_NONE";
const String SMOOTHS  = "SMOOTH_SINGLE";
const String SMOOTHD  = "SMOOTH_DOUBLE";

//--------------------------------------------------------------
/**
 * Class for a weighted vector of assets.
 */

class KVPToolIndicator : public KVPToolAtom {
public:
	/**
	 * Constructor.
	 */
	KVPToolIndicator(SharedPointer<KIndicator> ins, KVTree &vt);

	/**
	 * Destructor.
	 */
	~KVPToolIndicator();


	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPToolIndicator");};

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
		SharedPointerConvertTo(mIndicator, vp);
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


    void    setFormula( double    LoBarrier,
                        double    HiBarrier,
                        KKnockIO  IoWindow,
                        KSmooth   smooth);

    String& getFormula () 
        { return mFormula;   }

protected:
				/** Curve name used for the slice geometry. */
	String				mCurveName;

    String              mFormula;

				/** WBundle associated. */
	SharedPointer<KIndicator>		mIndicator;

				/** Slice value. */
	KTSlice				*mValue;

};





#endif // _vtlindicator_H



