/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtldefprotect.h
 * Function:	
 * Author:	David Liu
 *****************************************************************/
#ifndef	_vtldefprotect_H
#define	_vtldefprotect_H

#include "vtlbase.h"
#include "vpdefprotect.h"


//--------------------------------------------------------------
/**
 * Class for a defprotect component.
 */

class KVPToolDefProtect : public KVPToolAtom {
public:
	/**
	 * Constructor.
	 */
	KVPToolDefProtect(SharedPointer<KVPDefProtect> ins, KVTree &vt);

	/**
	 * Destructor.
	 */
virtual	~KVPToolDefProtect();


	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPToolDefProtect");};

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
		SharedPointerConvertTo(mDefProtect, vp);
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
				/** Pointer to inetrument. */
	SharedPointer<KVPDefProtect>	mDefProtect;
				/** Curve name for slice geometry. */
	String		        mCurveName;
				/** Current Value */
	KTSlice		        *mValue;

				/** Claim bank for underlying
                                  * Use idx instead of settleDate
                                  * to avoid ambiguity if some settleDates are
                                  * the same.
                                  * (Should change option and knockIO as well!!)
                                  */
	KMap(int,KTSlice*)	mUnder;

				/** Temp slice for default prob*/
	KTSlice                 *mTmpTs;

				/** Default defprotect value.
				  * This is the default value assigned to
				  * the defprotect before any defprotect period
				  * is encountered, where mValue=NULL.
				  * To avoid return error in GetValue.
				  */
	KTSlice		        *mDefValue;

                                /** Array of default zero resets. */
        KVector(KZeroReset*)    mDefZeros;

                                /** Value date of risky curve 
                                  * (defprotect base date) 
                                  */
        TDate                   mValueDate;

};


#endif




