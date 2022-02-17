/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtlfleg.h
 * Function:	
 * Author:	David Liu
 *****************************************************************/
#ifndef	_vtlfleg_H
#define	_vtlfleg_H

#include "vtlbase.h"
#include "vpfleg.h"


//--------------------------------------------------------------
/**
 * Class for a weighted vector of assets.
 */

class KVPToolFloatLeg : public KVPToolAtom {
public:
	/**
	 * Default Constructor
	 */
	KVPToolFloatLeg(KVTree &vt) : KVPToolAtom(vt) {};

	/**
	 * Constructor.
	 */
	KVPToolFloatLeg(
		SharedPointer<KVPFloatLeg> ins,	// (I) inetrument to value
		KVTree &vt)			// (I) virtual tree
	: KVPToolAtom(vt) 
	{CreateVPToolFloatLeg (ins, vt);}

	/**
	 * Destructor.
	 */
virtual	~KVPToolFloatLeg();

	/**
	 * Constructor.
	 */
void	CreateVPToolFloatLeg(
		SharedPointer<KVPFloatLeg> ins,	// (I) inetrument to value
		KVTree &vt);			// (I) virtual tree

	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPToolFloatLeg");};

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
		SharedPointerConvertTo(mFloatLeg, vp);
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
				/** Pointer to instrument. */
	SharedPointer<KVPFloatLeg>	mFloatLeg;
                     /** Array of model reset dates. */
    KVector(TDate)      mModelResetDates;
				 /** Array of full accrual dates. */
    KVector(TDate)      mFullAccDates;

				/** Array of IR zero reset definition. */
	KVector(KZeroReset*)	mIRZeroResets;
				/** Array of survival zero reset to accStartDate definition. */
	KVector(KZeroReset*)	mDefZeroAccStarts;
				/** Array of survival zero reset to accEndDate definition. */
	KVector(KZeroReset*)	mDefZeroAccEnds;


                    /** Is default accrual. True is both discount
                     *  curve is risky and mDayCc is on for default
                     *  accrual
                     */
    bool            mDefAccrual;

                    /** Default value
                      * This is the default value assigned to
                      * the component, to avoid return error in GetValue.
                      */
    KTSlice         *mDefValue;
					/** Unadjusted value */
	KTSlice			*mUnValue;
					/** Current clean value */
	KTSlice			*mClValue;
					/** Tmp slice used to store the pynt.*/
	KTSlice			*mTmpTs;
					/** Tmp slice used to store the rate.*/
	KTSlice			*mFRateTs;
					/** Pending resets. */
	KMap(int,KTSlice*)	mPResets;
					/** name for geometry of the slice */
	String			mCurveName;
					/** Value date of leg */
	TDate			mValueDate;

};




#endif




