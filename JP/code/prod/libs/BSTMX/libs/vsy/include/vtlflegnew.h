/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtlflegnew.h
 * Function:	
 * Author:	Changhong He
 *****************************************************************/
#ifndef	_vtlflegnew_H
#define	_vtlflegnew_H

#include "vtlbase.h"
#include "vpflegnew.h"


//--------------------------------------------------------------
/**
 * Class for a weighted vector of assets.
 */

class KVPToolFloatLegNew : public KVPToolAtom {
public:
	/**
	 * Default Constructor
	 */
	KVPToolFloatLegNew(KVTree &vt) : KVPToolAtom(vt) {};

	/**
	 * Constructor.
	 */
	KVPToolFloatLegNew(
		SharedPointer<KVPFloatLegNew> ins,	// (I) inetrument to value
		KVTree &vt)			// (I) virtual tree
	: KVPToolAtom(vt) 
	{CreateVPToolFloatLegNew (ins, vt);}

	/**
	 * Destructor.
	 */
virtual	~KVPToolFloatLegNew();

	/**
	 * Constructor.
	 */
void	CreateVPToolFloatLegNew(
		SharedPointer<KVPFloatLegNew> ins,	// (I) inetrument to value
		KVTree &vt);			// (I) virtual tree

	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPToolFloatLegNew");};

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

void setValueDate( KVTree& vt);

void setCurveName();

protected:
				/** Pointer to instrument. */
	SharedPointer<KVPFloatLegNew>	mFloatLeg;
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




