/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtlkio.h
 * Function:	
 * Author:	David Liu
 *****************************************************************/
#ifndef	_vtlkio_H
#define	_vtlkio_H

#include "vtlbase.h"
#include "vpkio.h"


//--------------------------------------------------------------
/**
 * Class for a weighted vector of assets.
 */

class KVPToolKnockIO : public KVPToolAtom {
public:
	/**
	 * Constructor.
	 */
	KVPToolKnockIO(SharedPointer<KVPKnockIO> ins, KVTree &vt);

	/**
	 * Destructor.
	 */
virtual	~KVPToolKnockIO();


	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPToolKnockIO");};

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
		SharedPointerConvertTo(mKnockIO, vp);
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

	/**
	 * Calculates and returns the results at time.
	 * In addition to "PV", also contains the fields
	 * "XTPROB" (probability of hit), "XTEXP" (expected hit time)
	 * "XTSDV" (standard dev of hit time) and "XTFUG" (fugit).
	 */
virtual	KMap(String,double)  GetResults()
	{
		return mResults;
	}


private:
				/** Pointer to inetrument. */
	SharedPointer<KVPKnockIO>	mKnockIO;
				/** Curve name for slice geometry. */
	String		mCurveName;
				/** Curve name for trigger rate.  */
	String		mRateIdxName;

                /** Default value
                  * This is the default value assigned to
                  * the component, to avoid return error in GetValue.
                  */
    KTSlice     *mDefValue;

				/** Current Value Slice */
	KTSlice		*mValue;
				/** Unsmoothed current Value Slice */
	KTSlice		*mValueUS;
				/** Current rebate slice */
	KTSlice		*mRebateValue;
				/** Current unsmoothed rebate slice */
	KTSlice		*mRebateValueUS;
				/** Used to compute the current KIO value. */
	KTSlice		*mGetValue;


				/** Model observation reset dates */
	KVector(TDate)	mModelObsDates;
				/** Claim bank for underlying settlements */
	KMap(TDate, KTSlice*)	mExer;
				/** Claim bank for rebate at settlement dates */
	KMap(TDate, KTSlice*)	mRebate;

				/** Run the fugit calculations, etc */
	int			mRunStat;
				/** Probabilty of hit */
	KTSlice			*mXT0;
				/** Expected hit time */
	KTSlice			*mXT1;
				/** Expected hit time^2 */
	KTSlice			*mXT2;

				/** Results (PV,UNDER,FUGIT,ETC) */
	KMap(String,double)	mResults;

friend  class   KVPToolKnockIO2Idx;

};




#endif




