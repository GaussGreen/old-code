/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vtloption.h
 * Function:	
 * Author:	Christian Daher, David Liu
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_vtloption_H
#define	_vtloption_H

#include "vtlbase.h"
#include "vpoption.h"


//--------------------------------------------------------------
/**
 * Class for a weighted vector of assets.
 */

class KVPToolOption : public KVPToolAtom {
public:
	/**
	 * Constructor.
	 */
	KVPToolOption(SharedPointer<KVPOption> ins, KVTree &vt);

	/**
	 * Destructor.
	 */
virtual	~KVPToolOption();


	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPToolOption");};

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
		SharedPointerConvertTo(mOption, vp);
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
	 * "XTPROB" (probability of exercise), "XTEXP" (expected exercise time)
	 * "XTSDV" (standard dev of exercise time) and "XTFUG" (fugit).
	 */
virtual	KMap(String,double)  GetResults()
	{
		return mResults;
	}


private:
				/** Pointer to inetrument. */
	SharedPointer<KVPOption>	mOption;
				/** Curve name for slice geometry. */
	String		mCurveName;
				/** Current Value */
	KTSlice		*mValue;

				/** Claim bank for settlements */
	KMap(TDate,KTSlice*)	mExer;
				/** Default option value.
				  * This is the default value assigned to
				  * the option before any exercise event
				  * encountered, where mValue=NULL.
				  * To avoid return error in GetValue.
				  */
	KTSlice		*mDefValue;

				/** Results (PV,UNDER,FUGIT,etc.) */
	KMap(String,double)	mResults;

				/** Calculate fugit, etc. */
	int		mRunStat;
				/** Probabilty of exercise */
	KTSlice		*mXT0;
				/** Expected exercise time */
	KTSlice		*mXT1;
				/** Expected exercise time^2 */
	KTSlice		*mXT2;

};




#endif




