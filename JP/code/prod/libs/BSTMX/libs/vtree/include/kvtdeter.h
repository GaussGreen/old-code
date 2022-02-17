/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	kvtdeter.h
 * Function:	
 * Author:	
 * Revision:	$Header$
 ***************************************************************/
#ifndef	_kvtdeter_H
#define	_kvtdeter_H
#include "kvtree.h"	// Standard definitions & error hadling

extern "C" {
#include "modlcons.h"                   /* GTO_MODEL_XXX constants */
#include "modltype.h"
#include "zerodate.h"                   /* GtoZeroDates  */
#include "virttree.h"                   /* TVirtualTree */

};


//--------------------------------------------------------------
/**
 * Class definition for pure deterministic model.
 */

class	KVTreeDeter : public KVTree {
public:


	/** Virtual method of KVTree. */
virtual			~KVTreeDeter();


	/** Virtual method of KVTree. */
virtual	void		Insert(TDate critDate);

	/** Virtual method of KVTree. */
virtual void		Insert(const KZeroReset&, bool isCrit=TRUE);

	/** Virtual method of KVTree. */
virtual TDate   	Insert(const KRateReset&, bool isCrit=TRUE);

	/** Virtual method of KVTree.  */
virtual	void	Get(KTSlice& ts, const KZeroReset&);

	/** Virtual method of KVTree. */
virtual	void	Get(KTSlice& ts, const KRateReset&);



	/** Virtual method of KVTree. */
virtual	void		Calibrate();

	/** Virtual method of KVTree. */
virtual void		Update(int tpIdx);



	/** Virtual method of KVTree. */
virtual KTSlice&	TSliceDev(KTSlice&, const String&);

	/** Virtual method of KVTree. */
virtual KTSlice&	TSliceEv(KTSlice&);

	/** Virtual method of KVTree. */
virtual	KTSlice&   	TSliceCreate(KTSlice&);

	/** Virtual method of KVTree. */
virtual	void		TSliceDestroy(KTSlice&);

	/** Virtual method of KVTree. */
virtual	double		TSliceGetCenter(KTSlice&);

	/** Virtual method of KVTree. */
virtual	bool		TSliceCompare(KTSlice&, KTSlice&, KTSComp);

	/** Virtual method of KVTree. */
virtual	KTSlice&	TSliceScalarOper(KTSlice&, double, KOper);

	/** Virtual method of KVTree. */
virtual	KTSlice&	TSliceUnaryOper(KTSlice&, const KTSlice&, KOper);

	/** Virtual method of KVTree. */
virtual	void		TSlicePut(KTSlice&, ostream& os, int minMaxOnly);

	/** Virtual method of KVTree. */
virtual	void		TSliceSpecialOper(KTSlice&, char* what, ...);

	/** Virtual method of KVTree. */
//virtual KTSlice&	TSliceInterp(const  KTSlice& ts1, 
//				     const  KTSlice& ts2, 
//				     double w1,
//				     double w2,
//				     int    interpType);



	/**
	 * Create and set the begining of slice node 
	 */
virtual	KTSliceNode&	TSliceNodeBegin(KTSlice&, KTSliceNode&);

	/**
	 * Test the end of of slice node 
	 */
virtual	bool		TSliceNodeEnd(KTSliceNode&);

	/**
	 * Move to next slice node 
	 */
virtual	KTSliceNode&	TSliceNodeNext(KTSliceNode&);

	/**
	 * Get the value of slice node 
	 */
virtual	double&		TSliceAccessNodeValue(KTSlice&, KTSliceNode&);

	/**
	 *  Get the maximum difference of node values of slice 
	 *  in nearest neighbors of specified node 
	 */
virtual	double		TSliceGetNodeStepMax(KTSlice&, KTSliceNode&);

	/** Prints a slice node on a stream. */
virtual	void		TSliceNodePut(KTSliceNode&, ostream& os);




	/** Virtual method of KVTree. */
virtual	int		TPNum();

	/** Virtual method of KVTree. */
virtual	int		TPIdxCurrent();

	/** Returns the date corresponding to the current timepoint. */
virtual	TDate		TPDateCurrent();


	/** Returns the year fraction from the start of the tree 
	 *  to the current timepoint. */
virtual	double		TPTimeCurrent();

	/** Virtual method of KVTree. */
virtual int             TPIdx(TDate date);

	/** Virtual method of KVTree. */
virtual	TDate		TPDate(int tpIdx);

	/** Returns the today date. */
virtual	TDate		TPToday()
	{
		return mToday;
	}

	/** Returns the value date. */
virtual	TDate		TPValueDate()
	{
		return mZcCurves[0]->fBaseDate;
	}


	/** Virtual method of KVTree. */
virtual	int		TPIsCriticalDate(int tpIdx);




	/** Constructor */
			KVTreeDeter();

	/**
	 * Initialisation.
	 *
	 * valueCurveName: curve to compute PV as value date
	 */
	void	Initialize(
		TDate todayDate,	   // (I) 
		int   numZCurves,	   // (I) 
		TCurve ** zcCurves,        // (I) [numZCurves] first = DISCOUNT 
		KVector(String) curveNames,// (I) [numZCurves]
		const String &valueCurveName);// (I) reference disc curve


	/**
	 * Returns the discount factor between today and the value
	 * date of the reference discount curve value date.
	 */
	double	DiscZeroShift()
			{ return mValueDiscZeroShift; }

private:
	double	FwdDisc(
		TDate date1,			// (I) 
		TDate date2,			// (I) 
		const String& discCurveName);	// (I) 

	TDate		mToday;
	int             mNumZcCurves;
	TCurve **       mZcCurves;       // [mNumZcCs] mZcCurves[0]=DISCOUNT
	int		mTpIdx;		

	TDateList	*criticalDL;

	String		mValueCurveName;	// reference disc curve
	double		mValueDiscZeroShift;
};













#endif /* _kvtdeter_H */
