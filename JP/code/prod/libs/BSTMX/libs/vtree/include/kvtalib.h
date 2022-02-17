/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	kvtalib.h
 * Function:	
 * Author:	
 * Revision:	$Header$
 ***************************************************************/
#ifndef	_kvtalib_H
#define	_kvtalib_H
#include "kvtree.h"	// Standard definitions & error hadling

#include "kmodpar.h"	// 

extern "C" {
#include "modlcons.h"                   /* GTO_MODEL_XXX constants */
#include "modltype.h"
#include "zerodate.h"                   /* GtoZeroDates  */
#include "virttree.h"                   /* TVirtualTree */

};


//--------------------------------------------------------------
/**
 * Pure virtual base class definition for tree.
 */

class	KVTreeAL : public KVTree {
public:


	/** Virtual method of KVTree. */
virtual			~KVTreeAL();


	/** Virtual method of KVTree. */
virtual	void		Insert(TDate critDate);

	/** Virtual method of KVTree. */
virtual void		Insert(const KZeroReset&, bool isCrit=TRUE);

	/** Virtual method of KVTree. */
virtual TDate   	Insert(const KRateReset&, bool isCrit=TRUE);

	/** Virtual method of KVTree. */
virtual TDate   	Insert(const KRateReset&, TDate, bool isCrit=TRUE);

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
virtual	TDate	GetValueDate(int curveIdx)
	{
		return mZcCurves[curveIdx]->fBaseDate;
	}


	/** Virtual method of KVTree. */
virtual	int		TPIsCriticalDate(int tpIdx);




	/** Constructor */
			KVTreeAL();


	/**
	 * Initialisation.
	 *
	 */
	void	Initialize(
		TDate todayDate,	   // (I) 
		int   numZCurves,	   // (I) 
		TCurve ** zcCurves,	   // (I) [numZCurves] first = DISCOUNT 
		KVector(String) curveNames,// (I) [numZCurves]
		const String& diffuseCurve,// (I) one of the above names
		int mNumFact,	           // (I) number of factors
		const double *mBeta,	   // (I) model param
		const double *mAlpha,	   // (I) model param
		const double *mRho,	   // (I) model param
		int mPpy,		   // (I) 
		double smoothFact,	   // (I) 
		double mNumStdevCut,	   // (I) 
		double mQ1,		   // (I) 
		double mQ2,		   // (I) 
		double mFwdSh,		   // (I) 
	
					//     VOLATILITY DATA 
		int numVolDates,	// (I) # of volatility dates 
		const TDate *volDates,	// (I) array of vol dates 
		const double *volMat,	// (I) array of vol fwd maturities 
		const int *volFreq,	// (I) array of vol frequency 
		const double *volRates);// (I) array of vol 





private:

	TDateList	*criticalDL;
	TZeroDates	*zeroDates;	// Used to set up zero bank 

	TModelInitFunc	*modelInitFunc;	/* Initializes model */
	OCalibInfo	calibInfo;	/* Calib info; w/ vol curve*/
	TTimeLineInfo	*tlInfo;	/* Periods/year, etc */
	TCalibFreeFunc	*calibFreeFunc;	// Used to free calibInfo 

	//TDate		today;		/* Where diffusion starts */

	TVirtualTree	vt;		/* Init all fields to 0 */
	OZeroBank	zeroBank;	/* Provides zeros */
	int		mTpIdx;		


					// Tree Input Data
	TDate		mToday;
	int             mNumZcCurves;
	TCurve **       mZcCurves;       // [mNumZcCs] mZcCurves[0]=DISCOUNT
	int             mDiffuseIdx;     // GTO_CURVE_DISCOUNT, GTO_CURVE_INDEX1 etc

	int		mNumFact;
	double		*mBeta;
	double		*mAlpha;
	double		*mRho;
	int		mPpy;
	double		mSmoothFact;
	double		mNumStdevCut;
	double		mQ1;
	double		mQ2;

	int		mNumVolDates;
	TDate		*mVolDates;
	double		*mVolMat;
	int		*mVolFreq;
	double		*mVolRates;




};













#endif /* _kvtalib_H */
