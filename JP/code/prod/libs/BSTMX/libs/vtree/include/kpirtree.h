/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	kpirtree.h
 * Function:	
 * Author:	David Liu
 ***************************************************************/
#ifndef	_kpirtree_H
#define	_kpirtree_H
#include "kstdinc.h"    // Standard definitions & error hadling 
#include "kmodpar.h"    // class KMrParam, kSmileParam
#include "kmktcrv.h"	// class KMarketCurves
#include "kvoldat.h"	// class KVolDiag

#include "kvtree.h"
#include "kmrntree.h"
#include "kzcurve.h"
#include "kzrbank.h"
#include "krstbank.h"	// Past rate resets information (KResetBank)


extern  "C" {
#include "datelist.h"
#include "date_sup.h"
 
#include "modltype.h"
#include "ir1fact.h"
#include "ir1fini.h"                    /* Gto1FIRModelInit */
#ifdef UNIX
#include "ir2fact.h"                    /* Gto2FIRCalibInfoNewWrap */
#include "ir2fini.h"                    /* Gto2FIRModelInit */
#endif
#include "modlcons.h"                   /* GTO_MODEL_XXX constants */
#include "calibir.h"                    /* TCalibInfoIR  */
 
#include "parfixed.h"                    /* TCalibInfoIR  */
 
#include "zerodate.h"                   /* TZeroDates, GtoZeroDatesBuild1 */ 
 
#include "drlio.h"              	/* FScanStruct */
#include "drltime.h"
#include "drlts.h"
#include "vnfmanly.h"
 
};

	class	KCetData;


//
// Accuracy on sum of state prices.
//
const   double  STATE_PRICE_TOL = 1e-4;


/**
 * A pure interest rate tree class based
 * on an N-factor mean reverting tree.
 */

class KPirTree : public KMrNTree {
public:
	//------------------------------------------------------
	// KVTree base class methods
	//------------------------------------------------------

	/** Default constructor. */
			KPirTree();

	/** Destructor. */
virtual			~KPirTree();

	/**
         * Insert extra zero bank dates, and then
	 * compute the tree timeline according to the ppy rule.
         */
virtual void		SetUpTimeline();


	/**
	 * Calibrate should be called after SetUpTimeline(),
	 * and it calls successively
	 *	TreeSetUp,
	 *	CheckTreeValid and
	 * 	CalibrateDrift.
	 *
	 * It does the following:
	 * 1. Call KMrNTree::Calibrate to compute vol dependent tree 
	 *    parameters such as jump size, orthogonal factors, 
	 *    and tree limits, etc.)
	 * 2. Compute the zero prices and forward rates at each time step.
	 * 3. Compute the the drift at each time step in the forward loop.
	 */
virtual	void		Calibrate();

	/** Update IR Tree. */
virtual	void		Update(int tpIdx);

	/** DEV a slice in IR Tree. */
virtual KTSlice&	TSliceDev(KTSlice&, const String&);

	/** Forward a slice. */
virtual KTSlice&	TSliceFwd(KTSlice&);

	/** Insert a critical date in the tree */
virtual void		Insert(TDate critDate)
			{ KMrNTree::Insert(critDate);}	

	/** Insert zero dates in the tree and ZBank */
virtual void		Insert(const KZeroReset &zeroReset, bool isCrit);

	/** Insert a KRateReset, return the rate reset date. */
virtual TDate		Insert(const KRateReset &rtReset, bool isCrit);

	/** Insert a KCplxRateReset, return the rate reset date. */
virtual TDate		Insert(const KCplxRateReset &rtCplxReset, bool isCrit);

	/** Insert a KRateReset between reset date and end date,
	 *  and return the rate reset date. 
	 */
virtual TDate	Insert(const KRateReset &rtReset, TDate endDate,  bool isCrit);

virtual TDate	Insert(const KCplxRateReset &rtReset, TDate endDate,  bool isCrit);

	/** Create a slice of dimension mIRDim */
virtual	KTSlice&	TSliceCreate(KTSlice&);

	/**
	 * Calculates the discount factor applicable between
	 * tpIdx and tpIdx+1.
	 */
virtual	void		CalcDiscount(
			int tpIdx,		// (I) time point index
			int nIRDim,		// (I) IR dimemsions
			double Zt,		// (I) center offset 
			KMap(int,double*) &mDiscountIR);// (O) discount slice

	/**
	 * Solves for the incremental offset to calibrate
	 * the zero price at tpIdx+1.
	 */
virtual	void		SolveOffset(
			int tpIdx,		// (I) time point index
			int nIRDim,		// (I) IR dimemsion
			double *discount,	// (I) slice disc bet t and t+1
			double *statePr,	// (I) slice state price at t
			double zeroPrice,	// (I) zero price at t+1
			double Zt,		// (I) intial offset
			double *DelZt);		// (O) incremental offset 


	/** Compute a zero reset.  The pure interest tree can
	 *  deal with two types of zeros:
	 *  1. rate out of diffuse curve.
	 *  2. rate with deterministic spread over diffuse curve.
	 */
virtual	void		Get(KTSlice&, const KZeroReset&);

	/** Compute a reset index rate  
	 *  The pure interest tree can deal with following 
	 *  two types of rate: 
	 *  1. rate out of diffuse curve.
	 *  2. rate with deterministic spread over diffuse curve.
	 */
virtual	void		Get(KTSlice&, const KRateReset&);

virtual	void		Get(KTSlice&, const KCplxRateReset&);

	/** Check tree validity */
virtual void		CheckTreeValid();


	//------------------------------------------------------
	// KPirTree specific public methods
	//------------------------------------------------------

	/** Get the KZCurve given a curve index */
virtual TDate		GetValueDate(int curveIdx); 

	/** Get the KZCurve given a curve index */
virtual KZCurve&	GetKZCurve(int	curveIdx); 

	/** Get the ZeroBank given a curve index */
virtual KZeroBank&	GetZeroBank(int curveIdx); 

	/** Get the discount factor at tpIdx, given a curve index */
virtual double*		GetDiscount(int curveIdx); 

	/** Get the zero coupon price on a date given a curve index */
virtual double		GetZeroPrice(int curveIdx, TDate dt); 

	/** Get the zero coupon price given a curve index */
virtual double*		GetZeroPrice(int curveIdx); 

	/** Get the forward rate given a curve index */
virtual double*		GetForwardRate(int curveIdx); 


	/** Compute the backbone factor at given time point */
virtual double		GetVolBbq(int t); 



	/**
	 * Initialize performs:<br>
         * 1. Set up model parameters
	 * 2. Set the 2Q mapping parameters.
	 * 3. Initialize zero curve and zero banks.
	 */
	void    Initialize(
		TDate 		   todayDt,	// (I) today's date
		KMrParam           &mrPar,	// (I) full dimension mr info
		KSmileParam        &irSmilePar,	// (I) ir smile info
		int		   nIRDim,	// (I) IR dimension
		KVector(int)       &cvTypes, 	// (I) cv types (KV_DIFF...)
		KMap(int, KZCurve) &cv,     	// (I) array of curves
		KMap(int, String)  &cvNames,	// (I) array of curve names
		KMap(int, TDate)   &cvValueDates, // (I) array of cv value dates
		KVector(TDate)     &volDates,   // (I) volatility dates
		KVector(KVector(double)) &factVol, // (I) spot vols
		KResetBank         &resetBank);	// (I) rate reset bank


	/**
	 * Update factor volatilities, including memory allocation for
	 * tree limits, probabilities, tree slices that depends on factor vols.
	 * Used by CET only.
	 */
	void    UpdateFactorVol(
		KVector(TDate) &volDates,           // (I) volatility dates
		KVector(KVector(double)) &factVol); // (I) factor spot vol
 
 
	/**
	 * Free tree memories (such as tree limits, and transitional
	 * probabilities, etc.) that depends on factor vols.
	 * This is used only by CET routine for effective
	 * tree memory management.
	 */
        void    ClearTreeVolMem();


	//------------------------------------------------------
	// Calibration enhancement tool facility
	//------------------------------------------------------

	/**
	 * Master CET calibration routine.
	 * Adjust the input "irVolDiag" volatility diagonal to fit
	 * the corresponding benchmarks options using the CET algorithm.
	 */
	friend	void	KPirTreeCet(
	KPirTree &pirTree,		// (B) tree
  	KMarketCurves& marketCurves,	// (I) curves and curve types
  	KVolDiag& irVolDiag,		// (B) IR volatility data.
  	KMrParam& irMrParam,		// (I) IR mr data.
  	KSmileParam& irSmileParam,	// (I) IR skew data.
  	KResetBank &resetBank);		// (I) rate reset bank


protected:
	//------------------------------------------------------
	// State price access facility
	//------------------------------------------------------

	/**
	 * Insert an express DEV date in the tree. The state prices
	 * are kept and can be accessed with the method GetStatePrice.
	 * WARNING: the state price slices are not automatically
	 * deleted and must be manually removed with the method FreeEDev.
	 */
virtual void		InsertStatePrice(TDate eDevDate)
			{
				Insert(eDevDate);
				mDevDates.push_back(eDevDate);
			}


	/**
	 * Frees the obsolete state price (requested by InsertEDev).
	 */
virtual void		FreeStatePrice(TDate eDevDate)
			{
				KMap(TDate, double*)::iterator it
					= mDevStatePr.find(eDevDate);
				sliceDelete((*it).second);
				mDevStatePr.erase(it);
			}

	/**
	 * Get the state price given a date.
	 */
virtual double*		GetStatePrice(TDate eDevDate); 




	//------------------------------------------------------
	// Protected methods and data
	//------------------------------------------------------
protected:

	/** Calibrates the drift to a ZC curve. 
	 */
virtual	void		CalibrateDrift();

	/** 
	 * TreeSetUp is called after all product related critical dates being
	 * inserted, and does the following:
	 * 1. Run the zero bank date optimation and insert the "optimal"
	 *    dates in the critical date list.
	 * 2. Call KMrNTree::Calibrate to set up the tree
	 *    (Jump size, orthogonal factors, and tree limits, etc.)
	 * 3. Initialize temp discount slices for each curve
	 *    (only after tree limit set up).
	 * 4. Compute the zero prices and forward rates at each time step.
	 * 5. Sort and merge mDevDates
	 */
virtual void		TreeSetUp();



	/**
	 * Free all the memory allocated during tree initialization
         * and calibration.
         * It is used only by CET tool to free the tree memory
         * during each iteration, so to allow re-initialization
         * of tree with new volatilities.
         * It is identical to destructor ~KPirTree(), but keep
         * the object intact.
         */
virtual void            DeleteMemory();



private:

	/**
	 * Set up zero curves and curve names
	 */
	void	InitializeCurves(
		KVector(int)       &cvTypes,// (I) array of cv types(KV_DIFF...)
		KMap(int, KZCurve) &cv,     // (I) array of curves
		KMap(int, String)  &cvNames,// (I) array of curve names
		KMap(int, TDate)   &cvValueDates); //(I) array of cv value dates


	/**
	 * Compute the zero prices and forward rates at each time step.
	 * Should be called after all the time points are set.
	 */
	void		ComputeZeroAndForward();

	/**
	 * Calculates the diffuse discount factor applicable between
	 * tpIdx and tpIdx+1.
	 */
	void		CalcDiscDiff(
			int tpIdx,		// (I) time point index
			int nIRDim,		// (I) IR dimensions
			double Zt,		// (I) center offset 
			double *discount);	// (O) discount slice

	void		SolveOffsetDiff(
			int tpIdx,		// (I) time point index
			int nIRDim,		// (I) IR dimensions
			double *discount,	// (I) slice disc bet t and t+1
			double *statePr,	// (I) slice state price at t
			double zeroPrice,	// (I) zero price at t+1
			double Zt,		// (I) intial offset
			double *DelZt);		// (O) incremental offset 




	/*
	 * Needs to access protected funcs
	 */
friend	void	KPirTreeCalcCet(
		KPirTree &pirTree,		// (I) tree
		KResetBank &resetBank,		// (I) rate reset bank
		KVector(KCetData) &cetData);	// (I) array of cet data



protected:
	//------------------------------------------------------
	// Interest rates specific data
	//------------------------------------------------------
				/**
				 * IR dimension (can be less
				 * than  actual tree dimension
				 */
	int		    mIRDim;
				/** Left q parameter. */
	double		mQLo;
				/** Left q parameter. */
	double		mQHi;
				/** Forward shift parameter. */
	double		mFSh;
				/** Backbone coef. (0=normal, 1-lognormal) */
	double		mBbq;
				/** Total normal model vol. */
	double		mVolNorm;
				/** Total lognormal model vol. */
	double		mVolLogn;




protected:
	//------------------------------------------------------
	// Zero curves and banks
	//------------------------------------------------------

				/** number of zero curves/banks. */
	KVector(int)         mCVTypes;
				/** zero banks [curveIdx]. */
	KMap(int, KZeroBank) mZBanks;
				/** zero curves [curveIdx]. */
	KMap(int, KZCurve)   mKZCurves;
				/** zero curve value dates. */
	KMap(int, TDate)     mValueDates;

				/** discount factor[cvIdx][tpIdx]. */
	KMap(int, double*)	mDiscountIR;

				/** zero coupon price[cvIdx][tpIdx] . */
	KMap(int, double*)   mTpZPrices;
				/** forward rates [cvIdx][tpIdx]. */
	KMap(int, double*)   mTpFRates;
				/** rate mapping offset [0..NbTP]. */
	double	  	    *mTpZCenter;

	//------------------------------------------------------
	// Express Dev tool
	//------------------------------------------------------
				/** Express Dev ON/OFF. */
	bool	  	     mDevOn;
				/** State prices [date][sliceIdx]. */
	KMap(TDate, double*) mDevStatePr;
				/** Express Dev dates. */
	KVector(TDate)       mDevDates;


protected:
				/** Rate reset bank */
	KResetBank	*mResetBank;
};



#endif /* _kpirtree_H */


