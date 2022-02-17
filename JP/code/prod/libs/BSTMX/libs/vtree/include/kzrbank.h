/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	kzrbank.h
 * Function:	
 * Author:	David Liu
 ***************************************************************/
#ifndef	_kzrbank_H
#define	_kzrbank_H
#include "kstdinc.h"	/* Standard definitions & error hadling */
#include "kvtree.h"
#include "krate.h"
#include "kcbank.h"

extern  "C" {
#include "interp.h"
};

class	KPirTree;

//--------------------------------------------------------------
/**
 * Zero bank based on claim bank.
 */
class	KZeroBank : public KCBank {
public:

	/** Default constructor */
			KZeroBank();

	/** Copy constructor. */
			KZeroBank(const KZeroBank &zb);

	/** Constructor */
			KZeroBank(KPirTree &vt,
			          const String &zbName = K_DEFAULT_NAME,
			          const String &curveName = K_DEFAULT_NAME,
			          const int    interpType = GTO_LINEAR_INTERP);

	/** Destructor. */
			~KZeroBank();

	/** Set curve index of the zero bank, given a curve name */
	void		SetCurveIdx(const String &curveName)
			{mCurveIdx = mVTree->GetCurveIdx(curveName);}

	/** Get curve index of the zero bank */
	const int	GetCurveIdx() const { return mCurveIdx;}
	
	/** Get curve name of the zero bank */
	const String&	GetCurveName() 
			const { return mVTree->GetCurveName(mCurveIdx);}

	/** Updates (DEV, discard unused) bank at time node "tpIdx". */
virtual void		Update(int   tpIdx);  	// (I) Index of time point 

	/** Get zero values on a time slice */
virtual void		GetZeroSlice(KTSlice &ts,
				     TDate   matDate);

	/** Initialize a set of optimized zero dates given all the
	 *  zero maturity and usage dates, and add the final zero dates
	 *  to the critical date list 
	 */
virtual void 		InitializeZeroDates(); 


	/** Compute a par yield on the tree on rollback process */
virtual void		ComputeParYieldAndAnnuity(
			int tpIdx,			// (I) tree time point
			const KRateReset &rateReset,	// (I) rate index
			KTSlice	&parYield,		// (O) par yield
			KTSlice	&annuity);		// (O) annuity

	/** Return interpolation type */
virtual int     ZeroInterpType() const { return mZeroInterpType;}


	/** Assignment operator */
	KZeroBank& operator=(const KZeroBank &zb);


private:

	/** Process zero dates */
void	ProcessDates(
        int    mNbEvDates,  /* (I) Nb of input mats          */
        TDate  *mEvDates,   /* (I) zero mat datelist         */
        TDate  *mErDates,   /* (I) zero usage list           */
        int    *NbOutMat,   /* (O) Nb of processed zero mats */
        TDate  **Mat,       /* (O) processed zero mats       */
        TDate  **LaMat,     /* (O) latest usage date list    */
        TDate  **ErMat);    /* (O) earliest usage date list */


	/** Optimize zero maturity and usage dates */
void	OptimizeDates(
        int     NbTP,           /* (I) Nb of critical dates        */
        TDate   *TPDates,       /* (I) critical dates              */
        TDate   today);         /* (I) today                       */


private:
					// --- unoptimized zero dates
	int	mCurveIdx;		/** curve index           */
	int	mZeroInterpType;	/** ZC interpolation type */

};

#endif	/* _kzrbank_H */
