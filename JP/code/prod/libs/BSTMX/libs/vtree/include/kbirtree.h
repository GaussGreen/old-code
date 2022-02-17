/***************************************************************
 * Module:	BasisTree
 * Submodule:	
 * File:	kbirtree.h
 * Function:	
 * Author:	
 ***************************************************************/
#ifndef	_kbstree_H
#define	_kbstree_H
#include "kutilios.h"    // Standard definitions & error hadling 
#include "kvtree.h"
#include "ktmxirtree.h"


//--------------------------------------------------------------
/**
 * A basis interest rate tree class based
 * on an pure interest rate tree.
 */

class KBirTree : public KTMXirTree {
public:
	//------------------------------------------------------
	// KVTree base class methods
	//------------------------------------------------------

	/** Default constructor. 
	 */
			KBirTree();

	/** Destructor. 
	 */
virtual			~KBirTree();

        /** Insert a critical date in the tree */
virtual void            Insert(TDate critDate)
			{ KPirTree::Insert(critDate);} 
 
	/** Insert zero dates in the tree and ZBank */
virtual void            Insert(const KZeroReset &zeroReset, bool isCrit)
			{ KPirTree::Insert(zeroReset, isCrit);}

        /** Insert a KRateReset, return the modelReset date as following:
	 *  1. if curve is NON-basis, return the rate reset date.
	 *  2. if curve is KV_BASIS, return the delayed libor reset date.
	 *  3. if curve is KV_SPREAD, return the rate reset date.
	 */
virtual TDate           Insert(const KRateReset &rtReset, bool isCrit);


	/** Insert a KRateReset between reset date and end date,
	 *  and return the rate reset date. 
	 */
virtual TDate   Insert(const KRateReset &rtReset, TDate endDate,  bool isCrit);

	/** Virtual function of BS Tree. 
	 */
virtual	void		Calibrate();

	/** Virtual function of BS Tree.
	 */
virtual	void		Update(int tpIdx);

	/** Create a slice based on curve type:
	 *  1. if curve is BASIS, allocate mIRDim+mBSDim dimension.
	 *  2. if curve is non-BASIS, allocate mIRDim dimension. 
	 */
virtual KTSlice&	TSliceCreate(KTSlice&);



	/** Slice unary operations.  If the two slices have different
	 *  dimensions, the resulting slice ts always conforms to the 
	 *  one with higher dimension.
	 */
virtual KTSlice&	TSliceUnaryOper(
				KTSlice &ts,
				const KTSlice &ts1,
				KOper oper);

	/** Slice DEV based on slice dimension: 
	 *  1. if dimension is mIRDim, same as KPirTree.
	 *  2. if dimension is mIRDim + mBSDim, expand the discount slice
	 *     from mIRDim to mIRDim + mBSDim first, then discount.
	 */
virtual KTSlice&	TSliceDev(KTSlice&, const String&);

	/** Get a zero reset at current time point. 
	 */
virtual void		Get(KTSlice &ts, const KZeroReset &zeroReset);

	/** Get a rate reset at current time point. 
	 *  1. if curve is NON-basis, return IR rate slice.
	 *  2. if curve is KV_BASIS, return basis rate slice.
	 *  3. if curve is KV_SPREAD, return basis spread slice.
	 */
virtual void		Get(KTSlice&, const KRateReset&);


        /**
         * Get IR discount curve name.
         *  a) for IR curve, just return itself,
         * b) for Basis curve, return the IR reference discount curve,
         */
virtual const String&   GetIRDiscCurveName(const String& curveName);



	//------------------------------------------------------
	// Basis tree specific public methods
	//------------------------------------------------------

	/** Get basis spread slice on a specified date. 
	 */
	virtual void	GetSpread(KTSlice &spread, TDate resetDate);


	/** Get basis reset date given a basis delay date. 
	 */
	TDate		GetBasisResetDate(TDate delayDate);

	/** Get basis accural end date given a basis delay date. 
	 */
	TDate		GetBasisAccrualEndDate(TDate delayDate);

	/** Get basis forward spread at specified basis delay date. 
	 */
	double		GetBasisFwdSpread(TDate delayDate, bool isInterp=false);

	/** Get basis forward rate at specified basis delay date. 
	 */
	double		GetBasisFwdRate(TDate delayDate);

	/** Get calibrated center offset of basis spread 
	 *  at specified basis delay date. 
	 */
	double		GetBasisSpreadCenter(TDate delayDate);


	/** Get underlying libor rate with the same maturity as the basis rate
	 *  at specified basis delay date. 
	 */
	KRateReset* 		GetLiborRate(TDate delayDate);

	/** Compute the spread backbone factor at given point */
	double		GetBasisSpreadVolBbq(int t);



	/**
	 * Initialize performs:<br>
	 * 1. Set up model parameters
	 * 2. Set the 2Q mapping parameters.
	 * 3. Initialize zero curve and zero banks.
	 */
	void    Initialize(
		KMarketCurves&  marketCurves,	// (I) curves and curve types

		KVolDiag&	irVolDiag,	// (I) IR volatility data.
		KMrParam&	irMrParam,	// (I) IR mr data.
		KSmileParam&	irSmileParam,	// (I) IR skew data.
		KVolDiag&	bsVolDiag,	// (I) Basis volatility data.
		KMrParam&	bsMrParam,	// (I) Basis mr data.
		KSmileParam&	bsSmileParam,	// (I) Basis skew data.
		double		irBsCorr,	// (I) IR basis correlation.

		KResetBank	&resetBank);	// (I) rate reset bank

    /**
     * Reset BS spot vols
     * This function is only used in TMX tree
     */
	void    ResetBsSpotVol(
		KMarketCurves&  marketCurves,	// (I) curves and curve types

		KVolDiag&	irVolDiag,	// (I) IR volatility data.
		KMrParam&	irMrParam,	// (I) IR mr data.
		KSmileParam&	irSmileParam,	// (I) IR skew data.
		KVolDiag&	bsVolDiag,	// (I) Basis volatility data.
		KMrParam&	bsMrParam,	// (I) Basis mr data.
		KSmileParam&	bsSmileParam,	// (I) Basis skew data.
		double		irBsCorr,	// (I) IR basis correlation.

		KResetBank	&resetBank);	// (I) rate reset bank


	/**
	 * Returns the IR spot volatility vector.
	 */
	KVector(double)	IrSpotVols();

	/**
	 * Returns the basis spot volatility vector.
	 */
	KVector(double)	BsSpotVols();

protected:
	//------------------------------------------------------
	// Basis tree specific methods
	//------------------------------------------------------

	/**
	 * TreeSetUp is called after all product related critical dates being
	 * inserted, and does the following:
	 * 1. Run the zero bank date optimation and insert the "optimal"
	 *    dates in the critical date list.
	 * 2. Call KMrNTree::Calibrate to set up the timeline according
	 *    to the ppy rule and tree parameters (jump size, orthogonal
	 *    factors, and tree limits, etc.)
	 * 3. Initialize temp discount slices for each curve
	 *    (only after tree limit set up).
	 * 4. Compute the zero prices and forward rates at each time step.
	 * 5. Sort and merge mDevDates
	 * 6. Initialize temporary slice for basis rate index.
	 * 7. Sort and merge mBSDelayDates.
	 */
virtual void            TreeSetUp();

	
	/** Check the tree validity */
virtual void		CheckTreeValid();


	/** Calibrate the drift of spread on a reset date
	 *  so that the 1-period basis forward would be priced exactly.
	 *  Store the current basis rate slice as well as drift for
	 *  later use by Get routine. 
	 */
virtual	void		CalibrateSpreadDrift(int tpIdx);// (I) time point index

	/**
	 * Calculate the basis spread at t using the 2-Q mapping function.
	 */
virtual	void		CalcSpread(
			int tpIdx,		// (I) time point index
			double Zt,		// (I) center offset 
			double *spread);	// (O) basis spreads

	/**
	 * Solves for the incremental offset to calibrate
	 * the zero price at tpIdx+1.
	 */
virtual	void		SolveBasisOffset(
			int    tpIdx,           // (I) time point index
			double fwdBasisRate,	// (I) forward basis rate 
			double *spread,         // (I) spread at tpIdx
			double *liborRateIR,    // (I) libor rates at tpIdx
			double *statePr,        // (I) state prices at tpIdx
			double *discountIR,     // (I) discount between t and T
			double zeroPrice,	// (I) discount zero at t+1
			double Zt,		// (I) base shift in expansion
			double *Eta);		// (O) incremental offset 

	/**
	 * Solves for the new spread offset to calibrate
	 * the zero price at tpIdx+1.
     * Use the offset at previous date as initial guess.
	 */
virtual	void		SolveBasisOffsetNew(
			int    tpIdx,           // (I) time point index
			double fwdBasisRate,	// (I) forward basis rate 
			double *liborRateIR,    // (I) libor rates at tpIdx
			double *statePr,        // (I) state prices at tpIdx
			double *discountIR,     // (I) discount between t and T
			double zeroPrice,	    // (I) discount zero at t+1
			double *Zt);		    // (I/O) drift 

	/**
	 * Solves for the incremental offset to calibrate
	 * the zero price at tpIdx+1.
	 * mBSType = SPREAD
	 */
virtual	void		SolveBasisOffsetSpread(
			int    tpIdx,           // (I) time point index
			double fwdBasisRate,	// (I) forward basis rate
			double *spread,         // (I) spread at tpIdx
			double *liborRateIR,    // (I) libor rates at tpIdx
			double *statePr,        // (I) state prices at tpIdx
			double *discountIR,     // (I) discount between t and T
			double zeroPrice,	// (I) discount zero at t+1
			double Zt,		// (I) base shift in expansion
			double *Eta);		// (O) incremental offset 

	/**
	 * Solves for the offset to calibrate
	 * the zero price at tpIdx+1.
	 * mBSType = SPREAD
	 */
virtual	void		SolveBasisOffsetSpreadNew(
			int    tpIdx,           // (I) time point index
			double fwdBasisRate,	// (I) forward basis rate
			double *liborRateIR,    // (I) libor rates at tpIdx
			double *statePr,        // (I) state prices at tpIdx
			double *discountIR,     // (I) discount between t and T
			double zeroPrice,	    // (I) discount zero at t+1
			double *Zt);            // (I/O) drift


	/**
	 * Solves for the incremental offset to calibrate
	 * the zero price at tpIdx+1.
	 * mBSType = PERCENT
	 */
virtual	void		SolveBasisOffsetPercent(
			int    tpIdx,           // (I) time point index
			double fwdBasisRate,	// (I) forward spread
			double *spread,         // (I) spread at tpIdx
			double *liborRateIR,    // (I) libor rates at tpIdx
			double *statePr,        // (I) state prices at tpIdx
			double *discountIR,     // (I) discount between t and T
			double zeroPrice,	// (I) discount zero at t+1
			double Zt,		// (I) base shift in expansion
			double *Eta);		// (O) incremental offset 

	/**
	 * Solves for the offset to calibrate
	 * the zero price at tpIdx+1.
	 * mBSType = PERCENT
	 */
virtual	void		SolveBasisOffsetPercentNew(
			int    tpIdx,           // (I) time point index
			double fwdBasisRate,	// (I) forward basis rate
			double *liborRateIR,    // (I) libor rates at tpIdx
			double *statePr,        // (I) state prices at tpIdx
			double *discountIR,     // (I) discount between t and T
			double zeroPrice,	    // (I) discount zero at t+1
			double *Zt);            // (I/O) drift

protected:

				/** true=basis; false=pure IR. */
	bool	mBSOn;
				/** basis spread dimension. */
	int	mBSDim;

				//--- Basis spread smile parameters 
				/** Basis smile left Q. */
	double	mBSQLo;
				/** Basis right left Q. */
	double	mBSQHi;	
				/** Basis smile forward shift. */
	double	mBSFSh;	
				/** Delay due to averaging. */
	double	mBSDelayShift;

				//--- Basis backbone parameters
				/** Backbone coef (0=normal, 1=lognormal) */
	double	mBSBbq;		
				/** Total normal model vol */
	double	mBSVolNorm;
				/** Total lognormal model vol */
	double	mBSVolLogn;

				/** Basis type: SPREAD or PERCENT. */
	KSpd	mBSType;
				/** Libor index curve name. */
	String	mBasisCVName;
				/** Libor index curve name. */
	String	mLiborCVName;
				/** Basis index rate day count convention. */
	KDayCc	mBasisDCC;
				/** Libor index rate day count convention. */
	KDayCc	mLiborDCC;
				/** Basis discount curve name. */
	String	mBSDiscCVName;

                /** Basis market vol inputs.  Used for estimating the 
                  * initial drift in the back end 
                  */
    KVolDiag mBSVolDiag;

private:

	//------------------------------------------------------
	// Basis reset dates and rates
	//------------------------------------------------------

				/** Delayed libor index reset dates . */
	KVector(TDate)      mBSDelayDates;

				/** Basis reset dates. */
	KMap(TDate, TDate)  mBSResetDates;

				/** Accural end dates. */
	KMap(TDate, TDate)  mBSAEDates;

				/** Basis fwd rates spreads . */
	KMap(TDate, double) mBSTpFSprds;

				/** Basis forward rates. */
	KMap(TDate, double) mBSTpFRates;

				/** Rate mapping offset. */
	KMap(TDate, double) mBSTpCenters;

				/** Libor index rates. */
	KMap(TDate, KRateReset*)	mLiborRates;

				// --- basis bank for spread calc
				/** True/false for spread calculation. */
	bool		    mBSSprdOn;

				/** Store basis rate slice. */
	KCBank		    mBSBank;

				/** Store the tmp basis rate slice. */
	KTSlice		    *mBSTmpRate;

public:

};




#endif /* _kbstree_H */


