/***************************************************************
 * Module:        Risky Tree
 * Submodule:        
 * File:        kcrxtree.h
 * Function:        
 * Author:        
 ***************************************************************/
#ifndef  _kcrtree_H
#define  _kcrtree_H
#include "kutilios.h"    // Standard definitions & error hadling 
#include "kvtree.h"
#include "kpirtree.h"

// forward declarations
class KCrxTreeCet;
class KCrxCetData;


//--------------------------------------------------------------
/**
 * A credit tree class based
 * on an pure interest rate tree.
 */

class KCrxTree : public KPirTree {
public:
        //------------------------------------------------------
        // KVTree base class methods
        //------------------------------------------------------

        /** Default constructor. 
         */
                        KCrxTree();

        /** Destructor. 
         */
virtual                ~KCrxTree();

        /** Insert a critical date in the tree */
virtual void            Insert(TDate critDate)
                        { KPirTree::Insert(critDate);} 
 
        /** Insert zero dates in the tree and ZBank */
virtual void            Insert(const KZeroReset &zeroReset, bool isCrit)
                        { KPirTree::Insert(zeroReset, isCrit);}

        /** Insert a KRateReset, return the modelReset date:
         */
virtual TDate           Insert(const KRateReset &rtReset, bool isCrit);


        /** Insert a KRateReset, return the modelReset date:
         */
virtual TDate           Insert(const KRateReset &rtReset, TDate endDate, bool isCrit);


        /** Virtual function of CIR Tree. 
         */
virtual        void                Calibrate();

        /** Virtual function of CIR Tree.
         */
virtual        void                Update(int tpIdx);

        /** Create a slice based on curve type:
         *  1. if curve is IR, allocate mIRDim dimension. 
         *  2. if curve is BASIS, allocate mIRDim+mBSDim dimension.
         *  2. if curve is CR, allocate mIRDim+mBSDim+mCRDim dimension.
         */
virtual KTSlice&        TSliceCreate(KTSlice&);


        /** Slice unary operations.  If the two slices have different
         *  dimensions, the resulting slice ts always conforms to the
         *  one with higher dimension.
         */
virtual KTSlice&    TSliceUnaryOper(
                KTSlice &ts,
                const KTSlice &ts1,
                KOper oper);



        /** Slice DEV based on slice dimension: 
         *  1. if dimension is mIRDim, same as KPirTree.
         *  2. if dimension is mIRDim + mCRDim, risky discount
         */
virtual KTSlice&        TSliceDev(KTSlice&, const String&);

    /**
     * Slice default dev.  Contains sum of 2 parts:
     * 1. riskless Dev of ts,
     * 2. 1-recovery payment conditional on default in [t, t+1].
     */
virtual KTSlice&     TSliceDefaultDev(
                        KTSlice &ts,
                        double defPayment,
                        const String &discCurveName);



        /** Get a zero reset at current time point. 
         */
virtual void         Get(KTSlice &ts, const KZeroReset &zeroReset);

        /** Get a rate reset at current time point. 
         *  1. if curve is NON-Credit, return IR rate slice.
         *  2. if curve is KV_CREDIT, return credit rate slice.
         */
virtual void         Get(KTSlice&, const KRateReset&);

        /**
         * Get IR discount curve name.
         * a) for IR curve, just return itself,
         * b) for CDS curve, return the IR reference discount curve,
         */
virtual const String&   GetIRDiscCurveName(const String& curveName);

        //------------------------------------------------------
        // Credit tree specific public methods
        //------------------------------------------------------

        /** Get credit forward spread at specified reset date. 
         */
        double       GetCRFwdSpread(TDate delayDate, bool isInterp=false);



        /** Get deterministic credit forward spread at specified reset date. 
         */
        double       GetForwardSpread(int tpIdx);




        /** Get calibrated center offset of credit spread 
         */
        double       GetCRSpreadCenter(TDate delayDate);

        /** Compute the spread backbone factor at given point */
        double       GetCRSpreadVolBbq(int t);



        /**
         * Initialize performs:<br>
         * 1. Set up model parameters
         * 2. Set the 2Q mapping parameters.
         * 3. Initialize zero curve and zero banks.
         */
        void    Initialize(
                KMarketCurves&  marketCurves,        // (I) curves and curve types

                KVolDiag&        irVolDiag,        // (I) IR volatility data.
                KMrParam&        irMrParam,        // (I) IR mr data.
                KSmileParam&        irSmileParam,        // (I) IR skew data.
                KVolDiag&        crVolDiag,        // (I) Credit volatility data.
                KMrParam&        crMrParam,        // (I) Credit mr data.
                KSmileParam&        crSmileParam,        // (I) Credit skew data.
                double                IrCrCorr,        // (I) IR credit correlation.

                KResetBank        &resetBank);        // (I) rate reset bank



        /**
         * Returns the IR spot volatility vector.
         */
        KVector(double)        IrSpotVols();

        /**
         * Returns the credit spot volatility vector.
         */
        KVector(double)        CrSpotVols();

        /*=====================================================================
         * CREDIT CET FRIENDS
         *===================================================================*/
        /**
	    * Master credit CET calibration routine.
	    * Adjust the input "crVolDiag" volatility diagonal to fit
	    * the corresponding benchmarks options using the CET algorithm.
	    */
	    friend	void	KCrxTreeCet(
	        KCrxTree &pirTree,		        // (B) tree
  	        KMarketCurves& marketCurves,	// (I) curves and curve types
  	        KVolDiag& irVolDiag,		    // (I) IR volatility data.
  	        KMrParam& irMrParam,		    // (I) IR mr data.
  	        KSmileParam& irSmileParam,	    // (I) IR skew data.
            KVolDiag& crVolDiag,		    // (B) CR volatility data.
  	        KMrParam& crMrParam,		    // (I) CR mr data.
  	        KSmileParam& crSmileParam,	    // (I) CR skew data.
            double irCrCorr,                // (I) credit/rates correlation
  	        KResetBank &resetBank);		    // (I) rate reset bank

        /*
	    * Needs to access protected funcs
	    */
        friend	void	KCrxTreeCalcCet(
		    KCrxTree &crxTree,		// (I) tree
		    KResetBank &resetBank,		// (I) rate reset bank
		    KVector(KCrxCetData) &cetData);	// (I) array of cet data

protected:
        //------------------------------------------------------
        // Credit tree specific methods
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
         * 6. Initialize temporary slice for credit rate index.
         * 7. Sort and merge mBSDelayDates.
         */
virtual void            TreeSetUp();

        
        /** Check the tree validity */
virtual void                CheckTreeValid();


        /** Calibrate the drift of spread on a reset date
         *  so that the 1-period credit forward would be priced exactly.
         *  Store the current credit rate slice as well as drift for
         *  later use by Get routine. 
         */
virtual        void        CalibrateDrift(); 

    /**
     * Calculates the discount factor applicable between
     * tpIdx and tpIdx+1.
     */
virtual void CalcCRSpreadDiscount(
        int tpIdx,          // (I) time point index
        double Zt,          // (I) center offset
        double *Discount);  // (O) discount slice



        /**
         * Solves for the incremental offset to calibrate
         * the zero price at tpIdx+1.
         */
virtual        void                SolveSpreadOffset(
                        int tpIdx,           // (I) time point index
                        int nDim,            // (I) Spread dimensions
                        double *StatePrRisky,// (I) slice risky state price at t
                        double *DiscountIR,  // (I) IR discount from t to t+1
                        double zeroPrice,    // (I) zero price at t+1
                        double *crZt);       // (I/O) spread shift at t


protected:

                                /** true=credit; false=pure IR. */
        bool          mCROn;
                                /** Credit spread dimension. */
        int           mCRDim;

                                //--- Credit spread smile parameters 
                                /** Credit smile left Q. */
        double        mCRQLo;
                                /** Credit right left Q. */
        double        mCRQHi;        
                                /** Credit smile forward shift. */
        double        mCRFSh;        

                                /** Recovery rate */
        double        mRecovery;

                                //--- Credit backbone parameters
                                /** Backbone coef (0=normal, 1=lognormal) */
        double        mCRBbq;                
                                /** Total normal model vol */
        double        mCRVolNorm;
                                /** Total lognormal model vol */
        double        mCRVolLogn;

                                /** CR discount factor[cvIdx][tpIdx]. */
        double        *mDiscountCR;

                                /** IR discount curve name. */
        String        mIRDiscCVName;

private:

                                /** CR mapping offset. */
        double        *mCRTpZCenter;


                                /** Store protection slice. */
        KZeroBank     mProtBank;


public:

};




#endif /* _kcrxtree_H */


