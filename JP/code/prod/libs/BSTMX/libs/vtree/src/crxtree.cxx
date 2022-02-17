/***************************************************************
 * Module:        CrxTree
        else if (curveIdx == KV_CREDIT_RISKY)   // CDS 
 * Submodule:        
 * File:        
 * Function:      Credit tree 
 * Author:        David Liu, April 2004
 ***************************************************************/
#include "kstdinc.h"        // Standard definitions & error hadling 
#include "kstlutil.h"       // KCArray 

#define        _kmrntree_SRC
#include "kmrntree.h"
#include "kpirtree.h"
#include "kcrxtree.h"

extern "C" {
#include "cgeneral.h"                   // Has stdlib.h
#include "bastypes.h"                   // TDateList
#include "cerror.h"                     // GtoErrMsg 
#include "cmemory.h"                    // MALLOC, FREE 
#include "extract.h"                    // GtoExtractArray 
#include "macros.h"                     // MAX 

#include "datelist.h"                   // TDateList routines 
#include "ldate.h"                      // GtoDtFwdAny 
#include "dateadj.h"                        // GtoDtBackAdj 
#include "convert.h"                    // GtoFormatDate 
#include "zerodate.h"                   // TZeroDates 
#include "termtype.h"                   // TFloatRateArray 
#include "cfinanci.h"                   // GtoRateToDiscount 
#include "interp.h"                     // GtoRateToDiscount 
#include "zr2simp.h"                    // GtoZerosToSimplePoint 

#include "drlmem.h"                     // DrlDoubleVectAlloc/Free 
#include "drlio.h"                      // 
#include "drltime.h"                    // DrlTDatePrint 
#include "drlvtype.h"                   // DrlVTypeVectAdd 
#include "drlsort.h"                    // DrlTDateArrayFloorIdx() 
#include "drlroot.h"                    // DrlNRPoly

#include "vnfmanly.h"                   // VnfmVolCalib1VArbitrary

};

//
// Q cutoff
//
static  const   double  QCUTOFF = 1e-5;


#ifndef IS_Q
#define IS_Q(q) (fabs((q)) > QCUTOFF)
#endif

#ifndef SHIFT_ZERO
#define SHIFT_ZERO(x) ((x) = (IS_ALMOST_ZERO(x) ? 1e-5 : (x)))
#endif


// >=1. Polynomial order for solving drift of spread.
const        int        NPOLY = 2;        


//--------------------------------------------------------------
//

KCrxTree::KCrxTree() : KPirTree()
{
        mCROn = false;

        mCRDim = 0;

        // Lognormal
        mCRQLo = 1e0;
        mCRQHi = 1e0;
        mCRFSh = 0e0;
        mCRBbq = 0e0;

        mCRTpZCenter = NULL;
        mDiscountCR  = NULL;

        mRecovery    = 0.0;
}




//--------------------------------------------------------------
//

KCrxTree::~KCrxTree()
{

    if (mCROn)
    {
        sliceDelete(mDiscountCR);
        delete [] mCRTpZCenter;
    }

}




//---------------------------------------------------------------
// Initialize credit tree.
 
void
KCrxTree::Initialize(
        KMarketCurves&   marketCurves,     // (I) curves and curve types
        KVolDiag&        irVolDiag,        // (I) IR volatility data.
        KMrParam&        irMrParam,        // (I) IR mr data.
        KSmileParam&     irSmileParam,     // (I) IR skew data.
        KVolDiag&        crVolDiag,        // (I) Credit volatility data.
        KMrParam&        crMrParam,        // (I) Credit mr data.
        KSmileParam&     crSmileParam,     // (I) Credit skew data.
        double           irCrCorr,         // (I) IR crdit correlation.
        KResetBank       &resetBank)// (I) rate reset bank
{
static  char            routine[] = "KCrxTree::Initialize";
        int             nDim,
                        irDim,
                        crDim;              // Full dimension of credit tree
        int             idx;
        KMrParam        treeMrParam;        // full tree mr parameters

        
 try{

        KVector(TDate)           factVolDates;       // spot vol dates
        KVector(KVector(double)) factVolCurves;      // factor spot vol curves


        // Check inputs
        if (!marketCurves.IsValid() ||
            !irVolDiag.IsValid()    ||
            !irMrParam.IsValid()    ||
            !irSmileParam.IsValid() ||
            !crVolDiag.IsValid()    ||
            !crMrParam.IsValid()    ||
            !crSmileParam.IsValid()) 
/*
            !IsConsistent(irVolDiag, crVolDiag))
*/
                throw KFailure("Invalid market or model inputs.\n");



        //----------------------------------------------
        // (1) Build the full mr model parameters
        //----------------------------------------------
        treeMrParam = Correlate(irMrParam, crMrParam, irCrCorr);
        nDim  = treeMrParam.mNumFact;
        crDim = crMrParam.mNumFact;
        irDim = nDim - crDim;

        

        //----------------------------------------------
        // (2) Turn credit flag ON if credit curve
        //     or spread is involved in the product
        //     Insert additional basis curve if only
        //     basis spread curve is specified in the env.
        //----------------------------------------------
 
        mCROn = marketCurves.IsCredit();


        //----------------------------------------------
        // (3) Perform IR volatility calibration
        //----------------------------------------------

        //
        // Get the diffused curve for the volatility calibration
        //
        KZCurve &irDiffCurve = marketCurves.GetIRDiffuse();

        // Calibrate IR spot vols in irDim dimensions
        DppCalibIRVolDiag(
                irMrParam,
                irSmileParam,
                irDim,
                irVolDiag,
                irDiffCurve,
                factVolDates,
                factVolCurves);

        
        //----------------------------------------------
        // (4) Assign credit model parameters
        //     if credit curve or spread is involved in the product
        //----------------------------------------------
        if (mCROn) {
                                
                mCRDim = crDim;
                mCRQLo = crSmileParam.mQ1;
                mCRQHi = crSmileParam.mQ2;
                mCRFSh = crSmileParam.mQF;

                // Credit Backbone parameters
                mCRBbq = crMrParam.mBackboneQ;
                
                // Compute reference vol constants
                double norm = 0e0;
                for (int idx = 0; idx<crDim; idx++)
                        norm += crMrParam.mAlpha[idx] * crMrParam.mAlpha[idx];

                norm = sqrt(norm);
                if (fabs(norm) < DBL_EPSILON)
                        throw KFailure("%s: total alpha too small %g.\n",
                                routine, norm);
 
                if (IS_ALMOST_ZERO(mCRBbq - 1e0)) {
                        mCRVolNorm = 0.;
                        mCRVolLogn = norm;
                } else if (IS_ALMOST_ZERO(mCRBbq - 0e0)) {
                        mCRVolNorm = norm;
                        mCRVolLogn = 0.;
                } else {
                        throw KFailure("%s: CR backbone must be 0 or 1 (%g).\n",
                                routine, mCRBbq);
                }

                // Check if mIRDiscCVName is in the table.
                // to be done ...
                //
                mIRDiscCVName = marketCurves.mIRDiscCVName;
                
                // Recovery
                mRecovery = marketCurves.mRecovery;

                //
                // Calibrate credit spot vols in crDim dimensions
                // And add credit spot vol curve to total factor vols
                //
                DppCalibCRSpreadVolDiag(
                     crMrParam,
                     irMrParam,
                     irCrCorr,
                     crSmileParam,
                     marketCurves.GetDiffuse().BaseDate(),
                     marketCurves.mCV[marketCurves.GetCurveType(mIRDiscCVName)],
                     marketCurves.mCV[KV_CREDIT_RISKY],
                     crVolDiag,
                     mRecovery,
                     factVolDates,
                     factVolCurves);
            
                
        }



        // Print the vol factors
        if (debugLevel > DEBUG_LEVEL_VOLATILITY) {
            dppLog << endl;
            dppLog << "============ FACTOR VOLS ============" << endl;
            dppLog << "Date        ";
            for (idx=1; idx<=nDim; ++idx)
                   dppLog << format("Factor %-3d", idx);

            dppLog << endl;

            for (idx=0; idx<=factVolDates.size()-1; idx++)
            {
                    dppLog << format("%-12s",
                        GtoFormatDate(factVolDates[idx]));

                    for (KVector(KVector(double))::iterator 
                        itVol=factVolCurves.begin();
                     itVol != factVolCurves.end();
                     ++itVol)
                        dppLog << format("%-10.6f", (*itVol)[idx]);

                    dppLog << endl;
            }
        }
        

        //----------------------------------------------
        // (5) Initialize IR tree
        //----------------------------------------------

        KPirTree::Initialize(
                        marketCurves.mToday,
                        treeMrParam,
                        irSmileParam,
                        irDim,
                        marketCurves.mCVTypes,
                        marketCurves.mCV,       
                        marketCurves.mCVNames, 
                        marketCurves.mValueDates, 
                        factVolDates,   
                        factVolCurves,
                        resetBank);

        
        // Update the tree engine for risky zero bank
        // Risky zero bank update and calculation of risky
        // annuity will be done automatically
        if (mCROn) {
            int zcInterpType = GTO_FLAT_FORWARDS;

            KMap(int, KZeroBank)::iterator where;
            where = mZBanks.find(KV_CREDIT_RISKY);
            if (where != mZBanks.end())
	    {
		zcInterpType = where->second.ZeroInterpType();
                mZBanks.erase(KV_CREDIT_RISKY);
	    }

            mZBanks.insert(KMap(int, KZeroBank)::value_type(
                                               KV_CREDIT_RISKY,
                                               KZeroBank(*this,
                                               format("Risky Zero Bank %d",
                                               KV_CREDIT_RISKY),
                                               GetCurveName(KV_CREDIT_RISKY),
					       zcInterpType)));

            //
            // Add 2 internal credit banks for default probability
            // and contingent payment calculation.

            //
            // 1. CR default bank for default probability calculation
            //
            
            // Add curve type KV_CREDIT_DEFPROB for convenience
            //
            MapCurveName("DEF_PROB_"+GetCurveName(KV_CREDIT_RISKY), 
                         KV_CREDIT_DEFPROB);


            mZBanks.insert(KMap(int, KZeroBank)::value_type(
                                            KV_CREDIT_DEFPROB,
                                            KZeroBank(*this,
                                            format("Default Prob Zero Bank %d",
                                            KV_CREDIT_DEFPROB),
                                            GetCurveName(KV_CREDIT_DEFPROB),
					    zcInterpType)));

            //
            // 2. CR protection bank for par spread calculation
            //
            
            // Add curve type KV_CREDIT_PROT_DEFRECOV to use default
            // recovery rate given in the mkt
            //
            MapCurveName("PROT_DEFRECOV_"+GetCurveName(KV_CREDIT_RISKY), 
                         KV_CREDIT_PROT_DEFRECOV);


            mZBanks.insert(KMap(int, KZeroBank)::value_type(
                                            KV_CREDIT_PROT_DEFRECOV,
                                            KZeroBank(*this,
                                            format("Protection Zero Bank %d",
                                            KV_CREDIT_PROT_DEFRECOV),
                                            GetCurveName(KV_CREDIT_PROT_DEFRECOV),
					    zcInterpType)));
            //
            // 3. CR protection bank for protection leg calculation
            //
            
            // Add curve type KV_CREDIT_PROT_BINRECOV to use binary 
            // recovery rate override given in the deal.
            //
            MapCurveName("PROT_BINRECOV_"+GetCurveName(KV_CREDIT_RISKY), 
                         KV_CREDIT_PROT_BINRECOV);


            mZBanks.insert(KMap(int, KZeroBank)::value_type(
                                            KV_CREDIT_PROT_BINRECOV,
                                            KZeroBank(*this,
                                            format("Protection Zero Bank %d",
                                            KV_CREDIT_PROT_BINRECOV),
                                            GetCurveName(KV_CREDIT_PROT_BINRECOV),
					    zcInterpType)));
        }

        
    }
    catch (KFailure) {
		
        throw KFailure("%s: failed.\n", routine);
    }
	
   
}




const String&
KCrxTree::GetIRDiscCurveName(const String& curveName)
{
        int discCurveIdx;

        discCurveIdx = GetCurveIdx(curveName);

        if (discCurveIdx < KV_CREDIT_RISKY)
                return curveName;
        else
                return mIRDiscCVName;

}



//---------------------------------------------------------------
//
//

KVector(double)
KCrxTree::IrSpotVols()
{
        return mFactVol[0];
}




//---------------------------------------------------------------
//


KVector(double)
KCrxTree::CrSpotVols()
{
        if (mCRDim <= 0) {
                throw KFailure("%s: no credit spread dimension %d (=%d-%d).\n",
                        "KCrxTree::CRSpotVols",
                        mCRDim, mNbFactor, mIRDim);
        }
        return mFactVol[mIRDim];
}



//---------------------------------------------------------------
// Check the tree.
//

void
KCrxTree::CheckTreeValid()
{
static  char    routine[] = "KCrxTree::CheckTreeValid";
 
        KPirTree::CheckTreeValid();

        if (mCROn) {
           // Dimension consistency
           if ((mIRDim+mCRDim) > 3 || (mIRDim+mCRDim) == 0)
                throw KFailure("%s: total dimension (%d+%d) of the tree is "
                               " either greater than 3 or equal to 0.\n",
                                routine, mIRDim, mCRDim);
        
           if (GetCurveIdx(mIRDiscCVName) == K_DEFAULT_IDX)
                throw KFailure("%s: IR discount curve %s does not exist.\n",
                                routine, mIRDiscCVName.c_str());

           if (debugLevel > DEBUG_LEVEL_TIMELINE) {
                dppLog << format("%s:\n", routine);
                dppLog << format("mCRDim = %d\n",         mCRDim);
                dppLog << format("mRecovery = %lf\n",     mRecovery);
                dppLog << format("mCRQLo = %lf\n",        mCRQLo);
                dppLog << format("mCRQHi = %lf\n",        mCRQHi);
                dppLog << format("mCRFSh = %lf\n",        mCRFSh);
                dppLog << format("mCRBbq = %lf\n",        mCRBbq);
                dppLog << format("mCRVolNorm = %lf\n",    mCRVolNorm);
                dppLog << format("mCRVolLogn = %lf\n",    mCRVolLogn);
           }


        }

}




//--------------------------------------------------------------
// Allocate a slice based on curve type
// DO NOT call directly for slice allocation. 
// Use KTSlice constructor instead
 
KTSlice&
KCrxTree::TSliceCreate(KTSlice &ts)
{
static  char    routine[] = "KCrxTree::TSliceCreate";
 
        int     sliceDim;
        int        curveIdx;
        int        nDim;
 
 try{
        // Check that we have anything to do
        if (!ts.IsEmpty()) return(ts);
 
        nDim = mIRDim + mCRDim;

        curveIdx = ts.GetCurveIdx();

        // Allocate slice dimension based on the curve type:
        // 1. if curve is CREDIT, allocate mIRDim+mCRDim dimension.
        // 2. if curve is non-CREDIT, allocate mIRDim dimension.
        if (curveIdx >= KV_CREDIT_RISKY)        
        {
                sliceDim = nDim;
        }
        else
                sliceDim = mIRDim;
 
        return TSliceDimCreate(ts, sliceDim);

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}




//--------------------------------------------------------------
// Performs a unary operation on a timeslice.
//

KTSlice&
KCrxTree::TSliceUnaryOper(
        KTSlice &ts,
        const KTSlice &ts1,
        KOper oper)
{
static        char        routine[] = "KCrxTree::TSliceUnaryOper";
        int        tpIdx = TPIdxCurrent();

        int        curveIdx;

 try{
        // Ensure slice allocated
        if (ts1.IsEmpty()) 
                throw KFailure("%s: while performing operation `%s' "
                        " on slice `%s', argument slice `%s' is empty.\n",
                        routine, printOper(oper), 
                        ts.GetSliceName().c_str(), 
                        ts1.GetSliceName().c_str());


        // Ensure target ts is allocated
        TSliceCreate(ts);


        // Compare the dimensions of two slices, and
        // always conform to the one with higher dimension.
        //
        if(ts.GetSliceDim() == ts1.GetSliceDim()) {
                KMrNTree::TSliceUnaryOper(
                                ts,
                                ts1,
                                oper);
        }
        else if (ts.GetSliceDim() > ts1.GetSliceDim()) {
                // temporary slice with dimension of ts.GetSliceDim()
                // Use curve idx (not curve name because name is not
                // always updated during slice operation) to 
                // get the right dimension info.
                //
                curveIdx = ts.GetCurveIdx();

                KTSlice tmpTS1(*this,
                               "Temporary slice for "
                               "dimension resize",
                               GetCurveName(curveIdx));

                // Expand dim of ts1 to dim of ts.
                //
                TSliceExpand(tpIdx,                
                             ts1,
                             tmpTS1);

                // Perform the operation
                //
                KMrNTree::TSliceUnaryOper(
                                ts,
                                tmpTS1,
                                oper);
        }
        else {                // ts.GetSliceDim() < ts1.GetSliceDim()
                // temporary slice with dimension of ts.GetSliceDim()
                //

                curveIdx = ts1.GetCurveIdx();

                KTSlice tmpTS(*this,
                              "Temporary slice for "
                              "dimension resize",
                              GetCurveName(curveIdx));

                TSliceExpand(tpIdx,                
                             ts,
                             tmpTS);


                // Reallocate the slice memory
                TSliceDestroy(ts);

                // Set the new slice dimension of ts.
                //
                // !!!!! Riskless value would inherit risky and contingent 
                // discount effects.  Make sense?!
                //
                ts.SetCurveIdx(curveIdx);
                TSliceCreate(ts);


                // Copy the expanded slice data.
                //
                sliceUnaryOper(
                                (double*) ts.mData,
                                (long) ts.GetSliceDim(),
                                (double*) tmpTS.mData,
                                COPY);

                
                // Perform the operation
                //
                KMrNTree::TSliceUnaryOper(
                                ts,
                                ts1,
                                oper);

        }

        return(ts);

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}



//--------------------------------------------------------------
// Slice dev.          Should be called after tree being updated, i.e.
// the discount factor and probabality are calculated. 
// The discount curve index can only be diffuse curve or 
// curves with deterministic spreads with respect to diffuse curve.
//
 
KTSlice&
KCrxTree::TSliceDev(KTSlice &ts, const String &discCurveName)
{
static  char    routine[] = "KCrxTree::TSliceDev";
 
        int     discCurveIdx;                // discount curve index
        int     irDiscountIdx;               // Risky IR discount curve index
        int     sliceDim;                    // slice dimension
        int     nDim;                        // full dimension of the tree
        int     t;
        int     tpIdxTree,                   // current time point of the slice
                tpIdxSlice;                  // current time point of the tree
        
        double  *tmpSlice = NULL;            // temp slice in full dimension
        double  *discountIR = NULL,
                *discountRisky = NULL;

 try{

        nDim = mIRDim + mCRDim;

        tpIdxTree  = this->TPIdxCurrent();
        tpIdxSlice = ts.GetTpIdx();


        // Dev empty slices does nothing
        if (ts.mData == NULL) return ts;
 
        // Don't do anything at last timept
        if (this->TPIdxCurrent() == TPNum()) return ts;
 
        // Check timepoint consistency
        if (tpIdxSlice != tpIdxTree+1)
                throw KFailure("%s: attempting to Ev slice `%s' "
                        "with TP %d on tree at TP %d.\n",
                         routine, ts.GetSliceName().c_str(),
                        tpIdxSlice, tpIdxTree);
 

        t = tpIdxTree;

        // Slice dimension
        //

        sliceDim = ts.GetSliceDim();


        // Discount curve index

        discCurveIdx = GetCurveIdx(discCurveName);
        
        // riskless IR discount
        if (discCurveIdx < KV_CREDIT_RISKY)
        {    

            //
            // Discount slice dimension is mIRDim, need to dev slice 
            // differently based on its dimension:
            // 1. pure IR slice, straightforward.
            // 2. risky slice, need to expand the discount slice in nDim.
            //

            if (sliceDim == mIRDim) {
                KPirTree::TSliceDev(ts, discCurveName);
            }
            else if(sliceDim == nDim) {

                tmpSlice = sliceNew(nDim);

                // Expand discount slice from mIRDim to nDim. 
                sliceExpand(t,
                            nDim,
                            mIRDim,
                            KPirTree::GetDiscount(discCurveIdx),
                            tmpSlice);

                KMrNTree::sliceEv(
                                  (double*) ts.mData,
                                  nDim,
                                  tmpSlice,
                                  t);

                sliceDelete(tmpSlice);
            }
            else {
                throw KFailure("%s: invalid slice dimension (%d) in the tree"
                               " where nDim=%d, mCRDim=%d.\n", 
                                routine, sliceDim, nDim, mCRDim);
            }
        }
        else if (discCurveIdx == KV_CREDIT_RISKY)   // Risky discount
        {

            // ts with lower dimension
            //
            if (sliceDim == mIRDim) {
                tmpSlice = sliceNew(nDim);

                // Expand slice from mIRDim to nDim. 
                sliceExpand(t+1,
                            nDim,
                            mIRDim,
                            (double *)(ts.mData),
                            tmpSlice);

                sliceDelete((double *)(ts.mData));

                //Reassign the slice dimension and value 
                ts.SetSliceDim(nDim);
                ts.SetCurveIdx(discCurveName);
                ts.mData = tmpSlice;
            }

            // IR discount curve
            irDiscountIdx = GetCurveIdx(mIRDiscCVName);

            discountIR = KPirTree::GetDiscount(irDiscountIdx),

            //
            // Compute risky discount factor
            //
            discountRisky = sliceNew(nDim);

            // Expand IR discount
            sliceExpand(t,
                        nDim,
                        mIRDim,
                        discountIR,
                        discountRisky);

            // Multiply by survival prob            
            sliceUnaryOper(discountRisky,
                           nDim,
                           mDiscountCR,
                           MULT);

            // Dev
            KMrNTree::sliceEv( (double*) ts.mData,
                               nDim,
                               discountRisky,
                               t);
            
            sliceDelete(discountRisky);

        }
        else if (discCurveIdx == KV_CREDIT_DEFPROB)   // Default (survival) prob
        {

            // ts with lower dimension
            //
            if (sliceDim == mIRDim) {
                tmpSlice = sliceNew(nDim);

                // Expand slice from mIRDim to nDim. 
                sliceExpand(t+1,
                            nDim,
                            mIRDim,
                            (double *)(ts.mData),
                            tmpSlice);

                sliceDelete((double *)(ts.mData));

                //Reassign the slice dimension and value 
                ts.SetSliceDim(nDim);
                ts.SetCurveIdx(discCurveName);
                ts.mData = tmpSlice;
            }

            // Dev - multiply by conditional survival probability
            KMrNTree::sliceEv( (double*) ts.mData,
                               nDim,
                               mDiscountCR,
                               t);
            

        }
        else if (discCurveIdx == KV_CREDIT_PROT_DEFRECOV)
        {
            // contingent discoun with default recovery rate
            TSliceDefaultDev(ts,
                             1e0 - mRecovery,
                             discCurveName);
        }
        else if (discCurveIdx == KV_CREDIT_PROT_BINRECOV)
        {
            // contingent discoun binary default recovery rate
            // given in the deal
            TSliceDefaultDev(ts,
                             1e0,
                             discCurveName);
        }



        // Set new TP for the slice
        //

        ts.SetTpIdx(t);


        return(ts);

    }
    catch (KFailure) {
        sliceDelete(discountRisky);
        throw KFailure("%s: failed.\n", routine);
    }
}




//--------------------------------------------------------------
// Slice default dev.  Contains sum of 2 parts:
// 1. riskless Dev of ts,
// 2. payment conditional on default in [t, t+1].
//
 
KTSlice&
KCrxTree::TSliceDefaultDev(
    KTSlice &ts, 
    double defPayment,
    const String &discCurveName)
{
static  char    routine[] = "KCrxTree::TSliceDefaultDev";
 
        int     discCurveIdx;            // discount curve index
        int     irDiscountIdx;           // Risky IR discount curve index
        int     sliceDim;                // slice dimension
        int     nDim;                    // full dimension of the tree
        int     t;
        int     tpIdxTree,               // current time point of the slice
                tpIdxSlice;              // current time point of the tree
        
        double  *tmpSlice      = NULL;        // temp slice in full dimension
        double  *discountIR    = NULL,
                *logCR         = NULL,
                *logRisky      = NULL,
                *discountCRTmp = NULL,
                *discountRisky = NULL;

 try{

        nDim = mIRDim + mCRDim;

        tpIdxTree  = this->TPIdxCurrent();
        tpIdxSlice = ts.GetTpIdx();


        // Dev empty slices does nothing
        if (ts.mData == NULL) return ts;
 
        // Don't do anything at last timept
        if (this->TPIdxCurrent() == TPNum()) return ts;
 
        // Check timepoint consistency
        if (tpIdxSlice != tpIdxTree+1)
                throw KFailure("%s: attempting to Ev slice `%s' "
                        "with TP %d on tree at TP %d.\n",
                         routine, ts.GetSliceName().c_str(),
                        tpIdxSlice, tpIdxTree);
 
        t = tpIdxTree;


        // Discount curve index
        discCurveIdx = GetCurveIdx(discCurveName);

        // Must be risky discount curve
        
        if (discCurveIdx != KV_CREDIT_PROT_DEFRECOV &&
            discCurveIdx != KV_CREDIT_PROT_BINRECOV ) 

        {    
            throw KFailure("%s: must be risky contingent curve (%s)\n",
                           routine,
                           discCurveName.c_str());
        }

        // Slice dimension
        //
        sliceDim = ts.GetSliceDim();

        // If slice is lower dimension, then expand to nDim
        //
        if (sliceDim == mIRDim)
        {
            if ((tmpSlice = sliceNew(nDim)) == NULL)
            {
                throw KFailure("%s: failed to allocte tmpSlice\n", routine);
            }

            // Expand slice from mIRDim to nDim. 
            sliceExpand(t+1,
                        nDim,
                        mIRDim,
                        (double *)(ts.mData),
                        tmpSlice);

            sliceDelete((double *)(ts.mData));

            // Reassign the slice dimension and value 
            ts.SetSliceDim(nDim);
            ts.SetCurveIdx(discCurveName);
            ts.mData = tmpSlice;
        }

        // IR discount curve
        irDiscountIdx = GetCurveIdx(mIRDiscCVName);

        // 1-period IR discount
        discountIR = KPirTree::GetDiscount(irDiscountIdx),

        // 1-period risky discount factor
        //
        discountRisky = sliceNew(nDim);
        discountCRTmp = sliceNew(nDim);
        logCR         = sliceNew(nDim),
        logRisky      = sliceNew(nDim),

        //
        // 1-period default probability in flat CR and IR assumption: 
        // (1-Z_CR*Z_IR)*ln(Z_CR)/(ln(Z_CR*Z_IR)
        // This is a bit overdone, but consistent with the 
        // the deterministic calculation, even when time step is large
        //
        // mDiscountCR computed by tree.Update() already
        //

        // IR discounting
        // Expand IR discount
        sliceExpand(t,
                    nDim,
                    mIRDim,
                    discountIR,
                    discountRisky);

        sliceUnaryOper(discountRisky,
                       nDim,
                       mDiscountCR,
                       MULT);


        // ln(Z_CR*Z_IR)
        sliceUnaryOper(logRisky,
                       nDim,
                       discountRisky,
                       COPY);

        sliceScalarOper(logRisky,
                        nDim,
                        1.,    // dummy
                        LOG);
        // ln(Z_CR)
        sliceUnaryOper(logCR,
                       nDim,
                       mDiscountCR,
                       COPY);

        sliceScalarOper(logCR,
                        nDim,
                        1.,    // dummy
                        LOG);
        

        // Use discountRisky to compute
        // (1-Z_CR*Z_IR)*ln(Z_CR)/(ln(Z_CR*Z_IR)
        //
        sliceScalarOper(discountRisky,
                        nDim,
                        1.,
                        SUB);

        sliceUnaryOper(discountRisky,
                       nDim,
                       logCR,
                       MULT);
        
        sliceUnaryOper(discountRisky,
                       nDim,
                       logRisky,
                       DIV);

        // Flip the sign
        sliceScalarOper(discountRisky,
                        nDim,
                        -1.,
                        MULT);


/*
        // SIMPLE calculation of protection payment in
        // given interval: (1-Z_CR)*(1-R)
        //
        sliceScalarOper(discountRisky,
                       nDim,
                       1.,
                       COPY);
        sliceUnaryOper(discountRisky,
                       nDim,
                       mDiscountCR,
                       SUB);
*/


        // Weighted by default payoff
        sliceScalarOper(discountRisky,
                        nDim,
                        defPayment,
                        MULT);



/*
        sliceScalarOper(discountRisky,
                        nDim,
                        1.,
                        COPY);

        sliceUnaryOper(discountRisky,
                       nDim,
                       mDiscountCR,
                       SUB);
                        

        // IR discounting
        // Expand IR discount
        sliceExpand(tpIdxSlice,
                    nDim,
                    mIRDim,
                    discountIR,
                    discountIR_nDim);

        // (1-Z_CR)*Z_IR
        sliceUnaryOper(discountRisky,
                       nDim,
                       discountIR_nDim,
                       MULT);

        // Weighted by default payoff
        sliceScalarOper(discountCRTmp,
                        nDim,
                        defPayment,
                        COPY);

        // 1-period payment conditional on default
        KMrNTree::sliceEv( discountCRTmp,
                           nDim,
                           discountRisky,
                           t);
*/

        // Riskless Dev
        String    crDiscCVName = GetCurveName(KV_CREDIT_RISKY);

        TSliceDev(ts, crDiscCVName);    


        // Sum of the two
        sliceUnaryOper((double *)(ts.mData),
                       nDim,
                       discountRisky,  //discountCRTmp,
                       ADD);
        
        // Free memeory
        sliceDelete(discountRisky);
        sliceDelete(discountCRTmp);
        sliceDelete(logRisky);
        sliceDelete(logCR);


        // Set new TP for the slice
        //

        ts.SetTpIdx(tpIdxTree);


        return(ts);

    }
    catch (KFailure) {
        sliceDelete(discountRisky);
        sliceDelete(discountCRTmp);
        sliceDelete(logRisky);
        sliceDelete(logCR);
        throw KFailure("%s: failed.\n", routine);
    }
}




//--------------------------------------------------------------
//
void
KCrxTree::Update(int tpIdx)
{
static  char    routine[] = "KCrxTree::Update";
        
        int t = tpIdx;          // Time step index

        double  irZt,
                crZt;

 try{

        // Don't do anything at last timept
        if (tpIdx != TPNum())
        {
            
            // Compute the transition probability between t and t+1
            KMrNTree::Update(t);

            // Calibrated center offset
            irZt = mTpZCenter[t];

            // Calculate the discount between t and t+1
            //
            KPirTree::CalcDiscount(t, mIRDim, irZt, mDiscountIR);

            if (mCROn) {
                // Calibrated spread center offset
                crZt = mCRTpZCenter[t];

                // Calculate the CR discount between t and t+1
                //
                CalcCRSpreadDiscount(t, crZt, mDiscountCR);
            }
        }


        // Update the zero bank
        for(KMap(int, KZeroBank)::iterator itZB  = mZBanks.begin();
                           itZB != mZBanks.end(); ++itZB)
            (*itZB).second.Update(t);

    }
                
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}



//--------------------------------------------------------------
// TreeSetUp is called after all product related critical dates being
// inserted, and does the following:
// 1. Run the zero bank date optimation and insert the "optimal"
//    dates in the critical date list.
// 2. Call KMrNTree::Calibrate to set up the timeline according to the ppy
//    rule and tree parameters (jump size, orthogonal factors,
//    and tree limits, etc.)
// 3. Initialize temp discount slices for each curve
// 4. Compute the zero prices and forward rates at each time step.
// 5. Sort and merge mDevDates
//
void
KCrxTree::TreeSetUp()
{
static  char    routine[] = "KCrxTree::TreeSetUp";

        
 try{
        if (debugLevel > DEBUG_LEVEL_DRIFT) 
        { 
                dppLog << endl;
                dppLog << format("============%s===========", routine) << endl;
                dppLog << endl;
        }

        KPirTree::TreeSetUp();

        // Allocate CR spread memory
        if (mCROn) {
            ASSERT_OR_THROW((mDiscountCR = sliceNew(mIRDim + mCRDim)) != NULL);

            ASSERT_OR_THROW((mCRTpZCenter = new double [NbTP+2]) != NULL);

        } // mCROn

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}



//--------------------------------------------------------------
//
void
KCrxTree::Calibrate()
{
static  char    routine[] = "KCrxTree::Calibrate";
        
 try{

        TreeSetUp();

        CheckTreeValid();

        /* Calibrate IR and spread drifts */
        CalibrateDrift();

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}



//--------------------------------------------------------------
// Calibrate both the IR and spread drift 
// and save the state price in nDim tree slice.

void
KCrxTree::CalibrateDrift()
{
static  char    routine[] = "KCrxTree::CalibrateDrift";

    int t;          // Time step index

    int     nDim;

    TDate   currDate;

    double  irFwdShift, 
            crFwdShift;

    double  *StatePrRisky    = NULL,
            *StatePrIR       = NULL,
            *StatePrRiskyL   = NULL,
            *StatePrIRL      = NULL,
            *DiscountIR      = NULL,
            *DiscountRisky   = NULL;


    double  *devStatePr      = NULL,
            *devStatePrL     = NULL;

    double  irZeroPrice,
            crZeroPrice,
            riskyZeroPrice,
            p;

    double  irZt,
            irDelZt,
            crZt;

    double  FwdSprd;

    KVector(TDate)::iterator iterDev;

    try {

    nDim   = mNbFactor; // Total dimension

    if (!mCROn)
    {
        KPirTree::CalibrateDrift();
        return;
    }

    // Check dimension valid
    if (nDim   < 1 || nDim   > 3 ||
        mIRDim < 1 || mIRDim > 3 ||
        mCRDim > 3)
        throw KFailure("%s: invalid factor numbers.\n"
                "Total dimension is %d, IR dimension is %d.\n",
                routine, nDim, mIRDim);

    irFwdShift = mFSh;
    crFwdShift = mCRFSh;

    //
    // Allocate tmp working space
    //
    StatePrIR       = sliceNew(mIRDim);
    StatePrRisky    = sliceNew(nDim);
    DiscountRisky   = sliceNew(nDim);

    // Pointer to IR diffuse discount slice
    DiscountIR      = KPirTree::GetDiscount(KV_DIFF);

    //
    // Set initial state price to 1 at TPIdx = 0.
    //
    StatePrIRL = StatePrIR + NodeOffset(mIRDim, 0, 0, 0);
    StatePrIRL[0] = 1e0;

    StatePrRiskyL = StatePrRisky + NodeOffset(nDim, 0, 0, 0);
    StatePrRiskyL[0] = 1e0;


    //
    // Set up ExpressDEV (r) tool.
    //
    mDevOn = false;

    if( !mDevDates.empty()){
        mDevOn = true;
        iterDev = mDevDates.begin();


        // If the start date is a express dev date
        if ( *iterDev == TPToday())
        {
            devStatePr = sliceNew(nDim);

            devStatePrL = devStatePr + NodeOffset(nDim, 0, 0, 0);
            devStatePrL[0] = 1e0;

            mDevStatePr.insert(KMap(TDate, double*)::value_type(
                       TPToday(), devStatePr));

            devStatePr = NULL;

            iterDev++;
        }
    }

    // initial guess of forward shift
    // It will be overwritten in CalcDiscDiff to
    // ensure that the rate distribution at X=0
    // corresponds to the forward rate
    irZt = irFwdShift;
    crZt = crFwdShift;

    //
    // Calibrate the drifts of both IR and CR at each point in the forward loop
    //
    for (t = 0; t < NbTP; t++)
    {
        tpIdxCurrent = t;   // set the current time point

        currDate = TPDate(t);

        //
        // Calculate the discount between t and t+1
        // using the previous offset as a first guess
        //
        KPirTree::CalcDiscount(t, mIRDim, irZt, mDiscountIR);

        if (debugLevel >= DEBUG_LEVEL_GEOMETRY) {
            dppLog << format("%s: discount at tpIdx %3d "
                "before solve.\n", routine, t);
            slicePrint(
                    mDiscountIR[KV_DIFF],
                    mIRDim,
                    t,
                    FALSE,
                    dppLog);
        }


        //
        // Solve for the offset to calibrate zero at t+1
        // using previous offset as initial guess
        //
        irZeroPrice = GetZeroPrice(KV_DIFF)[t+1];
        crZeroPrice = GetZeroPrice(KV_CREDIT_RISKY)[t+1];

        riskyZeroPrice = irZeroPrice * crZeroPrice;


        KPirTree::SolveOffset(
                t,
                mIRDim,
                DiscountIR,
                StatePrIR,
                irZeroPrice,
                irZt,
                &irDelZt);

        //
        // Store new drift offset
        //
        irZt += irDelZt;
        mTpZCenter[t] = irZt;


        //
        // Recompute the exact IR discount factor slice between (t and t+1)
        // with offset solved.
        //
        KPirTree::CalcDiscount(t, mIRDim, irZt, mDiscountIR);

        // 
        // Check forward spread can not be negative
        //
        FwdSprd = GetForwardSpread(t);
        if (FwdSprd < 0)
        {
            throw KFailure("%s: forward credit spread (%s) < 0.\n",
                           routine,
                           GtoFormatDate(currDate));
        }

              
        //
        // Calibrate spread drift
        //
        SolveSpreadOffset(
                t,
                nDim,
                StatePrRisky,
                DiscountIR,
                riskyZeroPrice,
                &crZt);
        
        mCRTpZCenter[t] = crZt;


        //
        // Recompute the exact CR discount factor slice between (t and t+1)
        // with offset solved.
        //
        CalcCRSpreadDiscount(t, crZt, mDiscountCR);


        //
        // Compute transition probabilities between t and t+1
        //
        KMrNTree::Update(t);


        //
        // Expand DiscountIR in nDim dimensions,
        // which only varies in the first mIRDim
        // dimensions and is constant in the
        // higher nDim-mIRDim dimensions
        //
        sliceExpand(t, nDim, mIRDim, DiscountIR, DiscountRisky);

        //
        // Combined risky 1-period discount factor
        //
        sliceUnaryOper(DiscountRisky,
                       nDim,
                       mDiscountCR,
                       MULT);
        

        //
        // Compute IR/risky state prices at t+1 by
        // forwarding the slice to t+1 (adjoint of discounting)
        //
        sliceFw(StatePrIR,    mIRDim, DiscountIR,    t);
        sliceFw(StatePrRisky, nDim,   DiscountRisky, t);

#ifndef __NO_CALIB__
        //
        // Test for correct zero coupon price:
        // sum of all state prices should equal zero bond
        //
        sliceSpecialOper(StatePrIR, mIRDim, "sum", &p);

        if (fabs (p - irZeroPrice) > STATE_PRICE_TOL) {
            throw KFailure("%s: at time point %d, %d-D IR "
                "state prices don't add up to Z (res=%lf).\n",
                routine, t, mIRDim, fabs (p - irZeroPrice));
        }

        sliceSpecialOper(StatePrRisky, nDim, "sum", &p);
        if (fabs (p - riskyZeroPrice) > STATE_PRICE_TOL) {
            throw KFailure("%s: at time point %d, EDev %d-D risky state "
                                   "prices don't add up to Z (res=%lf).\n",
                                   routine, t, nDim, fabs (p - riskyZeroPrice));
        }

#endif


        if (debugLevel > DEBUG_LEVEL_TIMELINE) {
            dppLog << format("%s[%4d]:  IR_drift=%14.10f,   CR_dirft=%14.10f "
                             "Z_ir=%12.8f, Z_risky=%12.8f\n",
                        GtoFormatDate(currDate), 
                        t, mTpZCenter[t], mCRTpZCenter[t], 
                        irZeroPrice, riskyZeroPrice);  
        }



        //
        // Set the current time point to t+1 after the forward
        //
        tpIdxCurrent = t+1;


        //
        // Express DEV tool: store state prices in full nDim dimension.
        //
        if (mDevOn && (iterDev != mDevDates.end()))
        {
            //
            // Store state price at special dev dates
            //
            if ( *iterDev == TPDates[t+1])
            {
                devStatePr = sliceNew(nDim);


                // Store corresponding state prices
                if(nDim > mIRDim)
                    sliceUnaryOper(devStatePr,
                                   nDim,
                                   StatePrRisky,
                                   COPY);
                else    // nDim=mIRDim
                    sliceUnaryOper(devStatePr,
                                   mIRDim,
                                   StatePrIR,
                                   COPY);

                mDevStatePr.insert(KMap(TDate, double*)::value_type(
                          *iterDev, devStatePr));

                devStatePr = NULL;

                iterDev++;
            }   // if iterDev
        }   // mDevOn

    }  // for t


    // Free tmp memory
    sliceDelete(StatePrRisky);
    sliceDelete(StatePrIR);
    sliceDelete(DiscountRisky);


    }
    catch (KFailure) {
    // Free tmp memory
    sliceDelete(StatePrRisky);
    sliceDelete(StatePrIR);
    sliceDelete(DiscountRisky);


    throw KFailure("%s: failed.\n", routine);
    }
}
 




//--------------------------------------------------------------
// CMCDS and IR rate are treated differently based on curve idx:
//
// For CMCDS, protection for the period is inserted in the bank,
// in addition to the risky annuity zeros. 
// 
//
TDate
KCrxTree::Insert(const KRateReset &rateReset, bool isCrit)
{
static  char    routine[] = "KCrxTree::Insert(RateReset)";

        int        curveIdx;
        int        nDim = mIRDim + mCRDim;
        TDate      resetDate,                   // reset date
                   resetEffDate,                // rate effective date
                   matDate;

        TDate      modelResetDate;

        KRateReset *irRateReset = NULL;

 try{


        // Fixed rate does NOT have a valid curve name
        // associated with it.  GetCurveIdx() would fail 
        // for fixed rate.

        if (!rateReset.Rate().IsFloating())
                curveIdx = KV_DIFF;
        else
                curveIdx = GetCurveIdx(rateReset.Rate().CurveName());
        

        // Insert zero - both riskless or risky zero to zerobank
        // Only curve index is required, no knowledge of dynamics.
        // TSliceDev will be responsible for risky dev.
        //
        modelResetDate = KPirTree::Insert(rateReset, isCrit);

        // Insert protection leg in the bank
        //
        if (curveIdx >= KV_CREDIT_RISKY)
        {
            resetDate    = rateReset.ResetDate();
            if (resetDate >= mTodayDate)
            {

                resetEffDate = rateReset.EffDate();
                matDate      = resetEffDate + rateReset.Rate().Maturity();

                // Insert two periods between reset, resetEffDate, 
                // and MatDate in the credit bank to be
                // used to calculate protection leg
                // at reset date

	        if (resetEffDate != resetDate)
            	    GetZeroBank(KV_CREDIT_PROT_DEFRECOV).InsertDates(
                                                                 resetEffDate,
               		                                         resetDate);

                GetZeroBank(KV_CREDIT_PROT_DEFRECOV).InsertDates(matDate,
                                       	                         resetDate);

                Insert(resetDate);
                Insert(resetEffDate);
                Insert(matDate);
            }   
        }

        return modelResetDate;

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}


//--------------------------------------------------------------
// CMCDS and IR rate are treated differently based on curve idx:
// Needs extra zero with maturity at "endDate + rateTenor"
// to zero and protection bank.  This represents the last possible
// rate reset for american type of excersize, and allow us to 
// compute the rate at any reset date between resetDate and endDate,
// with interpolation if needed.
// 
//
TDate
KCrxTree::Insert(const KRateReset &rateReset, TDate endDate, bool isCrit)
{
static  char    routine[] = "KCrxTree::Insert(RateReset, EndDate)";

        int        curveIdx;
        int        nDim = mIRDim + mCRDim;
        TDate      resetDate,                   // reset date
                   resetEffDate,                // rate effective date
                   matDate,
                   lastMatDate;

        TDate      modelResetDate;

        KRateReset *irRateReset = NULL;

 try{


        // Fixed rate does NOT have a valid curve name
        // associated with it.  GetCurveIdx() would fail 
        // for fixed rate.

        if (!rateReset.Rate().IsFloating())
                curveIdx = KV_DIFF;
        else
                curveIdx = GetCurveIdx(rateReset.Rate().CurveName());
        

        // Insert zero - both riskless or risky zero to zerobank
        // Only curve index is required, no knowledge of dynamics.
        // TSliceDev will be responsible for risky dev.
        //
        modelResetDate = KPirTree::Insert(rateReset, endDate, isCrit);

        // Insert protection leg in the bank
        //
        if (curveIdx >= KV_CREDIT_RISKY)
        {
            resetDate    = rateReset.ResetDate();
            if (resetDate >= mTodayDate)
            {
                resetEffDate = rateReset.EffDate();
                matDate      = resetEffDate + rateReset.Rate().Maturity();

                lastMatDate  = endDate 
                             + rateReset.Rate().SpotOffset()
                             + rateReset.Rate().Maturity();

                // Insert two periods between reset, resetEffDate, 
                // and MatDate in the credit bank to be
                // used to calculate protection leg
                // at reset date

	        if (resetEffDate != resetDate)
            	    GetZeroBank(KV_CREDIT_PROT_DEFRECOV).InsertDates(
                                                                 resetEffDate,
               		                                         resetDate);

                GetZeroBank(KV_CREDIT_PROT_DEFRECOV).InsertDates(matDate,
                                       	                         resetDate);

                GetZeroBank(KV_CREDIT_PROT_DEFRECOV).InsertDates(lastMatDate,
                             	                                 resetDate);
                Insert(resetDate);
                Insert(resetEffDate);
                Insert(matDate);
                Insert(lastMatDate);
            }   
        }

        return modelResetDate;

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}




//--------------------------------------------------------------
// Get the zero reset. Discount curve can only be diffuse, or
// deterministic spread curve
//
void
KCrxTree::Get(KTSlice &ts, const KZeroReset &zeroReset)
{
static  char        routine[] = "KCrxTree::Get(ZeroReset)";

        int     curveIdx;

        // Curve for zero reset has to be either diffuse
        // or deterministic spread curve
        //
        // Risky discount and default probability are allowed
        //
        curveIdx = GetCurveIdx(zeroReset.mCurveName);
        if (curveIdx >= KV_BASIS &&
            curveIdx < KV_CREDIT_RISKY)
                throw KFailure("%s: Discount curve can NOT be basis curve.\n",
                                routine);

        KPirTree::Get(ts, zeroReset);

}



//--------------------------------------------------------------
// Called AFTER tree update.  Check the curve type:
// 1. If curve is NON-CREDIT, then use KPirTree:Get
// 2. If curve is KV_CREDIT_RISKY, then compute the rate from
//   
// Stub is needed if reset date is in the past, and today falls 
// between resetDate and effAEDate.  The reset value is contained in
// the reset bank.
//
// For basis rate or spread, the reset value in the bank is assumed
// to be the average of past discret reset rates.  The result is a 
// weighted average of past reset rate and future rate.
//
// The dimension of rateTS is always nDim.
//
void
KCrxTree::Get(KTSlice &rateTS, const KRateReset &rateReset)
{
static  char    routine[] = "KCrxTree::Get(RateReset)";
        
        int        t = TPIdxCurrent();
        TDate      resetDate, effDate, matDate;
        TDate      valueDate;
        TDate      currDate = TPDateCurrent();

        double     cdsResetValue;
        double     resetFwd, currFwd;

        int        nDim = mIRDim + mCRDim;                // The full dimension

        String     rateTSCVName;

        int        curveIdx,
                   irDiscountIdx;


try {

        // Fixed rate does NOT have a valid curve name
        // associated with it.  GetCurveIdx() would fail 
        // for fixed rate.

        if (!rateReset.Rate().IsFloating())
                curveIdx = KV_DIFF;   // so to fall under the pir tree category
        else
                curveIdx = GetCurveIdx(rateReset.Rate().CurveName());
        
        resetDate = rateReset.ResetDate();
        effDate   = rateReset.EffDate();
        matDate   = effDate + rateReset.Rate().Maturity();

        valueDate = GetValueDate(curveIdx);

        // Save the curve name if available
        rateTSCVName = rateTS.GetCurveName();

        // Set the time point of the slice
        rateTS.SetTpIdx(t);

        // check rate type
        if(curveIdx < KV_CREDIT_RISKY)    // IR rate
        {

            KPirTree::Get(rateTS, rateReset);

        }
        else if (curveIdx == KV_CREDIT_RISKY)   // CDS 
        {
            KRateReset rateResetNS = rateReset;
            rateResetNS.Rate().SetSpread(0e0);
    
 
            //
            //
            // Get from reset bank if has already been manually
            // reset and past accural end date
            //
            if(mResetBank->Get(rateResetNS, &cdsResetValue))
            {
                    rateTS = cdsResetValue;
            }
            else 
            {
                if (resetDate < mTodayDate) {
                        throw KFailure("%s: (today %s) no available reset "
                                "value on %s.\n", routine,
                                GtoFormatDate(mTodayDate),
                                GtoFormatDate(resetDate));

                }

                // Check slice dimension for consistency
                if(rateTS.GetSliceDim() != nDim)
                    throw KFailure("%s: dimension of cds spread slice "
                                   "(%s) %d != total nDim %d.\n", 
                                   routine, 
                                   rateTS.GetSliceName().c_str(),
                                   rateTS.GetSliceDim(),
                                   nDim);

                // Temp slice for annuity
                KTSlice annuity(*this,
                                format("Annuity reset on %s",
                                        GtoFormatDate(resetDate)),
                                rateTSCVName);

                KTSlice protTmpSlice(*this,
                                     format("Tmp protection slice on %s",
                                             GtoFormatDate(resetDate)),
                                     rateTSCVName);

                //
                // Reset strictly in the past, need to adjust by the forward
                // ratio.
                // If resetDate < currDate <= resetEffDate,
                // use (currDate, resetEFfDate) as reset and reset effective
                // dates.  This would use correct zeros to produce
                // correct rate, rather than approximation.
                //
                if (currDate <= resetDate)
                {
                    // curveIdx == KV_CREDIT_RISKY
                    GetZeroBank(curveIdx).ComputeParYieldAndAnnuity(t,
                                                                rateReset,
                                                                rateTS,
                                                                annuity);

                }
                else if (currDate <= effDate)
                {
                    // Same KRate reset at (currDate, effDate)
                    //
                    KRateReset currRateReset(rateReset.Rate().CurveName(),
                                             currDate,
                                             effDate,
                                             rateReset.Rate());
                    GetZeroBank(curveIdx).ComputeParYieldAndAnnuity(t,
                                                                currRateReset,
                                                                rateTS,
                                                                annuity);
                }
                else  // (effDate < currDate). truely reset in the past
                {
                    // Same KRate reset at current date
                    //
                    KRateReset approxRateReset(
                                        rateReset.Rate().CurveName(),
                                        currDate,
                                        currDate+rateReset.Rate().SpotOffset(),
                                        rateReset.Rate());


                    GetZeroBank(curveIdx).ComputeParYieldAndAnnuity(t,
                                                                approxRateReset,
                                                                rateTS,
                                                                annuity);


                    // IR discount curve
                    irDiscountIdx = GetCurveIdx(mIRDiscCVName);
                    resetFwd = rateReset.Rate().ForwardCDS(
                                                  GetKZCurve(irDiscountIdx),
                                                  GetKZCurve(curveIdx),
                                                  mRecovery,
                                                  resetDate);

                    currFwd  = rateReset.Rate().ForwardCDS(
                                                  GetKZCurve(irDiscountIdx),
                                                  GetKZCurve(curveIdx),
                                                  mRecovery,
                                                  currDate);

                    // Scale the forward by ratio.
                    annuity *= (currFwd/resetFwd);

                    // change the maturity date of index for
                    // protection leg calculation.
                    //
                    matDate  = currDate 
                             + rateReset.Rate().SpotOffset()
                             + rateReset.Rate().Maturity();
                    
                }

                //
                // Protection leg from [t, mat]
                // Reuse and overwrite the rateTS slice for protection 
                //
                GetZeroBank(KV_CREDIT_PROT_DEFRECOV).GetZeroSlice(
                                                             rateTS, 
                                                             matDate);

    
		// Substract protection leg from [t, resetEff]
		// if needed.
		if (effDate > currDate)
		{
                    GetZeroBank(KV_CREDIT_PROT_DEFRECOV).GetZeroSlice(
						protTmpSlice, 
                                                effDate);
                    rateTS -= protTmpSlice;
		}


                // Recovery rate 1-R.  Recovery is now included
                // in the protection DEV
                // CDS spread = Prot/Annuity
                rateTS /= annuity;

            }    // end of computing rateTS


            //-------------------------------------------------
            // !!!!  WARNING   !!!!!   WARNING !!!
            // We finally add the spread
            //-------------------------------------------------
            rateTS += rateReset.Rate().Spread();

        }  // end of if reset is true
        else
            throw KFailure("%s: invalid rate curve name (%s).\n",
                           routine,
                           rateReset.Rate().CurveName().c_str());
               

        // Restore the slice curve name if is NON-default
        if (rateTSCVName != K_DEFAULT_NAME)
            rateTS.SetCurveIdx(rateTSCVName);


        // Debug print for the rate 
        if(debugLevel > DEBUG_LEVEL_GET) {
                dppLog << format("==================== %s ==================", 
                                 routine) << endl;
                dppLog << endl;
                dppLog << rateTS << endl;
        }

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}






//--------------------------------------------------------------
//
 
double
KCrxTree::GetForwardSpread(int t)
{
static  char    routine[] = "KCrxTree::GetForwardSpread";
 
        // Credit spread forward is computed and stored as
        // part of IR curve initialization
        //
        KMap(int, double*)::iterator it = mTpFRates.find(KV_CREDIT_RISKY);
 
        if (it == mTpFRates.end())
                throw KFailure("Forward rate for credit curve [%d] "
                               "does not exist.\n",
                               routine, KV_CREDIT_RISKY);
 
        return (*it).second[t];
 
}





//--------------------------------------------------------------
//
double
KCrxTree::GetCRSpreadVolBbq(int t)
{
static  char    routine[] = "KCrxTree::GetCRSpreadVolBbq";
 
        TDate   currDate;

        double  QMid;
 
        double  fwdSpread;        // Fwd spread 
        double  FwdSprdA;        // Fwd spread adjusted
 
        double  VolBbq;
 

        currDate = TPDate(t);
 
        // Forward spread
        //
        fwdSpread = GetForwardSpread(t);

        // Avoid zero spread, which is allowed in normal case but
        // canceled out after the mapping.  Set the minimum to 0.01bp
        SHIFT_ZERO(fwdSpread);

        //
        // Check on mCRFSh non-singular, i.e. mCRFSh != -1,
        // QMid * mCRFSh != -1, and QRight * mCRFSh > -1, etc.
        // are done in modpar.cxx
        //
        QMid = (mCRQLo + mCRQHi) / 2;

        FwdSprdA  = fwdSpread / (1. + mCRFSh);


        VolBbq  = mCRVolLogn * mCRBbq
                + mCRVolNorm * (1. - mCRBbq) / fwdSpread;
 
        VolBbq *= (1. + mCRFSh) / (1. + QMid * mCRFSh);
 
        if (debugLevel >= DEBUG_LEVEL_DRIFT) {
                dppLog << format("%s: TPIDX=%4d VolBbq=%14.10f\n",
                        routine, t, VolBbq);
        }

        return VolBbq;
}





//--------------------------------------------------------------
// Solves for the offset to calibrate the zero price at tpIdx+1.
//

void
KCrxTree::SolveSpreadOffset(
        int tpIdx,                // (I) time point index
        int nDim,                 // (I) Spread dimensions
        double *StatePrRisky,     // (I) slice risky state price at t
        double *DiscountIR,       // (I) IR discount from t to t+1
        double zeroPrice,         // (I) zero price at t+1
        double *crZt)             // (I/O) spread shift at t 

{
static        char        routine[] = "KCrxTree::SolveSpreadOffset";

        int     t = tpIdx,              // more convenient ...
                i,               // dimension 1, 2, 3 indices
                j0, j1,
                Mid;

        TDate   currDate = TPDate(tpIdx);

        double  du;
        double  pJump00,
                pJump10, pJump11,
                pJump20, pJump21, pJump22;

        double  QSh;            // Used to calculate initial Zt
        double  baseSh;         // Base shift with distribution
                                // at X=0 corresponds to forward spread

        double  xSwitch;        // Actual 2q switch point in X spac

        double  QLeft    = this->mCRQLo;
        double  QRight   = this->mCRQHi;
        double  FwdShift = this->mCRFSh;

        double  FwdSprd;        // Fwd spread 
        double  FwdSprdA;       // Fwd spread adjusted 
        double  VolBbq;         // Sigma used in bone mapping

        double  Zidx;           // Zt index i,j,k adjusted

        double  MLo, MHi,       // constants used in grid calc
                SLo, SHi,
                Grid,           // grid step
                dGrid,
                discount,       // 1-period discount
                DQ,
                rJumpHi,        // rate jump size high and low
                rJumpLo,
                rJumpHi2,       // rate jump size high and low
                rJumpLo2,
                rJumpHi5,       // rate jump size high and low
                rJumpLo5;



        int     *Top1         = mTop1;
        int     *Bottom1      = mBottom1;
        int     **Top2        = mTop2;
        int     **Bottom2     = mBottom2;
        int     ***Top3       = mTop3;
        int     ***Bottom3    = mBottom3;
        int     *OutTop1      = mOutTop1;
        int     *OutBottom1   = mOutBottom1;
        int     **OutTop2     = mOutTop2;
        int     **OutBottom2  = mOutBottom2;
        int     ***OutTop3    = mOutTop3;
        int     ***OutBottom3 = mOutBottom3;


        double  *statePr  = NULL,
                *statePrL = NULL;

        double  eProt;          // expected value of 1-period protection
        double  crZt0 = *crZt;  // Initial guestt

        double  F, dF;          // Root finding of F and dF/dx.
        double  Eta;            // Change of crZt

        int     iterNR;

        const   int     MAXITER = 100;
        const   double  MAXXERR = 1e-7;


    try {

        // StatePrRisky * DiscountIR
        //
        statePr = sliceNew(nDim);


        //
        // Expand DiscountIR in nDim dimensions,
        // which only varies in the first mIRDim
        // dimensions and is constant in the
        // higher nDim-mIRDim dimensions
        //
        sliceExpand(t, nDim, mIRDim, DiscountIR, statePr);


        sliceUnaryOper( statePr,
                        nDim,
                        StatePrRisky,
                        MULT);


        //
        // Precompute jumps 
        //

        du = sqrt (JUMPCOEFF * LengthJ[t-1]);
        
        switch (nDim) {
        case 3:
                pJump20 = mAweight[t-1][3] * du;
                pJump21 = mAweight[t-1][4] * du;
                pJump22 = mAweight[t-1][5] * du;
        case 2:
                pJump10 = mAweight[t-1][1] * du;
                pJump11 = mAweight[t-1][2] * du;
        case 1:
                pJump00 = mAweight[t-1][0] * du;
        }


        //
        // Calc needed constants
        //
        FwdSprd = GetForwardSpread(t);

        // Special case for zero spread, which needs no calibration
        //   
        if (IS_ALMOST_ZERO(FwdSprd))
        {
            goto done;
        }

        // Avoid special treatment for normal.
        SHIFT_ZERO(QLeft);
        SHIFT_ZERO(QRight);


        FwdSprdA  = FwdSprd / (1. + FwdShift);


        // Vol backbone factor
        //
        VolBbq = GetCRSpreadVolBbq(t);


        //
        // The base shift in the spread distribution is defined in
        // such a way that spread distribution at X=0 ALWAYS corresponds
        // to the forward spread.
        //
        QSh = (FwdShift > 0) ? QRight : QLeft;
        if (fabs(QSh) > QCUTOFF)
        {
                baseSh = log(1. + QSh * FwdShift) / (QSh * VolBbq);
        }
        else
        {
                baseSh = FwdShift / VolBbq;
        }

        // Actual 2q switch point in the X space
        // This will be used as initial guess for Zt
        xSwitch = crZt0 + baseSh;


        // Calculate the expected value of 1-period protection
        eProt = zeroPrice * (1. - 1. / (1.+FwdSprd));



        //
        // Solve for the offset shift Zt for basis spread
        // Using Newton-Ralphson method
        //

        for (iterNR = 0; iterNR < MAXITER; iterNR++)
        {

            //
            // Compute F(X) and F'(X) given xSwitch
            //

            // Reset F and F'
            //
            F = dF = 0e0;

            MLo = FwdSprdA / QLeft;
            MHi = FwdSprdA / QRight;
            SLo = 1. + FwdSprdA - FwdSprdA / QLeft;
            SHi = 1. + FwdSprdA - FwdSprdA / QRight;
            DQ  = FwdSprdA * VolBbq;

            switch (nDim) {
            case 1:
                rJumpLo = exp(QLeft  * VolBbq * pJump00);
                rJumpHi = exp(QRight * VolBbq * pJump00);

                // 1-D calibration
                Zidx = xSwitch;
                Mid  = (int) ceil(-Zidx / pJump00) - 1;
                Mid  = Min( Max(Mid, mBottom1[t] - 1), mTop1[t]);

                // Compute coefficients

                statePrL  = statePr  + NodeOffset(1, 0, 0, t);

                // LEFT part of distribution
                Zidx += (pJump00) * mBottom1[t];
                Grid  = MLo * exp (QLeft * VolBbq * Zidx);
                dGrid = DQ  * exp (QLeft * VolBbq * Zidx);

                for (j0 = mBottom1[t]; j0 <= Mid; j0++) {
                        discount = 1. / (SLo + Grid);

                        Grid  *= rJumpLo;
                        dGrid *= rJumpLo;

                        F  += statePrL[j0] * discount;
                        dF += statePrL[j0] * dGrid * discount * discount;

                }

                // RIGHT part of distribution
                Zidx  = xSwitch + (pJump00) * (Mid + 1);
                Grid  = MHi * exp (QRight * VolBbq * Zidx);
                dGrid = DQ  * exp (QRight * VolBbq * Zidx);

                for (j0 = Mid + 1; j0 <= mTop1[t]; j0++) {
                        discount = 1. / (SHi + Grid);

                        Grid  *= rJumpHi;
                        dGrid *= rJumpHi;

                        F  += statePrL[j0] * discount;
                        dF += statePrL[j0] * dGrid * discount * discount;
                }

                break;

            case 2:   
                rJumpLo2 = exp (QLeft  * VolBbq * pJump11);
                rJumpHi2 = exp (QRight * VolBbq * pJump11);

                switch (mCRDim) {
                case 1:  // mIRDim = 1, mCRDim = 1;  THIS IS THE MODE SUPPORTED

                    for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                    {
                        // Find switch point index
                        Zidx  = xSwitch + j0 * pJump10;
                        Mid   = (int) ceil(-Zidx / pJump11) - 1;
                        Mid   = Min(Max(Mid, mBottom2[t][j0] - 1),
                                        mTop2[t][j0]);

                        statePrL  = statePr  + NodeOffset(2, j0, 0, t);

                        // LEFT part of distribution
                        Zidx += (pJump11) * mBottom2[t][j0];
                        Grid  = MLo * exp (QLeft * VolBbq * Zidx);
                        dGrid = DQ  * exp (QLeft * VolBbq * Zidx);

                        for (i = mBottom2[t][j0]; i <= Mid; i++) {
                                discount = 1. / (SLo + Grid);

                                Grid  *= rJumpLo2;
                                dGrid *= rJumpLo2;

                                F  += statePrL[i] * discount;
                                dF += statePrL[i] * dGrid * discount * discount;
                        }


                        // RIGHT part of distribution
                        Zidx  = xSwitch + j0 * pJump10 + (Mid + 1) * pJump11;
                        Grid  = MHi * exp (QRight * VolBbq * Zidx);
                        dGrid = DQ  * exp (QRight * VolBbq * Zidx);

                        for (i = Mid + 1; i <= mTop2[t][j0]; i++) {
                                discount = 1. / (SHi + Grid);

                                Grid  *= rJumpHi2;
                                dGrid *= rJumpHi2;

                                F  += statePrL[i] * discount;
                                dF += statePrL[i] * dGrid * discount * discount;
                        }
                    }  // for j0

                    break;

                case 2:  // mIRDim = 0, mCRDim = 2;
                    for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                    {
                        // Find switch point index
                        Zidx  = xSwitch + j0 * (pJump00 + pJump10);
                        Mid   = (int) ceil(-Zidx / pJump11) - 1;
                        Mid   = Min(Max(Mid, mBottom2[t][j0] - 1),
                                        mTop2[t][j0]);

                        statePrL  = statePr  + NodeOffset(2, j0, 0, t);


                        // LEFT part of distribution
                        Zidx += (pJump11) * mBottom2[t][j0];
                        Grid  = MLo * exp (QLeft * VolBbq * Zidx);
                        dGrid = DQ  * exp (QLeft * VolBbq * Zidx);

                        for (i = mBottom2[t][j0]; i <= Mid; i++) {
                                discount = 1. / (SLo + Grid);

                                Grid  *= rJumpLo2;
                                dGrid *= rJumpLo2;

                                F  += statePrL[i] * discount;
                                dF += statePrL[i] * dGrid * discount * discount;
                        }


                        // RIGHT part of distribution
                        Zidx  = xSwitch + j0 * (pJump00 + pJump10)
                                        + (Mid + 1) * pJump11;
                        Grid  = MHi * exp (QRight * VolBbq * Zidx);
                        dGrid = DQ  * exp (QRight * VolBbq * Zidx);

                        for (i = Mid + 1; i <= mTop2[t][j0]; i++) {
                                discount = 1. / (SHi + Grid);

                                Grid  *= rJumpHi2;
                                dGrid *= rJumpHi2;

                                F  += statePrL[i] * discount;
                                dF += statePrL[i] * dGrid * discount * discount;
                        }

                    }

                    break;

                default:
                    throw KFailure("%s: invalid factor dimensions."
                               " nDim=%d, nCRDim=%d.\n", routine, nDim, mCRDim);
                }

                break;

            case 3:

                rJumpLo5 = exp (QLeft  * VolBbq * pJump22);
                rJumpHi5 = exp (QRight * VolBbq * pJump22);

                switch (mCRDim) {
                case 1:     // mIRDim = 2, nCRDim = 1;

                    for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                    for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
                    {
                        Zidx  = xSwitch
                              + j0 * (pJump20)
                              + j1 * (pJump21);
                        Mid   = (int) ceil(-Zidx / pJump22) - 1;
                        Mid   = Min(Max(Mid, mBottom3[t][j0][j1] - 1),
                                        mTop3[t][j0][j1]);

                        statePrL  = statePr  + NodeOffset(3, j0, j1, t);

                        // LEFT part of distribution
                        Zidx += (pJump22) * mBottom3[t][j0][j1];
                        Grid  = MLo * exp (QLeft * VolBbq * Zidx);
                        dGrid = DQ  * exp (QLeft * VolBbq * Zidx);

                        for (i = mBottom3[t][j0][j1]; i <= Mid; i++) {
                                discount = 1. / (SLo + Grid);

                                Grid  *= rJumpLo5;
                                dGrid *= rJumpLo5;

                                F  += statePrL[i] * discount;
                                dF += statePrL[i] * dGrid * discount * discount;
                        }

                        // RIGHT part of distribution
                        Zidx = xSwitch
                             + j0 * (pJump20)
                             + j1 * (pJump21)
                             + (Mid+1) * (pJump22);
                        Grid  = MHi * exp (QRight * VolBbq * Zidx);
                        dGrid = DQ  * exp (QRight * VolBbq * Zidx);

                        for (i = Mid + 1; i <= mTop3[t][j0][j1]; i++) {
                                discount = 1. / (SHi + Grid);

                                Grid  *= rJumpHi5;
                                dGrid *= rJumpHi5;

                                F  += statePrL[i] * discount;
                                dF += statePrL[i] * dGrid * discount * discount;
                        }
                    }

                    break;

                case 2:     // mIRDim = 1, nCRDim = 2;

                    for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                    for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
                    {
                        Zidx  = xSwitch
                              + j0 * (pJump10 + pJump20)
                              + j1 * (pJump11 + pJump21);
                        Mid   = (int) ceil(-Zidx / pJump22) - 1;
                        Mid   = Min(Max(Mid, mBottom3[t][j0][j1] - 1),
                                        mTop3[t][j0][j1]);

                        statePrL  = statePr  + NodeOffset(3, j0, j1, t);

                        // LEFT part of distribution
                        Zidx += (pJump22) * mBottom3[t][j0][j1];
                        Grid  = MLo * exp (QLeft * VolBbq * Zidx);
                        dGrid = DQ  * exp (QLeft * VolBbq * Zidx);

                        for (i = mBottom3[t][j0][j1]; i <= Mid; i++) {
                                discount = 1. / (SLo + Grid);

                                Grid  *= rJumpLo5;
                                dGrid *= rJumpLo5;

                                F  += statePrL[i] * discount;
                                dF += statePrL[i] * dGrid * discount * discount;
                        }

                        // RIGHT part of distribution
                        Zidx = xSwitch
                             + j0 * (pJump10 + pJump20)
                             + j1 * (pJump11 + pJump21)
                             + (Mid+1) * (pJump22);
                        Grid  = MHi * exp (QRight * VolBbq * Zidx);
                        dGrid = DQ  * exp (QRight * VolBbq * Zidx);

                        for (i = Mid + 1; i <= mTop3[t][j0][j1]; i++) {
                                discount = 1. / (SHi + Grid);

                                Grid  *= rJumpHi5;
                                dGrid *= rJumpHi5;

                                F  += statePrL[i] * discount;
                                dF += statePrL[i] * dGrid * discount * discount;
                        }
                    }

                    break;

                case 3:     // mIRDim = 0, nCRDim = 3;
                    for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                    for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
                    {
                        // Find switch point index
                        Zidx  = xSwitch
                              + j0 * (pJump00 + pJump10 + pJump20)
                              + j1 * (pJump11 + pJump21);
                        Mid   = (int) ceil(-Zidx / pJump22) - 1;
                        Mid   = Min(Max(Mid, mBottom3[t][j0][j1] - 1),
                                        mTop3[t][j0][j1]);

                        statePrL  = statePr  + NodeOffset(3, j0, j1, t);

                        // LEFT part of distribution
                        Zidx += (pJump22) * mBottom3[t][j0][j1];
                        Grid  = MLo * exp (QLeft * VolBbq * Zidx);
                        dGrid = DQ  * exp (QLeft * VolBbq * Zidx);

                        for (i = mBottom3[t][j0][j1]; i <= Mid; i++) {
                                discount = 1. / (SLo + Grid);

                                Grid  *= rJumpLo5;
                                dGrid *= rJumpLo5;

                                F  += statePrL[i] * discount;
                                dF += statePrL[i] * dGrid * discount * discount;

                        }

                        // RIGHT part of distribution
                        Zidx = xSwitch
                             + j0 * (pJump00 + pJump10 + pJump20)
                             + j1 * (pJump11 + pJump21)
                             + (Mid+1) * (pJump22);
                        Grid  = MHi * exp (QRight * VolBbq * Zidx);
                        dGrid = DQ  * exp (QRight * VolBbq * Zidx);

                        for (i = Mid + 1; i <= mTop3[t][j0][j1]; i++) {
                                discount = 1. / (SHi + Grid);

                                Grid  *= rJumpHi5;
                                dGrid *= rJumpHi5;

                                F  += statePrL[i] * discount;
                                dF += statePrL[i] * dGrid * discount * discount;
                        }
                    }   // for j0, j1

                    break;

                default:
                    throw KFailure("%s: invalid factor dimensions."
                               " nDim=%d, nCRDim=%d.\n", routine, nDim, mCRDim);
                }

                break;    // nDim = 3

            default:
                throw KFailure("%s: N/A for nDim=%d.\n", routine, nDim);

            }    // end of switch nDim


            // Add constant terms and flip sign
            //
            //F -= eProt;
            F -= zeroPrice;
            dF *= -1.;

            //
            // Test the relative error
            // This is needed in case of zero fwd rate
            //
            if ( fabs(F) < MAXXERR)    // 1/1000 bps
            {
                *crZt = xSwitch - baseSh;
                goto done;
            }


            // Check derivative
            //
            if (IS_ALMOST_ZERO(dF))
            {
                throw KFailure("%s: derivative in NR solver should not "
                               "vanish at t=%s.\n",
                               routine,
                               GtoFormatDate(currDate));
            }


            Eta = F / dF;
            xSwitch -= Eta;

            //
            // Test the relative error
            //
            if ( fabs(F) < MAXXERR)    // 1/1000 bps
            {
                *crZt = xSwitch - baseSh;
                goto done;
            }



        }    /* for iterNR */


        throw KFailure("%s: failed spread drift calibration on date %s (F=%lf)."
                       " Maximum number of iterations (%d) exceeded.\n",
                       routine,
                       GtoFormatDate(currDate),
                       F,
                       MAXITER);


   done:

        if (debugLevel > DEBUG_LEVEL_DRIFT) {
                dppLog << endl;
                dppLog << "Credit center offset: " <<  *crZt - crZt0 << endl;
                dppLog << endl;
        }


        sliceDelete(statePr);
        statePrL  = NULL;

    }
    catch (KFailure) {
        sliceDelete(statePr);
        statePrL = NULL;
        throw KFailure("%s: failed.\n", routine);
    }

}





//--------------------------------------------------------------
// Calculates the discount factor applicable between
// tpIdx and tpIdx+1.
//

void
KCrxTree::CalcCRSpreadDiscount(
        int tpIdx,                        // (I) time point index
        double Zt,                        // (I) center offset 
        double *Discount)                // (O) discount slice
{
static        char        routine[] = "KCrxTree::CalcCRSpreadDiscount";


        int     t = tpIdx,                // more convenient ...
                i, j, k,                // dimension 1, 2, 3 indices
                Mid;

        int     nDim;                        // (I) CR dimensions

        double  du;
        double  pJump00,
                pJump10, pJump11,
                pJump20, pJump21, pJump22;


        double  QLeft   = this->mCRQLo;
        double  QRight  = this->mCRQHi;
        double  FwdShift = this->mCRFSh;

        double  FwdSprd;        // Fwd rate 
        double  FwdSprdA;        // Fwd rate adjusted 
        double  MLeft, SLeft;        // Multiple coeff and shift for grid pt
        double  MRight, SRight;
        double  VolBbq;                // Sigma used in bone mapping 
        
        double  QSh;                // q parameter for shift adjustment 
        
        double  Zidx;                // Zt index i,j,k adjusted 
        double  Grid;                // Grid points  

        double  RateJumpLeft,        // Rate space jump size for last dimension
                RateJumpRight,
                RateJumpLeft2,
                RateJumpRight2,
                RateJumpLeft5,
                RateJumpRight5;

        int     *Top1         = mTop1;
        int     *Bottom1      = mBottom1;
        int     **Top2        = mTop2;
        int     **Bottom2     = mBottom2;
        int     ***Top3       = mTop3;
        int     ***Bottom3    = mBottom3;
        int     *OutTop1      = mOutTop1;
        int     *OutBottom1   = mOutBottom1;
        int     **OutTop2     = mOutTop2;
        int     **OutBottom2  = mOutBottom2;
        int     ***OutTop3    = mOutTop3;
        int     ***OutBottom3 = mOutBottom3;


        double  *DiscountL;        // local slice pointer

    try {

    
        nDim = mIRDim + mCRDim;


        //
        // Precompute jumps 
        //

        du = sqrt (JUMPCOEFF * LengthJ[t-1]);
        
        switch (nDim) {
        case 3:
                pJump20 = mAweight[t-1][3] * du;
                pJump21 = mAweight[t-1][4] * du;
                pJump22 = mAweight[t-1][5] * du;
        case 2:
                pJump10 = mAweight[t-1][1] * du;
                pJump11 = mAweight[t-1][2] * du;
        case 1:
                pJump00 = mAweight[t-1][0] * du;
        }



        //
        // Calc needed constants
        //
        FwdSprd = GetForwardSpread(t);

        //
        // Special case for zero spread, which needs no calibration
        //   
        if (IS_ALMOST_ZERO(FwdSprd))
        {
            sliceScalarOper(
                    Discount,
                    nDim,
                    1.0,
                    COPY);
            return;    
        }


        FwdSprdA  = FwdSprd / (1. + FwdShift);


        // Vol backbone factor
        //
        VolBbq = GetCRSpreadVolBbq(t);


        // Set initial Zt
        // Make sure that the rate distribution at X=0 
        // corresponds to the forward rate
        if (t == 0)
        {
            QSh = (FwdShift > 0) ? QRight : QLeft;
            if (IS_Q(QSh))
            {
                Zt = log(1. + QSh * FwdShift) / (QSh * VolBbq);
            }
            else
            {
                Zt = FwdShift / VolBbq;
            }
        }

        //
        // Dimension specific 
        //

        //----------------------------------------------
        // 1 F
        //----------------------------------------------

        if (nDim == 1) {

        DiscountL = Discount + NodeOffset(1, 0, 0, t);


        if (IS_Q(QLeft))
        {
            MLeft        = FwdSprdA / QLeft;
            SLeft        = 1 + FwdSprdA - FwdSprdA / QLeft;     // 1 + r 
            RateJumpLeft = exp(QLeft * VolBbq * pJump00);
        }
        else
        {
            MLeft        = FwdSprdA;
            SLeft        = 1 + FwdSprdA;
            RateJumpLeft = FwdSprdA * VolBbq * pJump00;
        }
        if (IS_Q(QRight))
        {
            MRight        = FwdSprdA / QRight;
            SRight        = 1 + FwdSprdA - FwdSprdA / QRight;  // 1 + r      
            RateJumpRight = exp(QRight * VolBbq * pJump00);
        }
        else
        {
            MRight        = FwdSprdA;
            SRight        = 1 + FwdSprdA;
            RateJumpRight = FwdSprdA * VolBbq * pJump00;
        }


        /* 
        *   Set up the grid points. 
        */

        /* LEFT part of distribution */

        Zidx  = Zt; 
        Mid   = (int) ceil(-Zidx / pJump00) - 1;
        Mid   = MIN ( MAX (Mid, Bottom1[t] - 1), Top1[t]);
        Zidx += (pJump00) * Bottom1[t];



        if (IS_Q(QLeft))
        {
            Grid = MLeft * exp (QLeft * VolBbq * Zidx);
        
            for (i = Bottom1[t]; i <= Mid; i++)
            {
                DiscountL[i] = 1. / (SLeft + Grid);
                
                Grid *= RateJumpLeft;
            }           
        }
        else
        {
            Grid = MLeft * VolBbq * Zidx;

            for (i = Bottom1[t]; i <= Mid; i++)
            {
                DiscountL[i] = 1. / (SLeft + Grid);
                
                Grid += RateJumpLeft;
            }
        }

        /* Right part of distribution */

        Zidx = Zt 
             + (pJump00) * (Mid + 1);

        if (IS_Q(QRight))
        {
            Grid = MRight * exp (QRight * VolBbq * Zidx);
        
            for (i = Mid + 1; i <= Top1[t]; i++)
            {
                DiscountL[i] = 1. / (SRight + Grid);
                
                Grid *= RateJumpRight;
            }           
        }
        else
        {
            Grid = MRight * VolBbq * Zidx;

            for (i = Mid +1; i <= Top1[t]; i++)
            {
                DiscountL[i] = 1. / (SRight + Grid);
                
                Grid += RateJumpRight;
            }
        }


        //----------------------------------------------
        // 2 F
        //----------------------------------------------
        } else if (nDim == 2) {


        if (IS_Q(QLeft))
        {
            MLeft         = FwdSprdA / QLeft;
            SLeft         = 1 + FwdSprdA - FwdSprdA / QLeft;      
            RateJumpLeft2 = exp(QLeft * VolBbq * pJump11);
        }
        else
        {
            MLeft         = FwdSprdA;
            SLeft         = 1 + FwdSprdA;
            RateJumpLeft2 = FwdSprdA * VolBbq * pJump11;
        }
        if (IS_Q(QRight))
        {
            MRight         = FwdSprdA / QRight;
            SRight         = 1 + FwdSprdA - FwdSprdA / QRight;      
            RateJumpRight2 = exp(QRight * VolBbq * pJump11);
        }
        else
        {
            MRight         = FwdSprdA;
            SRight         = 1 + FwdSprdA;
            RateJumpRight2 = FwdSprdA * VolBbq * pJump11;
        }



        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {
            /* LEFT part of distribution */

            Zidx  = Zt + pJump10 * i;
            Mid   = (int) ceil(-Zidx / pJump11) - 1;
            Mid   = MIN ( MAX (Mid, Bottom2[t][i] - 1), Top2[t][i]);
            Zidx += (pJump11) * Bottom2[t][i];

            DiscountL = Discount + NodeOffset (2, i, 0, t);
    
            if (IS_Q(QLeft))
            {
                Grid = MLeft * exp (QLeft * VolBbq * Zidx);
            
                for (j = Bottom2[t][i]; j <= Mid; j++)
                {
                    DiscountL[j] = 1. / (SLeft + Grid);
                
                    Grid *= RateJumpLeft2;
                }
            }
            else
            {
                Grid = MLeft * VolBbq * Zidx;

                for (j = Bottom2[t][i]; j <= Mid; j++)
                {
                    DiscountL[j] = 1. / (SLeft + Grid);
                
                    Grid += RateJumpLeft2;
                }
            }

            /* Right part of distribution */

            Zidx = Zt 
                 + (pJump10) * i
                 + (pJump11) * (Mid + 1);

            if (IS_Q(QRight))
            {
                Grid = MRight * exp (QRight * VolBbq * Zidx);
            
                for (j = Mid + 1; j <= Top2[t][i]; j++)
                {
                    DiscountL[j] = 1. / (SRight + Grid);
                
                    Grid *= RateJumpRight2;
                }
            }
            else
            {
                Grid = MRight * VolBbq * Zidx;

                for (j = Mid + 1; j <= Top2[t][i]; j++)
                {
                    DiscountL[j] = 1. / (SRight + Grid);
                
                    Grid += RateJumpRight2;
                }
            }

        }  /* for i */



        //----------------------------------------------
        // 3 F
        //----------------------------------------------

        } else if (nDim == 3) {


        if (IS_Q(QLeft))
        {
            MLeft         = FwdSprdA / QLeft;
            SLeft         = 1 + FwdSprdA - FwdSprdA / QLeft;      
            RateJumpLeft5 = exp(QLeft * VolBbq * pJump22);
        }
        else
        {
            MLeft         = FwdSprdA;
            SLeft         = 1 + FwdSprdA;
            RateJumpLeft5 = FwdSprdA * VolBbq * pJump22;
        }
        if (IS_Q(QRight))
        {
            MRight         = FwdSprdA / QRight;
            SRight         = 1 + FwdSprdA - FwdSprdA / QRight;      
            RateJumpRight5 = exp(QRight * VolBbq * pJump22);
        }
        else
        {
            MRight         = FwdSprdA;
            SRight         = 1 + FwdSprdA;
            RateJumpRight5 = FwdSprdA * VolBbq * pJump22;
        }



        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {
            for (j = Bottom2[t][i]; j <= Top2[t][i]; j++)
            {
                /* LEFT part of distribution */

                Zidx  = Zt 
                      + (pJump00 + pJump10 + pJump20) * i
                      + (pJump11 + pJump21) * j;
                Mid   = (int) ceil(-Zidx / pJump22) - 1;
                Mid   = MIN ( MAX (Mid, Bottom3[t][i][j] - 1), Top3[t][i][j]);
                Zidx += (pJump22) * Bottom3[t][i][j];
            
                DiscountL = Discount + NodeOffset (3, i, j, t);
                
                if (IS_Q(QLeft))
                {
                    Grid = MLeft * exp (QLeft * VolBbq * Zidx);
    
                    for (k = Bottom3[t][i][j]; k <= Mid; k++)
                    {
                        DiscountL[k] = 1. / (SLeft + Grid);
                
                        Grid *= RateJumpLeft5;
                    }
                }
                else
                {
                    Grid = MLeft * VolBbq * Zidx;

                    for (k = Bottom3[t][i][j]; k <= Mid; k++)
                    {
                        DiscountL[k] = 1. / (SLeft + Grid);
                
                        Grid += RateJumpLeft5;
                    }
                }

                /* Right part of distribution */

                Zidx = Zt 
                     + (pJump00 + pJump10 + pJump20) * i
                     + (pJump11 + pJump21) * j
                     + (pJump22) * (Mid + 1);

                if (IS_Q(QRight))
                {
                    Grid = MRight * exp (QRight * VolBbq * Zidx);
    
                    for (k = Mid + 1; k <= Top3[t][i][j]; k++)
                    {
                        DiscountL[k] = 1. / (SRight + Grid);
                
                        Grid *= RateJumpRight5;
                    }
                }
                else
                {
                    Grid = MRight * VolBbq * Zidx;

                    for (k = Mid + 1; k <= Top3[t][i][j]; k++)
                    {
                        DiscountL[k] = 1. / (SRight + Grid);
                
                        Grid += RateJumpRight5;
                    }
                }

            }  /* for j */     
        
        }  /* for i */
       
 
        } // if (nDim == 
                



    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}

