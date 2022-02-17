/***************************************************************
 * Module:        BasisTree
 * Submodule:        
 * File:        
 * Function:        Basis tree 
 * Author:        David Liu, June 1999
 ***************************************************************/
#include "kstdinc.h"    // Standard definitions & error hadling 
#include "kstlutil.h"        // KCArray 

#define        _kmrntree_SRC
#include "kmrntree.h"
#include "ktmxirtree.h"
#include "kbirtree.h"

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
#include "cfinanci.h"                        // GtoRateToDiscount 
#include "interp.h"                        // GtoRateToDiscount 
#include "zr2simp.h"                        // GtoZerosToSimplePoint 

#include "drlmem.h"                        // DrlDoubleVectAlloc/Free 
#include "drlio.h"                        // 
#include "drltime.h"                        // DrlTDatePrint 
#include "drlvtype.h"                        // DrlVTypeVectAdd 
#include "drlsort.h"                        // DrlTDateArrayFloorIdx() 
#include "drlroot.h"                        // DrlNRPoly
#include "drloptio.h"                       // ConvexityC_BS2Q
#include "drlinter.h"                       // DrlTDateLinearInterp1d

#include "vnfmanly.h"                        // VnfmVolCalib1VArbitrary

};


#ifndef SHIFT_ZERO
#define SHIFT_ZERO(x) ((x) = (IS_ALMOST_ZERO(x) ? 1e-5 : (x)))
#endif


// >=1. Polynomial order for solving drift of spread.
const        int        NPOLY = 2;        

#ifndef QCUTOFF
const        double        QCUTOFF = 1e-5;
#endif


//--------------------------------------------------------------
//

KBirTree::KBirTree() : KTMXirTree()
{
        mBSOn = false;

        mBSTmpRate = NULL;
        mBSSprdOn = false;

        mBSDim = 0;
        mBSType = SUB_SPREAD;        

        // Lognormal
        mBbq   = 1e0;
        mBSQLo = 1e0;
        mBSQHi = 1e0;
        mBSFSh = 0e0;

        mBSDelayShift = 0e0;

        mBSDiscCVName= K_DEFAULT_NAME;
        mLiborCVName = K_DEFAULT_NAME;
        mBasisDCC    = GTO_ACT_365;
        mLiborDCC    = GTO_ACT_360;
        mBasisCVName = K_DEFAULT_NAME;

}




//--------------------------------------------------------------
//

KBirTree::~KBirTree()
{

        delete mBSTmpRate;

        for(KMap(TDate, KRateReset*)::iterator iterRate=mLiborRates.begin();
            iterRate != mLiborRates.end(); ++iterRate) {
                delete (*iterRate).second;
        }

}




//---------------------------------------------------------------
// Initialize basis tree.
 
void
KBirTree::Initialize(
        KMarketCurves&   marketCurves,     // (I) curves and curve types

        KVolDiag&        irVolDiag,        // (I) IR volatility data.
        KMrParam&        irMrParam,        // (I) IR mr data.
        KSmileParam&     irSmileParam,     // (I) IR skew data.
        KVolDiag&        bsVolDiag,        // (I) Basis volatility data.
        KMrParam&        bsMrParam,        // (I) Basis mr data.
        KSmileParam&     bsSmileParam,     // (I) Basis skew data.
        double           irBsCorr,         // (I) IR basis correlation.

        KResetBank       &resetBank)       // (I) rate reset bank
{
static  char            routine[] = "KBirTree::Initialize";
        int             nDim,
                        irDim,
                        bsDim;              // Full dimension of basis tree
        int             idx;
        KMrParam        treeMrParam;        // full tree mr parameters
 try{

        KVector(TDate)           factVolDates;        // spot vol dates
        KVector(KVector(double)) factVolCurves;       // factor spot vol curves


        // Check inputs
        if (!marketCurves.IsValid() ||
            !irVolDiag.IsValid()    ||
            !irMrParam.IsValid()    ||
            !irSmileParam.IsValid() ||
            !bsVolDiag.IsValid()    ||
            !bsMrParam.IsValid()    ||
            !bsSmileParam.IsValid() ||
            !IsConsistent(irVolDiag, bsVolDiag))
                throw KFailure("Invalid market or model inputs.\n");



        //----------------------------------------------
        // (1) Build the full mr model parameters
        //----------------------------------------------
        treeMrParam = Correlate(irMrParam, bsMrParam, irBsCorr);
        nDim  = treeMrParam.mNumFact;
        bsDim = bsMrParam.mNumFact;
        irDim = nDim - bsDim;


        //----------------------------------------------
        // (2) Turn basis flag ON if basis curve
        //     or spread is involved in the product
        //     Insert additional basis curve if only
        //     basis spread curve is specified in the env.
        //----------------------------------------------
 
        mBSOn = marketCurves.IsBasis();


        //----------------------------------------------
        // (3) Perform IR volatility calibration
        //----------------------------------------------

        //
        // Get the diffused curve for the volatility calibration
        //
        KZCurve &diffCurve = marketCurves.GetDiffuse();


        // Calibrate IR spot vols in irDim dimensions
        DppCalibIRVolDiag(
                irMrParam,
                irSmileParam,
                irDim,
                irVolDiag,
                diffCurve,
                factVolDates,
                factVolCurves);



        //----------------------------------------------
        // (4) Assign basis model parameters
        //     if basis curve or spread is involved in the product
        //----------------------------------------------
        if (mBSOn) {
                mBSType = marketCurves.mBSType;

                mBSDelayShift = marketCurves.mBSDelayShift;
                                
                mBSDim = bsDim;
                mBSQLo = bsSmileParam.mQ1;
                mBSQHi = bsSmileParam.mQ2;
                mBSFSh = bsSmileParam.mQF;

                // Basis Backbone parameters
                mBSBbq = bsMrParam.mBackboneQ;
                

                // Store the basis market vol for estimating the 
                // initial drift in the back end
                mBSVolDiag = bsVolDiag;

                // Compute reference vol constants
                double norm = 0e0;
                for (int idx = 0; idx<bsDim; idx++)
                        norm += bsMrParam.mAlpha[idx] * bsMrParam.mAlpha[idx];

                norm = sqrt(norm);
                if (fabs(norm) < DBL_EPSILON)
                        throw KFailure("%s: total alpha too small %g.\n",
                                routine, norm);
 
                if (IS_ALMOST_ZERO(mBSBbq - 1e0)) {
                        mBSVolNorm = 0.;
                        mBSVolLogn = norm;
                } else if (IS_ALMOST_ZERO(mBSBbq - 0e0)) {
                        mBSVolNorm = norm;
                        mBSVolLogn = 0.;
                } else {
                        throw KFailure("%s: backbone must be 0 or 1 (%g).\n",
                                routine, mBbq);
                }

                // Check if mBSDiscCVName and mLiborCVName are in the table.
                // to be done ...
                //
                mBSDiscCVName = marketCurves.mBSDiscCVName;
                mLiborCVName  = marketCurves.mLiborCVName;
                mBasisCVName  = marketCurves.mCVNames[KV_BASIS];
                
                mLiborDCC     = marketCurves.mLiborDCC;
                mBasisDCC     = marketCurves.mBasisDCC;

                //
                // Basis bank for par spread calculation
                //
                mBSBank = KCBank(*this,        
                                  format("CBank For Par Spread Curve %d",
                                        KV_PAR_SPREAD));

                //
                // Calibrate Basis spot vols in bsDim dimensions
                // And add basis spot vol curve to total factor vols
                //
                DppCalibSpreadVolDiag(
                        bsMrParam,
                        bsSmileParam,
                        bsDim,
                        marketCurves.GetDiffuse().BaseDate(),
                        mBSType,
                        marketCurves.mCV[marketCurves.GetCurveType(mLiborCVName)],
                        marketCurves.mCV[KV_BASIS],
                        mLiborDCC,
                        mBasisDCC,
                        bsVolDiag,
                        factVolDates,
                        factVolCurves);
        }


        // Collect spot vols (for testing)
        //



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

        KTMXirTree::Initialize(
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


    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
        
}

/* Reset Bs Spot Vols*/
void
KBirTree::ResetBsSpotVol(
        KMarketCurves&   marketCurves,     // (I) curves and curve types

        KVolDiag&        irVolDiag,        // (I) IR volatility data.
        KMrParam&        irMrParam,        // (I) IR mr data.
        KSmileParam&     irSmileParam,     // (I) IR skew data.
        KVolDiag&        bsVolDiag,        // (I) Basis volatility data.
        KMrParam&        bsMrParam,        // (I) Basis mr data.
        KSmileParam&     bsSmileParam,     // (I) Basis skew data.
        double           irBsCorr,         // (I) IR basis correlation.

        KResetBank       &resetBank)       // (I) rate reset bank
{
    static  char            routine[] = "KBirTree::ResetBsSpotVol";
    if ( TMXFlag == FALSE )
    {
        return;
    }
    if (!mBSOn)
    {
        return;
    }
    try
    {
        int             i;
        int             bsDim;
        TDate           VolMatDate;
        TDateInterval   dateInterval;
        double          bsVolInterp;

        int             numBSVols = bsVolDiag.Size();
        TDate          *bsVolDates = new TDate[numBSVols];
        double         *bsVolRates = new double[numBSVols];

        KVolDiag        newBsVolDiag;

        KVector(TDate)::iterator itDate;
        KVector(double)::iterator itRate;
        
        bsDim = bsMrParam.mNumFact;
        newBsVolDiag.mVolType = bsVolDiag.mVolType;

        itDate = bsVolDiag.mVolDates.begin();
        itRate = bsVolDiag.mVolRates.begin();

        for (i = 0; i < numBSVols; i++)
        {
            bsVolDates[i] = *itDate;
            bsVolRates[i] = *itRate;
            itDate++;
            itRate++;
        }

        IF_FAILED_THROW( GtoMakeDateInterval(
                            1,
                            'Q',
                            &dateInterval));

        for ( itDate = mVolDates.begin(); 
              itDate != mVolDates.end();
              itDate++)
              {
                  IF_FAILED_THROW( DrlTDateLinearInterp1d(
                                    bsVolDates,
                                    bsVolRates,
                                    numBSVols,
                                    *itDate,
                                    &bsVolInterp));
                  newBsVolDiag.mVolDates.push_back(*itDate);
                  newBsVolDiag.mVolMats.push_back(1.0/4.0);
                  newBsVolDiag.mVolFreqs.push_back(0);
                  newBsVolDiag.mVolRates.push_back(bsVolInterp);

              }

         // ReCalib spread vol
         KVector(TDate)           volDates;
         KVector(KVector(double)) factVolCurves;
         KVector(KVector(double))::iterator itFact;
         DppCalibSpreadVolDiag(
                        bsMrParam,
                        bsSmileParam,
                        bsDim,
                        marketCurves.GetDiffuse().BaseDate(),
                        mBSType,
                        marketCurves.mCV[marketCurves.GetCurveType(mLiborCVName)],
                        marketCurves.mCV[KV_BASIS],
                        mLiborDCC,
                        mBasisDCC,
                        newBsVolDiag,
                        volDates,
                        factVolCurves);

        if ( volDates.size() > mVolDates.size())
        {
            // remove base date from factVolCurves
            for (itFact = factVolCurves.begin();
                 itFact != factVolCurves.end();
                 itFact++)
                 {
                     (*itFact).erase( (*itFact).begin());
                 }
        }
        for (itFact = factVolCurves.begin();
             itFact != factVolCurves.end();
             itFact++)
             {
                mFactVol.insert(mFactVol.end(), bsDim, *itFact);
             }

        if (debugLevel > DEBUG_LEVEL_TIMELINE) 
        {
            dppLog << format("----------------%s-----------", routine) <<endl;
            for ( itDate = newBsVolDiag.mVolDates.begin(), 
                  itRate = newBsVolDiag.mVolRates.begin(); 
                  itDate != newBsVolDiag.mVolDates.end();
                  itDate++, itRate++)
            {
                dppLog << "VolDate " << GtoFormatDate(*itDate) << '\t'
                       << "BsVol " << (*itRate) << endl;
            }
            dppLog << "--------- New MR FactVols ---------" << endl;
            KVector(KVector(double))::iterator itVol1 = mFactVol.begin();
            KVector(KVector(double))::iterator itVol2 = itVol1;
            itVol2++;
            for (itDate = mVolDates.begin(), i = 0;
                 itDate != mVolDates.end();
                 itDate++, i++)
                 {
                     dppLog << "VolDate = " << GtoFormatDate(*itDate) << '\t'
                            << "spotVol = "<< (*itVol1)[i] <<  '\t' 
                            << "spreadVol = " << (*itVol2)[i] << endl;
                 }
      }
        delete [] bsVolDates;
        delete [] bsVolRates;

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}

const String&   
KBirTree::GetIRDiscCurveName(const String& curveName)
{
        int discCurveIdx;

        discCurveIdx = GetCurveIdx(curveName);

        if (discCurveIdx < KV_BASIS)
                return curveName;
        else
                return mBSDiscCVName;
}




//---------------------------------------------------------------
//
//

KVector(double)
KBirTree::IrSpotVols()
{
        return mFactVol[0];
}




//---------------------------------------------------------------
//

KVector(double)
KBirTree::BsSpotVols()
{
        if (mBSDim <= 0) {
                throw KFailure("%s: no basis dimension %d (=%d-%d).\n",
                        "KBirTree::BsSpotVols",
                        mBSDim, mNbFactor, mIRDim);
        }
        return mFactVol[mIRDim];
}




//---------------------------------------------------------------
// Check the tree.
//

void
KBirTree::CheckTreeValid()
{
static  char    routine[] = "KBirTree::CheckTreeValid";
 
        KTMXirTree::CheckTreeValid();

        if (mBSOn) {
           // Dimension consistency
           if ((mIRDim+mBSDim) > 3 || (mIRDim+mBSDim) == 0)
                throw KFailure("%s: total dimension (%d+%d) of the tree is "
                               " either greater than 3 or equal to 0.\n",
                                routine, mIRDim, mBSDim);
        
           // Libor curve and basis curve validity
           if (GetCurveIdx(mLiborCVName) == K_DEFAULT_IDX)
                throw KFailure("%s: curve %s does not exist.\n",
                                routine, mLiborCVName.c_str());

           if (GetCurveIdx(mBasisCVName) == K_DEFAULT_IDX)
                throw KFailure("%s: curve %s does not exist.\n",
                                routine, mBasisCVName.c_str());

           if (debugLevel > DEBUG_LEVEL_TIMELINE) {
                dppLog << format("%s:\n", routine);
                dppLog << format("mBSDim = %d\n",         mBSDim);
                dppLog << format("mBSQLo = %lf\n",        mBSQLo);
                dppLog << format("mBSQHi = %lf\n",        mBSQHi);
                dppLog << format("mBSFSh = %lf\n",        mBSFSh);
                dppLog << format("mBSBbq = %lf\n",        mBSBbq);
                dppLog << format("mBSVolNorm = %lf\n",    mBSVolNorm);
                dppLog << format("mBSVolLogn = %lf\n",    mBSVolLogn);
           }

        }

}




//--------------------------------------------------------------
// Allocate a slice based on curve type
// DO NOT call directly for slice allocation. 
// Use KTSlice constructor instead
 
KTSlice&
KBirTree::TSliceCreate(KTSlice &ts)
{
static  char    routine[] = "KBirTree::TSliceCreate";
 
        int     sliceDim;
        int        curveIdx;
        int        nDim;
 
 try{
        // Check that we have anything to do
        if (!ts.IsEmpty()) return(ts);
 
        nDim = mIRDim + mBSDim;

        curveIdx = ts.GetCurveIdx();

        // Allocate slice dimension based on the curve type:
        // 1. if curve is BASIS, allocate mIRDim+mBSDim dimension.
        // 2. if curve is non-BASIS, allocate mIRDim dimension.
        if (curveIdx >= KV_BASIS)        
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
KBirTree::TSliceUnaryOper(
        KTSlice &ts,
        const KTSlice &ts1,
        KOper oper)
{
static        char        routine[] = "KBirTree::TSliceUnaryOper";
        int        tpIdx = TPIdxCurrent();


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
                //
                KTSlice tmpTS1(*this,
                               "Temporary slice for "
                               "dimension resize",
                               ts.GetCurveName());

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
                KTSlice tmpTS(*this,
                              "Temporary slice for "
                              "dimension resize",
                              ts1.GetCurveName());

                TSliceExpand(tpIdx,                
                             ts,
                             tmpTS);

                // Reallocate the slice memory
                TSliceDestroy(ts);

                // Set the new slice dimension of ts.
                ts.SetCurveIdx(tmpTS.GetCurveName());
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
KBirTree::TSliceDev(KTSlice &ts, const String &discCurveName)
{
static  char    routine[] = "KBirTree::TSliceDev";
 
        int        discCurveIdx;            // discount curve index
        int        sliceDim;                // slice dimension
        int        nDim;                    // full dimension of the tree
        int        tpIdxTree,               // current time point of the slice
                   tpIdxSlice;              // current time point of the tree
        
        double     *discount = NULL;        // discount slice in full dimension

 try{

        nDim = mIRDim + mBSDim;

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
 

        // Slice dimension
        //

        sliceDim = ts.GetSliceDim();


        // Discount curve index

        // If discCurveName = mBasisCVName, use basis reference
        // discount curve mBSDiscCVName for discouting.  This  
        // makes sure that basis slice discouting in basis
        // spread and calibration calculation would be done 
        // correctly.
        //
        
        if (discCurveName == mBasisCVName)
                discCurveIdx = GetCurveIdx(mBSDiscCVName);
        else
                discCurveIdx = GetCurveIdx(discCurveName);
        
        if (discCurveIdx >= KV_BASIS)
                throw KFailure("%s: invalid discount curve %s with curve "
                               "index %d.\n"  
                               "Discount curve can only be diffuse "
                               "curve or curves with deterministic spreads "
                               "with respect to diffuse curve.\n",
                               routine,
                               discCurveName.c_str(),
                               discCurveIdx);

        //
        // Discount slice dimension is mIRDim, need to dev slice 
        // differently based on its dimension:
        // 1. pure IR slice, straightforward.
        // 2. basis slice, need to expand the discount slice in nDim.
        //

        if (sliceDim == mIRDim) {
                KTMXirTree::TSliceDev(ts, discCurveName);
        }
        else if(sliceDim == nDim) {
               
                if ( GetTmxFlag()== TRUE)
                {
                    // Tmx3
                    double *unitSlice   = NULL;
                    double *NmrInvSlice = NULL;

                    unitSlice = sliceNew(nDim);
                    NmrInvSlice = sliceNew(nDim);
                    sliceScalarOper ( unitSlice,
                                      nDim,
                                      1.,
                                      COPY);

                    if ( GetCcyToNmr() == TRUE)
                    {
                        sliceExpand(tpIdxTree,
                                    nDim,
                                    mIRDim,
                                    GetNmrInvLag(discCurveIdx),
                                    NmrInvSlice);
                        sliceUnaryOper( (double*)ts.mData,
                                        nDim,
                                        NmrInvSlice,
                                        MULT);
                    }

                    KMrNTree::sliceEv( (double*) ts.mData,
                                       nDim,
                                       unitSlice,
                                       tpIdxTree);

                    if( GetNmrToCcy() == TRUE) 
                    {
                        sliceExpand(tpIdxTree,
                                    nDim,
                                    mIRDim,
                                    GetNmrInv(discCurveIdx),
                                    NmrInvSlice);

                        sliceUnaryOper( (double*) ts.mData,
                                        nDim,
                                        NmrInvSlice,
                                        DIV);
                    }

                    sliceDelete(unitSlice);
                    sliceDelete(NmrInvSlice);

                }
                else
                {
                    // Fix3
                discount = sliceNew(nDim);

                // Expand discount slice from mIRDim to nDim. 
                sliceExpand(tpIdxTree,
                            nDim,
                            mIRDim,
                            GetDiscount(discCurveIdx),
                            discount);

                KMrNTree::sliceEv(
                            (double*) ts.mData,
                            nDim,
                            discount,
                            tpIdxTree);

                sliceDelete(discount);
                }
               
        }
        else {
                throw KFailure("%s: invalid slice dimension (%d) in the tree"
                               " where nDim=%d, mBSDim=%d.\n", 
                                routine, sliceDim, nDim, mBSDim);
        }



        // Set new TP for the slice
        //

        ts.SetTpIdx(tpIdxTree);


        return(ts);

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}




//--------------------------------------------------------------
//
void
KBirTree::Update(int tpIdx)
{
static  char    routine[] = "KBirTree::Update";
        
        double  *DiscountIR = NULL,
                *Discount   = NULL;

        TDate   currDate = TPDate(tpIdx);
 try{

        // Update the transition prob and discount factor 
        // between tpIdx and tpIdx+1.
        KTMXirTree::Update(tpIdx);

        //
        // Update mBSBank if is active
        // (Slice insertion is done in CalibrateSpreadDrift)
        //
        if (mBSSprdOn) 
                mBSBank.Update(tpIdx);

        //
        // Calibrate the spread drift if is critical date
        // and store the basis rate slice.
        //
        KVector(TDate)::iterator iterDelayDate = 
                find(mBSDelayDates.begin(), mBSDelayDates.end(), currDate);

        if (iterDelayDate != mBSDelayDates.end()) {
                CalibrateSpreadDrift(tpIdx);

        }
        
                
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
// 6. Initialize temporary slice for basis rate index. 
// 7. Sort and merge mBSDelayDates.
//
void
KBirTree::TreeSetUp()
{
static  char    routine[] = "KBirTree::TreeSetUp";
        
 try{
        if (debugLevel > DEBUG_LEVEL_DRIFT) 
        { 
                dppLog << endl;
                dppLog << format("============%s===========", routine) << endl;
                dppLog << endl;
        }

 
        KTMXirTree::TreeSetUp();

        if (mBSOn) {
           // Initialize tmp slice for basis rate index
           mBSTmpRate = new KTSlice(*this,
                                 "Temporary slice for basis rate",
                                 mBasisCVName);

           // Sort and merge mBSDelayDates in ascending order
           sort(mBSDelayDates.begin(), mBSDelayDates.end());
           KVector(TDate)::iterator iterBadDateSt
                  = unique(mBSDelayDates.begin(), mBSDelayDates.end());
           mBSDelayDates.erase(iterBadDateSt, mBSDelayDates.end());

           if (debugLevel > DEBUG_LEVEL_DRIFT) 
           { 
                dppLog << "Reset Date\tDelay Date\tMaturity Date\t"
                              << "Spread" << endl;
                for (KVector(TDate)::iterator it=mBSDelayDates.begin();
                                it!=mBSDelayDates.end(); ++it)
                {
                        dppLog << format("%10s\t%10s\t%10s\t%lf",
                                GtoFormatDate(GetBasisResetDate(*it)),
                                GtoFormatDate(*it),
                                GtoFormatDate(GetBasisAccrualEndDate(*it)),
                                GetBasisFwdSpread(*it)) << endl;
                }

           }

        } // mBSOn


    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}



//--------------------------------------------------------------
//
void
KBirTree::Calibrate()
{
static  char    routine[] = "KBirTree::Calibrate";
        
 try{

        TreeSetUp();

        CheckTreeValid();

        KTMXirTree::CalibrateDrift();

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}




//--------------------------------------------------------------
// Basis rate and basis spread are treated differently based on curve idx:
// 1. Model reset date for basis RATE is at delayed reset date. 
//    curveIdx = KV_BASIS.
// 2. Model reset date for basis SPREAD is at delayed reset date.
//    curveIdx = KV_SPREAD.
// 3. Model reset date for basis par SPREAD is at reset date.
//    curveIdx = KV_PAR_SPREAD.
//
// Stub is needed if reset date is in the past, and today falls 
// between resetDate and effAEDate.  The reset value is contained in
// the reset bank.
//
// For basis rate or spread, reset value in the bank is assumed
// to be the average of past discret reset rates.  The result is a 
// weighted average of past reset rate and future rate. 
// The modelResetDate where the future rate is computed is assumed
// to be delay adjusted by the same proportional ammount between 
// valueDate and AEDate.
// 
//
TDate
KBirTree::Insert(const KRateReset &rateReset, bool isCrit)
{
static  char    routine[] = "KBirTree::Insert(RateReset)";

        int        curveIdx;
        int        basisIdx, liborIdx, discIdx;
        int        nDim = mIRDim + mBSDim;
        TDate      resetDate,                // reset date
                   effDate,                  // rate effective date
                   AEDate,                   // rate maturity date
                   valueDate,                // curve value date
                   modelResetDate,           // 
                   delayResetDate,           // reset date for eq libor rate
                   liborEffDate,             // effective date for eq libor rate
                   liborEndDate;             // maturity date for eq libor rate
        double     delayShift;
        double     resetValue;

        KDateInterval delayShiftIntvl;
        KRateReset    *liborRateReset = NULL;
        
        double        liborFRate,
                      basisFRate,
                      bsSpread;

 try{

        // Fixed rate does NOT have a valid curve name
        // associated with it.  GetCurveIdx() would fail 
        // for fixed rate.

        if (!rateReset.Rate().IsFloating())
                curveIdx = KV_DIFF;
        else
                curveIdx = GetCurveIdx(rateReset.Rate().CurveName());
        

        delayShift = (rateReset.Rate().Maturity().Years())*mBSDelayShift;
        KDateInterval delayShiftIntvl = KDateInterval(delayShift);

        // Always insert start and end date as critical dates
        // if it is BASIS rate.
        //
        if (curveIdx < KV_BASIS){        // IR rate
                modelResetDate = KTMXirTree::Insert(rateReset, isCrit);
        }
        else {         // Basis rate or basis spread

                //
                // Check if is a valid basis rate
                //
                if(!rateReset.Rate().IsSimple())
                        throw KFailure("%s: invalid basis rate. "
                                       "Only simple rate allowed.\n",
                                        routine);

                //
                // Check for consistent DCC defined by rate
                // and defined by curve.
                //
                if (long(mBasisDCC) != long(rateReset.Rate().DayCc()))
                        throw KFailure("%s: basis rate dcc (%s) is inconsistent"
                        " with basis curve dcc (%s).\n",
                        routine,
                        GtoFormatDayCountConv(rateReset.Rate().DayCc()),
                        GtoFormatDayCountConv(mBasisDCC));


                //-------------------------------------------------
                // !!!!  WARNING   !!!!!   WARNING !!!
                // We strip the spread out, it will be added back
                // during the Get.
                //-------------------------------------------------
                KRateReset rateResetNS = rateReset;
                rateResetNS.Rate().SetSpread(0e0);

 
                //
                // Initialize the basis curve name.
                //
                if (mBasisCVName.empty())
                        mBasisCVName = rateResetNS.Rate().CurveName();

                //
                // Check dimension valid
                //
                if (nDim   < 1 || nDim   > 3 ||
                    mBSDim < 1 || mBSDim > 3)
                        throw KFailure("%s: invalid factor numbers.\n"
                           "IR dimension is %d, basis dimension is %d.\n",
                           routine, mIRDim, mBSDim);

                basisIdx = KV_BASIS;
                liborIdx = GetCurveIdx(mLiborCVName);
                discIdx  = GetCurveIdx(mBSDiscCVName);



                //     |<--------------Muturity------------->|
                //  |<--delay-->|
                //==r==e========dr==lEf=====================AE=========lEnd===>>
                //              |<------------- LIBOR ----------------->|
                //--------------c----c----------------------c-----------c-----


                //
                // Need stub for basis rate if reset date is before today
                //

                AEDate = rateResetNS.EffDate() + rateResetNS.Rate().Maturity();
                valueDate = GetValueDate(curveIdx);

                if (rateResetNS.ResetDate() >= mTodayDate)
                {
                    delayShift = 
                          (rateResetNS.Rate().Maturity().Years())*mBSDelayShift;
                    delayShiftIntvl = KDateInterval(delayShift);

                    resetDate = rateResetNS.ResetDate();
                    effDate   = rateResetNS.EffDate();
                    delayResetDate = resetDate + delayShiftIntvl;
                    liborEffDate = delayResetDate +
                                        rateResetNS.Rate().SpotOffset();
                    liborEndDate = liborEffDate + rateResetNS.Rate().Maturity();

                }
                else if(AEDate > valueDate)         // reset < today < eff end date
                {
                    // Reset in the past. Check if reset bank
                    // contains the reset rate value.
                
                    // For basis spread, the reset value in the
                    // bank is basis index rate, rather than spread.
                    //
                    if(curveIdx == KV_BASIS)
                    {
                        if(!mResetBank->Get(rateResetNS, &resetValue))
                           throw KFailure("%s: invalid basis rate reset "
                                    "in the past (reset %s < today %s).\n"
                                    "Could not find rate reset value from "
                                    "the reset bank. \n",
                                    routine,
                                    GtoFormatDate(rateResetNS.ResetDate()),
                                    GtoFormatDate(mTodayDate));
                    }
                    else if (curveIdx == KV_SPREAD) // store basis rate
                    {
                            // Basis index rate with same reset date
                               KRateReset basisReset ( mBasisCVName,
                                                          rateResetNS.ResetDate(),
                                                       rateResetNS.EffDate(),
                                                       rateResetNS.Rate());

                            if(!mResetBank->Get(basisReset, &resetValue))
                           throw KFailure("%s: invalid basis rate reset "
                                    "in the past (reset %s < today %s).\n"
                                    "Could not find rate reset value from "
                                    "the reset bank. \n",
                                    routine,
                                    GtoFormatDate(rateResetNS.ResetDate()),
                                    GtoFormatDate(mTodayDate));
                    }
                    else        // stored basis par spread in the bank 
                    {
                            if(!mResetBank->Get(rateResetNS, &resetValue))
                           throw KFailure("%s: invalid basis par spread "
                              "reset in the past (reset %s < today %s).\n"
                              "Could not find basis par spread from "
                              "the reset bank. \n",
                              routine,
                              GtoFormatDate(rateResetNS.ResetDate()),
                              GtoFormatDate(mTodayDate));
                    }


                    resetDate = mTodayDate;
                    effDate   = resetDate + rateResetNS.Rate().SpotOffset();

                    KDateInterval rtMaturity(AEDate, effDate);
                    delayShift = (rtMaturity.Years())*mBSDelayShift;
                    delayShiftIntvl = KDateInterval(delayShift);

                    delayResetDate = resetDate + delayShiftIntvl;
                    liborEffDate   = delayResetDate 
                                   + rateResetNS.Rate().SpotOffset();
                    liborEndDate   = liborEffDate + rateResetNS.Rate().Maturity();

                }


                //
                // model reset date where the rate slice is computed
                // on the rollback process.
                //
                modelResetDate = delayResetDate;

                //
                // Insert critical dates and compute fwd rate and spread
                //
                Insert(delayResetDate);            // basis calibration
                Insert(AEDate);                    // basis index discounting
                Insert(liborEffDate);            // calculation of libor index
                Insert(liborEndDate);            // calculation of libor index

                
                //
                // Insert libor rate reset on delayResetDate in the tree.
                //
                liborRateReset = new KRateReset(mLiborCVName,
                                                    delayResetDate,
                                                    liborEffDate,
                                                    rateResetNS.Rate());

                // Use the same rate tenor 
                // but reset to Libor DCC
                liborRateReset->Rate().SetDayCc(mLiborDCC);

                KTMXirTree::Insert(*liborRateReset, isCrit);
                
                //
                // Insert zero reset between AEDate and delayResetDate
                // for discounting
                //
                GetZeroBank(discIdx).InsertDates(AEDate,
                                                 delayResetDate);

                //
                // Insert delayResetDate in express DEV dates 
                //
                //mDevDates.push_back(delayResetDate);
                KTMXirTree::InsertStatePrice(delayResetDate);

                //
                // Insert all dates info
                //
                mBSDelayDates.push_back(delayResetDate);
                    mBSResetDates.insert(TDate_TDate(delayResetDate, 
                                                 resetDate));
                    mBSAEDates.insert(TDate_TDate(delayResetDate, AEDate));

                    mLiborRates.insert(
                        KMap(TDate, KRateReset*)::value_type(delayResetDate,
                                                             liborRateReset));
                //
                // Compute forward 1-period rate
                //
                liborFRate = liborRateReset->Rate().Forward(GetKZCurve(liborIdx),
                                                     liborEffDate);
                basisFRate = rateResetNS.Rate().Forward(GetKZCurve(basisIdx),
                                                     effDate);
                //
                // Store the basis forward rate and spread
                //
                    mBSTpFRates.insert(TDate_Double(delayResetDate, 
                                                basisFRate));


                //----------------------------------------------
                // Check sign of spread for ln distribution
                //----------------------------------------------

                if (mBSType == SUB_SPREAD)
                {
                        bsSpread = liborFRate - basisFRate;

                        // Check positivity of spreads for lognormal process
                        if ((IS_ALMOST_ZERO(mBSQLo-1e0) &&
                              IS_ALMOST_ZERO(mBSQHi-1e0)) 
                             && 
                             !(bsSpread > 0e0))        
                                throw KFailure("%s: basis spread (%f) between "
                                        "%s %s and %s forward rates on %s <= 0 "
                                        "for lognormal process.\n",
                                        routine, 
                                        bsSpread,
                                        DrlTDateIntervalPrint(NULL,
                                         (TDateInterval)(liborRateReset->Rate().Maturity())),
                                        GetCurveName(liborIdx).c_str(),
                                        GetCurveName(basisIdx).c_str(),
                                        GtoFormatDate(effDate));

                            mBSTpFSprds.insert(TDate_Double(delayResetDate, 
                                                        bsSpread));
                }
                else if (mBSType == ADD_SPREAD)
                {
                        bsSpread = basisFRate - liborFRate;

                        // Check positivity of spreads for lognormal process
                        if ((IS_ALMOST_ZERO(mBSQLo-1e0) &&
                              IS_ALMOST_ZERO(mBSQHi-1e0)) 
                             && 
                             !(bsSpread > 0e0))        
                                throw KFailure("%s: basis spread (%f) between "
                                        "%s %s and %s forward rates on %s <= 0 "
                                        "for lognormal process.\n",
                                        routine, 
                                        bsSpread,
                                        DrlTDateIntervalPrint(NULL,
                                         (TDateInterval)(liborRateReset->Rate().Maturity())),
                                        GetCurveName(liborIdx).c_str(),
                                        GetCurveName(basisIdx).c_str(),
                                        GtoFormatDate(effDate));

                            mBSTpFSprds.insert(TDate_Double(delayResetDate, 
                                                        bsSpread));
                }
                else        // PERCENTAGE
                {
                        if(liborFRate > 0e0 && 
                           basisFRate > 0e0)
                                bsSpread = basisFRate/liborFRate;
                        else
                                throw KFailure("%s: Cannot compute basis "
                                        "percetage spread.  Either libor "
                                        "forward rate (%f) on %s <= 0, or "
                                        "basis forward rate (%f) on %s <= 0.\n",
                                        routine,
                                        liborFRate,
                                        GtoFormatDate(liborEffDate),
                                        basisFRate,
                                        GtoFormatDate(effDate));

                        mBSTpFSprds.insert(TDate_Double(delayResetDate, 
                                                        bsSpread));
                }

                //
                // For KV_PAR_SPREAD curve, need to insert libor reset
                // on reset date,  where the "PAR" basis spread 
                // will be calculated..
                //
                // If reset is in the past, and today is before
                 // accural end date, then the modelResetDate is
                // still at delayedResetDate so to maintain the remaining
                // time value due to averaging/compounding.        
                //
                if (curveIdx == KV_PAR_SPREAD)
                {
                        mBSSprdOn = true;

                        Insert(resetDate);                // basis reset 
                        Insert(effDate);                // basis reset

                        KRateReset liborReset(mLiborCVName,
                                              resetDate,
                                              effDate,
                                              rateResetNS.Rate());

                        // Set to Libor DCC
                        liborReset.Rate().SetDayCc(mLiborDCC);

                        // Insert libor at reset date 
                        KTMXirTree::Insert(liborReset, isCrit);
                
                        //
                        // Insert zero reset between AEDate and resetDate
                        //
                        GetZeroBank(discIdx).InsertDates(AEDate,
                                                           resetDate);

                        // Insert zero reset between delayResetDate 
                        // and resetDate in the basis bank to be
                        // used to calculate basis rate discounted
                        // at reset date
                        mBSBank.InsertDates(delayResetDate,
                                            resetDate);

                        //
                        // Redefine modelResetDate to be reset date.
                        //
                        modelResetDate = resetDate;
                }
                
        }        // if curveIdx

        return modelResetDate;
                
    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}



//--------------------------------------------------------------
// For basis rate, we need to ensure the libor rate can be
// computed on any date between reset date and maturity date.
//
TDate
KBirTree::Insert(const KRateReset &rateReset, TDate endDate, bool isCrit)
{
static  char    routine[] = "KBirTree::Insert(RateReset, EndDate)";

        int     curveIdx;
        TDate   modelResetDate, lastZeroDate;


 try{

        // Fixed rate does NOT have a valid curve name
        // associated with it.  GetCurveIdx() would fail 
        // for fixed rate.

        if (!rateReset.Rate().IsFloating())
                curveIdx = KV_DIFF;
        else
                curveIdx = GetCurveIdx(rateReset.Rate().CurveName());
        
        // Always insert start and end date as critical dates
        // if it is BASIS rate.
        //
        if (curveIdx < KV_BASIS){        // IR rate
                modelResetDate = KTMXirTree::Insert(rateReset, 
                                                  endDate,
                                                  isCrit);
        }
        else {         // Basis rate or basis spread

                // insert the first basis reset rate
                modelResetDate = Insert(rateReset, isCrit);

                // Insert zeroReset between modelResetDate and
                // endDate + delayInterval + rateReset.SpotOffset + rateReset.Maturity 
                // or modelResetDate + 2*Maturity
                // in the tree and zero bank to ensure that 
                // reference libor rate calculation prior to reaching 
                // modelResetDate can be computed successfully.
                //
                lastZeroDate = modelResetDate + 2 * rateReset.Rate().SpotOffset()
                             + 2 * rateReset.Rate().Maturity();

                KZeroReset extraZeroReset(mLiborCVName,
                                              modelResetDate,
                                              lastZeroDate);

                KTMXirTree::Insert(extraZeroReset, isCrit);
                
        }        // if curveIdx

        return modelResetDate;
                
        
    }
    catch (KFailure) {
        throw KFailure("%s: failed inserting between %s and %s.\n", 
                        routine,
                        GtoFormatDate(rateReset.ResetDate()),
                        GtoFormatDate(endDate));
    }

}





//--------------------------------------------------------------
// Called AFTER tree update.  Check the curve type:
// 1. If curve is NON-basis, then use KTMXirTree:Get
// 2. If curve is KV_BASIS, then compute the rate from
//    libor index and spread based on following:
//    a) If tpIdx is a critical reset date, retrieve from 
//       the stored basis rate slice after tree update.
//    b) If tpIdx is NOT a critical reset date, use extrapolation
//       from previous spreads, and add/multiply back to libor. 
// 3. If curve is KV_SPREAD, then compute the spread from spread curve
//    using GetSpread.
// 4. If curve is KV_PAR_SPREAD, then compute the spread from mBSBank.
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
KBirTree::Get(KTSlice &rateTS, const KRateReset &rateReset)
{
static  char    routine[] = "KBirTree::Get(RateReset)";
        
        int        t = TPIdxCurrent();
        TDate      resetDate, effDate, AEDate;
        TDate      valueDate;
        TDate      currDate = TPDateCurrent();
        TDate      delayResetDate = -1; 

        double     tPast;
        double     liborResetValue, basisResetValue;

        double     *spread = NULL;
        double     *liborRateIR = NULL,
                   *liborRate   = NULL;

        double     resetFwd, currFwd;

        int        nDim = mIRDim + mBSDim;               // The full dimension

        String     rateTSCVName;

        int        curveIdx;

        double     drift;                                // drift of spread

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
        AEDate    = rateReset.EffDate() + rateReset.Rate().Maturity();
        valueDate = GetValueDate(curveIdx);

        // Save the curve name if available
        rateTSCVName = rateTS.GetCurveName();

        // Set the time point of the slice
        rateTS.SetTpIdx(t);

        // check rate type
        if(curveIdx < KV_BASIS)
        {

            KPirTree::Get(rateTS, rateReset);

        }
        else    // Basis rate or spread
        {
            //-------------------------------------------------
            // !!!!  WARNING   !!!!!   WARNING !!!
            // We strip the spread out, it will be added back
            // during the Get.
            //-------------------------------------------------
            KRateReset rateResetNS = rateReset;
            rateResetNS.Rate().SetSpread(0e0);
    
 
            //
            //
            // Get from reset bank if has already been manually
            // reset and past accural end date
            //
            if(valueDate >= AEDate &&
               mResetBank->Get(rateResetNS, &basisResetValue))
            {

                rateTS = basisResetValue;

            }
            else if (curveIdx == KV_BASIS)
            {
                // Check slice dimension for consistency
                if(rateTS.GetSliceDim() != nDim)
                    throw KFailure("%s: dimension of basis rate slice (%s) %d "
                                          "!= total nDim %d.\n", 
                                   routine, 
                                   rateTS.GetSliceName().c_str(),
                                   rateTS.GetSliceDim(),
                                   nDim);

                spread    = sliceNew (nDim);
                liborRate = sliceNew (nDim);

                KTSlice liborTS(*this,
                                format("Libor slice on %s.", 
                                GtoFormatDate(currDate)),
                                mLiborCVName);

                KVector(TDate)::iterator iterResetDate = 
                  find(mBSDelayDates.begin(), mBSDelayDates.end(), currDate);

                ////////////////////////////////////////////////////////////
                //
                // Interpolation is need if current date is not on 
                // critical modelResetDate.  
                // 1. For basis with SIMPLE/PAR stub, since the stub rate 
                //    will be simple rate reset at current date, only 0 
                //    delay shift is allowed.  Otherwise, the calibration 
                //    of spread will not be consistent with forward rate.
                // 2. For basis with BOND/NONE sub, resetDate < currDate. 
                //    We need to adjust the rate slice by the ratio of
                //    forward rates.
                //
                // Stub might be needed if current date is on today's date.
                // Since we calibrate the spread based on each underlying
                // rate(both basis and libor), as stored in mLiborRates,
                // it will ensure the stub rate be handled correctly.
                //
                ////////////////////////////////////////////////////////////
                if (iterResetDate == mBSDelayDates.end()) 
                {
                    // No exact match of date, need extrapolation of 
                    // calibrated drift on the right side 
                    // (flat extrapolation only).

                    KVector(TDate)::iterator iterResetDateRight = 
                        lower_bound(mBSDelayDates.begin(), mBSDelayDates.end(), 
                        currDate);
        
                    if (iterResetDateRight == mBSDelayDates.end()) //right most
                        throw KFailure("%s: failed to interpolate the "
                                       "spread center offset.\n"
                                       "The model reset date of %s: %s  >  "
                                       "the last date of basis delayed "
                                       "dates %s.\n",
                                        routine,
                                        rateTS.GetSliceName().c_str(),
                                        GtoFormatDate(currDate),
                                        GtoFormatDate(mBSDelayDates.back()));

                    // flat extroplation
                    drift = GetBasisSpreadCenter(*iterResetDateRight);

                    CalcSpread(t, drift, spread);
                
                    // Get the libor index rate slice on current date 
                    // Assuming libor index is the same as basis index 
                    // (same maturity and dcc, etc.)
                    effDate = currDate + rateResetNS.Rate().SpotOffset();
                    KRateReset liborDelayReset ( mLiborCVName,
                                                      currDate,
                                                      effDate,
                                                      rateResetNS.Rate());

                    // Set to Libor DCC
                    liborDelayReset.Rate().SetDayCc(mLiborDCC);

                    KPirTree::Get(liborTS, liborDelayReset);
                    liborRateIR = (double*)liborTS.mData;
                    sliceExpand(t,
                                    nDim,
                                    mIRDim,
                                    liborRateIR,
                                    liborRate);

                    sliceUnaryOper( (double*)rateTS.mData,
                                    nDim,
                                    liborRate,
                                    COPY);

                    if (mBSType == SUB_SPREAD)
                        sliceUnaryOper((double*)rateTS.mData,
                                               nDim,
                                        spread,
                                        SUB);        
                    else if (mBSType == ADD_SPREAD)
                        sliceUnaryOper((double*)rateTS.mData,
                                               nDim,
                                        spread,
                                        ADD);        
                    else
                        sliceUnaryOper((double*)rateTS.mData,
                                               nDim,
                                        spread,
                                        MULT);        

                    //
                    // rateTS is a basis rate reset at current date.
                    // If resetDate < currDate, we need to adjust it by 
                    // the ratio of forwards.
                    //
                    if (resetDate < currDate &&
                        resetDate >= mTodayDate) {
                        resetFwd = rateResetNS.Rate().Forward(GetKZCurve(KV_BASIS),
                                                     resetDate);
 
                        currFwd  = rateResetNS.Rate().Forward(GetKZCurve(KV_BASIS),
                                                     currDate);

                        rateTS *= (resetFwd / currFwd);
                    }
                    else if (resetDate < currDate &&
                             resetDate < mTodayDate) {        
                        //
                        // reset in the past, calculate the weighted
                        // average of past and future rates.
                        // The future rate start from today is again estimated
                        // from the ratio of forwards
                        //
                        resetFwd = rateResetNS.Rate().Forward(GetKZCurve(KV_BASIS),
                                                     mTodayDate);
 
                        currFwd  = rateResetNS.Rate().Forward(GetKZCurve(KV_BASIS),
                                                     currDate);

                        rateTS *= (resetFwd / currFwd);

                        if(!mResetBank->Get(rateResetNS, &basisResetValue))
                        throw KFailure("%s: invalid KRateReset for reset "
                                         "in the past (reset %s < today %s).\n"
                                         "Could not find rate reset value from "
                                         "the reset bank. \n",
                                         routine,
                                         GtoFormatDate(resetDate),
                                         GtoFormatDate(mTodayDate));

                        // The result is a weighted average of past 
                        // reset rate and future rate.
                        tPast   = (double)(mTodayDate - resetDate)
                                        /(double)(AEDate - effDate);

                        basisResetValue *= tPast;
                        
                        rateTS *= (1. - tPast);
                        rateTS += basisResetValue;
                    }
                }
                else     // falls on calibrated date  
                {
                    if (resetDate >= mTodayDate) 
                    {
                        rateTS = *mBSTmpRate;
                    }
                    else  
                    {
                        // Stub for reset in the past, and today falls
                        // between resetDate and effAEDate. 
                        // Check if reset bank contains the reset 
                        // rate value.
                            if(!mResetBank->Get(rateResetNS, &basisResetValue))
                            throw KFailure("%s: invalid KRateReset for reset "
                                         "in the past (reset %s < today %s).\n"
                                         "Could not find rate reset value from "
                                         "the reset bank. \n",
                                         routine,
                                         GtoFormatDate(resetDate),
                                         GtoFormatDate(mTodayDate));


                        // The result is a weighted average of past 
                        // reset rate and future rate.
                            tPast   = (double)(mTodayDate - resetDate)
                                        /(double)(AEDate - effDate);

                            basisResetValue *= tPast;

                            rateTS = *mBSTmpRate;
                            rateTS *= (1. - tPast);
                            rateTS += basisResetValue;

                    }

                }

                //-------------------------------------------------
                // !!!!  WARNING   !!!!!   WARNING !!!
                // We finally add the spread
                //-------------------------------------------------
                rateTS += rateReset.Rate().Spread();


                // Free used memory 
                sliceDelete(spread);
                sliceDelete(liborRate);

            }
            else if (curveIdx==KV_SPREAD)          // basis spread
            {

                // 1) If reset in the FUTURE, then the modelResetDate
                //    for spread is at delayed reset date.  The spread
                //    is computed from spread curve using GetSpread()
                //
                // 2) If reset in the PAST, and today is BEFORE the 
                //    accural end date, then the modelResetDate for
                //    spread is at delayed reset date between today
                //    and accural end date so to keep the time value
                //    due to averaging/compounding.  The spread is
                //    computed from reset libor stored in the reset bank
                //    and weighted average basis index rate from past reset
                //    and future value computed at modelResetDate.
                //
                // 3) If reset in the PAST, and today is AFTER the
                //    accural end date, then retrieve the reset value 
                //    of spread from the reset bank.


                // Check slice dimension for consistency
                if(rateTS.GetSliceDim() != nDim)
                    throw KFailure("%s: dimension of basis spread slice (%s) "
                                          "%d != total nDim %d.\n", 
                                   routine, 
                                   rateTS.GetSliceName().c_str(),
                                   rateTS.GetSliceDim(),
                                   nDim);

                // Check date consistency
                ASSERT_OR_THROW(resetDate <= currDate);
                ASSERT_OR_THROW(currDate  <= AEDate);



                //  Case 3). valueDate > AEDate.
                if (valueDate >= AEDate)
                {
                    // Reset in the past. Check if reset bank
                    // contains the reset rate value.
                    if(!mResetBank->Get(rateResetNS, &basisResetValue))
                        throw KFailure("%s: invalid KRateReset for basis "
                            "spread reset in the past (reset %s < today %s).\n"
                            "Could not find rate reset value from the "
                            "reset bank. \n",
                            routine,
                            GtoFormatDate(resetDate),
                            GtoFormatDate(mTodayDate));

                    rateTS = basisResetValue;
                }
                else if (resetDate >= mTodayDate) //case 1), reset in the future
                {

                    GetSpread(rateTS, currDate);

                }
                else   // Case 2), today falls between reset date and AEDate.
                {

                    // Libor index rate with same reset date
                    KRateReset liborReset ( mLiborCVName,
                                               resetDate,
                                               effDate,
                                               rateResetNS.Rate());

                    // Set to Libor DCC
                    liborReset.Rate().SetDayCc(mLiborDCC);

                    // Basis index rate with same reset date
                    KRateReset basisReset ( mBasisCVName,
                                                  resetDate,
                                               effDate,
                                               rateResetNS.Rate());

                    // Get the past libor reset rate from reset bank
                    //
                    if(!mResetBank->Get(liborReset, &liborResetValue))
                        throw KFailure("%s: invalid KRateReset for libor "
                            "reset in the past (reset %s < today %s).\n"
                            "Could not find libor reset value from the "
                            "reset bank. \n",
                            routine,
                            GtoFormatDate(resetDate),
                            GtoFormatDate(mTodayDate));

                    // Get the past basis reset rate from reset bank
                    //
                    if(!mResetBank->Get(basisReset, &basisResetValue))
                        throw KFailure("%s: invalid KRateReset for basis "
                            "spread reset in the past (reset %s < today %s).\n"
                            "Could not find rate reset value from the "
                            "reset bank. \n",
                            routine,
                            GtoFormatDate(resetDate),
                            GtoFormatDate(mTodayDate));


                    tPast   = (double)(mTodayDate - resetDate)
                                    /(double)(AEDate - effDate);

                    basisResetValue *= tPast;

                    // Check if current date is a valid model reset date
                    //
                    KVector(TDate)::iterator iterMRDate = 
                    find(mBSDelayDates.begin(), mBSDelayDates.end(), currDate);

                    if (iterMRDate != mBSDelayDates.end()) 
                    {
                        rateTS = *mBSTmpRate;
                        rateTS *= (1. - tPast);
                        rateTS += basisResetValue;
                    }
                    else
                    {
                        throw KFailure("%s: can not get basis rate that "
                                "reset in the past (%s) at current date (%s), "
                                "which is NOT a valid model reset date.\n",
                                routine,
                                GtoFormatDate(resetDate),
                                GtoFormatDate(currDate));
                    }


                    if (mBSType == SUB_SPREAD)
                    {
                        rateTS *= -1e0;
                        rateTS += liborResetValue;
                    }
                    else if (mBSType == ADD_SPREAD)
                    {
                        rateTS -= liborResetValue;
                    }
                    else
                        rateTS /= liborResetValue;


                }        // end of case 2)

            }   
            else if (curveIdx==KV_PAR_SPREAD)          // basis par spread
            {

                // The modelResetDate is at the libor reset date
                //
                // 1) If reset in the FUTURE, then the spread
                //    is computed from Libor and basis index rate both
                //    discounted to reset date.
                //
                // 2) If reset in the PAST, then retrieve the reset value 
                //    of spread from the reset bank.
                //

                // Check slice dimension for consistency
                if(rateTS.GetSliceDim() != nDim)
                    throw KFailure("%s: dimension of basis spread slice (%s) "
                                   "%d != total nDim %d.\n", 
                                    routine, 
                                    rateTS.GetSliceName().c_str(),
                                    rateTS.GetSliceDim(),
                                    nDim);

                // Check date consistency
                ASSERT_OR_THROW(resetDate <= currDate);
                ASSERT_OR_THROW(currDate  <= AEDate);


                // Check if it's already in the reset bank
                //
                if(mResetBank->Get(rateResetNS, &basisResetValue))
                    rateTS = basisResetValue;
                else if (mTodayDate > resetDate)  //  Case 2). reset in the past
                {
                    throw KFailure("%s: invalid KRateReset for basis par "
                          "spread reset in the past(reset %s< today %s).\n"
                          "Could not find basis par spread reset value "
                          "from the reset bank. \n",
                          routine,
                          GtoFormatDate(resetDate),
                          GtoFormatDate(mTodayDate));
                }
                else          // case 1) reset in the future        
                {
                    // Libor slice with nIRDim dimension
                    KTSlice liborIR(*this,
                                    format("Tmp slice on %s.", 
                                    GtoFormatDate(currDate)),
                                    mLiborCVName);

                    KRateReset liborReset_p( mLiborCVName,
                                             resetDate,
                                             effDate,
                                             rateResetNS.Rate());

                    // Set to Libor DCC
                    liborReset_p.Rate().SetDayCc(mLiborDCC);

                    // Get the libor rate slice
                    // If reset in the past, will get from the reset bank
                    //
                    KPirTree::Get(liborIR, liborReset_p);


                    // Slice for discount between reset and accural end
                    // dates with nIRDim dimension
                    //
                    KTSlice discIR( *this,
                                    format("Tmp slice on %s.", 
                                            GtoFormatDate(currDate)),
                                    mBSDiscCVName);

           
                    KZeroReset discZero (mBSDiscCVName,
                                         resetDate,
                                         AEDate);

                    KTMXirTree::Get(discIR,  discZero);
        

                    // Discounted libor is L/(1+L*dcf)
                    liborIR *= discIR;


                    // Find the discounted basis index rate from basis 
                    // bank.
                    // Find out the delay reset date corresponeing to
                    // reset date
                    for(KMap(TDate, TDate)::iterator 
                             iter_p = mBSBank.mErDates.begin();
                             iter_p != mBSBank.mErDates.end(); ++iter_p)
                    {
                        if (currDate == (*iter_p).second)
                            delayResetDate = (*iter_p).first;
                    }

                    if (delayResetDate != -1)          // match
                        mBSBank.GetSlice(rateTS, delayResetDate);
                    else
                        throw KFailure("%s: reset date (%s) is not found "
                                              "in the basis bank mErDates.\n", 
                                       routine, 
                                       GtoFormatDate(currDate));


                    if (mBSType == SUB_SPREAD)
                    {
                        rateTS *= -1e0;
                        rateTS += liborIR;
                    }
                    else if (mBSType == ADD_SPREAD)
                    {
                        rateTS -= liborIR;
                    }
                    else
                        rateTS /= liborIR;


                }        // end if stub

            }   // if KV_PAR_SPREAD
            else        // curveIdx
                throw KFailure("%s: invalid curve index %d.\n", 
                                routine, curveIdx);

        }  // end of if reset is true

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
        sliceDelete(spread);
        sliceDelete(liborRate);
        throw KFailure("%s: failed.\n", routine);
    }
}



//--------------------------------------------------------------
// Get the zero reset. Discount curve can only be diffuse, or
// deterministic spread curve
//
void
KBirTree::Get(KTSlice &ts, const KZeroReset &zeroReset)
{
static  char        routine[] = "KBirTree::Get(ZeroReset)";
 
        int     curveIdx;
 
        // Curve for zero reset has to be either diffuse
        // or deterministic spread curve
        curveIdx = GetCurveIdx(zeroReset.mCurveName);
        if (curveIdx >= KV_BASIS)
                throw KFailure("%s: Discount curve can NOT be basis curve.\n",
                                routine);
 
        KTMXirTree::Get(ts, zeroReset);
}




//--------------------------------------------------------------
// Get the basis spread slice on a delayed reset date.
//
void
KBirTree::GetSpread(KTSlice &spreadTS, TDate resetDate)
{
static  char    routine[] = "KBirTree::GetSpread";
        
        int     t = TPIdxCurrent();

        double  *spread = NULL;

        int     nDim = mNbFactor;                // The full dimension

        double  drift;                           // drift of spread

try {
        
        // Check slice dimension for consistency
        if(spreadTS.GetSliceDim() != nDim)
                throw KFailure("%s: basis spread slice dimension (%d) "
                               "!= nDim(%d).\n", 
                                routine, 
                                spreadTS.GetSliceDim(),
                                nDim);

        spread = sliceNew (nDim);

        KVector(TDate)::iterator iterResetDateRight = 
                lower_bound(mBSDelayDates.begin(), mBSDelayDates.end(), 
                            resetDate);

        // Flat interpolate the drift from the nearest calibrated date
        // on the right side.
        if (iterResetDateRight == mBSDelayDates.end()) 
                throw KFailure("%s: model reset date (%s) is larger than "
                               "the maximum basis benchmark date(%s).\n",
                               routine,
                               GtoFormatDate(resetDate),
                               GtoFormatDate(mBSDelayDates.back()));


        // Center offset
        drift = GetBasisSpreadCenter(*iterResetDateRight);

        CalcSpread(t, drift, spread);

        sliceUnaryOper( (double*)spreadTS.mData,
                        nDim,
                        spread,
                        COPY);
                
        // Set the time point of the slice
        spreadTS.SetTpIdx(t);

        sliceDelete(spread);

    }
    catch (KFailure) {
        sliceDelete(spread);
        throw KFailure("%s: failed.\n", routine);
    }
}



//--------------------------------------------------------------
//
 
TDate
KBirTree::GetBasisResetDate(TDate delayDate)
{
static  char    routine[] = "KBirTree::GetBasisResetDate";
 
        KMap(TDate, TDate)::iterator it = mBSResetDates.find(delayDate);
 
        if (it == mBSResetDates.end())
                throw KFailure("%s: no reset date corresponding to delay "
                               "date %s stored in the table.\n",
                               routine, GtoFormatDate(delayDate));
 
        return (*it).second;
 
}



//--------------------------------------------------------------
//
 
TDate
KBirTree::GetBasisAccrualEndDate(TDate delayDate)
{
static  char    routine[] = "KBirTree::GetBasisAccrualEndDate";
 
        KMap(TDate, TDate)::iterator it = mBSAEDates.find(delayDate);
 
        if (it == mBSAEDates.end())
                throw KFailure("%s: no accural end date corresponding to delay "
                               "date %s stored in the table.\n",
                               routine, GtoFormatDate(delayDate));
 
        return (*it).second;
 
}




//--------------------------------------------------------------
//
 
KRateReset*
KBirTree::GetLiborRate(TDate delayDate)
{
static  char    routine[] = "KBirTree::GetLiborRate";
 
        KMap(TDate, KRateReset*)::iterator it = mLiborRates.find(delayDate);
 
        if (it == mLiborRates.end())
                throw KFailure("%s: no reset date corresponding to delay "
                               "date %s stored in the table.\n",
                               routine, GtoFormatDate(delayDate));
 
        return (*it).second;
 
}



//--------------------------------------------------------------
// Retrieve basis forward spread from stored table.
// If isInterp is true, linear interpolate from two adjacent dates
//
 
double
KBirTree::GetBasisFwdSpread(TDate delayDate, bool isInterp)
{
static  char    routine[] = "KBirTree::GetBasisFwdSpread";
 
        TDate        date_L, date_R;                // dates for interpolation
        double        sprd_L, sprd_R;                // rates for interpolation
        double        spread = -100e0;

        KMap(TDate, double)::iterator it = mBSTpFSprds.find(delayDate);
 
        if (it != mBSTpFSprds.end())         // exact match
                spread = (*it).second;
        else if (isInterp)                // linear interpolate
        {
                KVector(TDate)::iterator iterResetDateRight = 
                        lower_bound(mBSDelayDates.begin(), mBSDelayDates.end(), 
                        delayDate);
                        
                // Check if beyond the last benchmark date
                //
                if (iterResetDateRight == mBSDelayDates.end()) // right most
                        throw KFailure("%s: failed to interpolate the "
                                       "basis forward spread.\n"
                                       "The delay date of %s  >  "
                                       "the last date of basis delayed "
                                       "dates %s.\n",
                                        routine,
                                        GtoFormatDate(delayDate),
                                        GtoFormatDate(mBSDelayDates.back()));
                
                date_R = *iterResetDateRight;

                // spreads on the right and left adjacent dates
                //
                KMap(TDate, double)::iterator itR = mBSTpFSprds.find(date_R);

                // Check if forward spreads exist on these dates.
                //
                if (itR == mBSTpFSprds.end())
                        throw KFailure("%s: no forward spread corresponding "
                                       "to delay date %s stored in the "
                                       "table.\n",
                                               routine, 
                                        GtoFormatDate(date_R));

                sprd_R = (*itR).second;

                // Check if at the first date on left boundary
                //
                if (iterResetDateRight == mBSDelayDates.begin())
                        spread = sprd_R;        // flat 
                else                // interpolate between two dates. Finally!
                {
                        date_L = *(iterResetDateRight-1);
                        
                        KMap(TDate, double)::iterator itL =
                                mBSTpFSprds.find(date_L);

                        // Check if forward spread on this date exists
                        //
                        if (itL == mBSTpFSprds.end())
                                throw KFailure("%s: no forward spread "
                                       "corresponding to delay date %s "
                                       "stored in the table.\n",
                                               routine, 
                                        GtoFormatDate(date_L));
        
                        sprd_L = (*itL).second;
                                
                        // Linear interpolate spread
                        spread = sprd_L + (delayDate - date_L)*
                                        (sprd_R - sprd_L)/(date_R - date_L);
                }
                
        }
        else                                // no match, no interpolation
                throw KFailure("%s: no forward spread corresponding to delay "
                               "date %s stored in the table.\n",
                               routine, GtoFormatDate(delayDate));
 
        return spread;
 
}




//--------------------------------------------------------------
//
 
double
KBirTree::GetBasisFwdRate(TDate delayDate)
{
static  char    routine[] = "KBirTree::GetBasisFwdRate";
 
        KMap(TDate, double)::iterator it = mBSTpFRates.find(delayDate);
 
        if (it == mBSTpFRates.end())
                throw KFailure("%s: no forward basis rate corresponding to "
                               "delay date %s stored in the table.\n",
                               routine, GtoFormatDate(delayDate));
 
        return (*it).second;
 
}




//--------------------------------------------------------------
//
 
double
KBirTree::GetBasisSpreadCenter(TDate delayDate)
{
static  char    routine[] = "KBirTree::GetBasisSpreadCenter";
 
        KMap(TDate, double)::iterator it = mBSTpCenters.find(delayDate);
 
        if (it == mBSTpCenters.end())
                throw KFailure("%s: no spread center offset corresponding to "
                               "delay date %s in the table.\n",
                               routine, GtoFormatDate(delayDate));
 
        return (*it).second;
 
}



//--------------------------------------------------------------
//
double
KBirTree::GetBasisSpreadVolBbq(int t)
{
static  char    routine[] = "KTMXirTree::GetVolBbq";
 
        TDate        currDate;

        double  QMid;
 
        double  fwdSpread;        // Fwd spread 
        double  FwdSpreadA;        // Fwd spread adjusted
 
        double  VolBbq;
 

        currDate = TPDate(t);
 
        // Forward spread
        // Linear interpolate from mBSTpFSprds
        //
        fwdSpread = GetBasisFwdSpread(currDate, true);

        // Avoid zero spread, which is allowed in normal case but
        // canceled out after the mapping.  Set the minimum to 0.01bp
        SHIFT_ZERO(fwdSpread);

        //
        // Check on mBSFSh non-singular, i.e. mBSFSh != -1,
        // QMid * mBSFSh != -1, and QRight * mBSFSh > -1, etc.
        // are done in modpar.cxx
        //
        QMid = (mBSQLo + mBSQHi) / 2;

        FwdSpreadA  = fwdSpread / (1. + mBSFSh);


        VolBbq  = mBSVolLogn * mBSBbq
                + mBSVolNorm * (1. - mBSBbq) / fwdSpread;
 
        VolBbq *= (1. + mBSFSh) / (1. + QMid * mBSFSh);
 
        if (debugLevel >= DEBUG_LEVEL_DRIFT) {
                dppLog << format("%s: TPIDX=%4d VolBbq=%14.10f\n",
                        routine, t, VolBbq);
        }

        return VolBbq;
}




//--------------------------------------------------------------
// Calibrate the drift of spread on a reset date
// so that the 1-period basis forward would be priced exactly.
// Store the current basis rate slice as well as drift for
// later use by Get routine.
//
void
KBirTree::CalibrateSpreadDrift(int  tpIdx)        // (I) time point index
{
static  char       routine[] = "KBirTree::CalibrateSpreadDrift";

        double     Zt;                        // calibrated spread drift

        TDate      delayResetDate = TPDate(tpIdx);

        TDate      AEDate,                        // basis Accrual end date
                   nextResetDate;

        int        t = tpIdx;                // for convenience

        int        bsDiscCurve;                // discount curve index

        int        nDim, nIRDim;
        
        double     fwdBasisRate;

        double     *discountIR = NULL,
                   *liborRateIRTS = NULL,
                   *liborRateTS   = NULL,
                   *spread   = NULL,
                   *statePr   = NULL;

        double     zeroPrice;

        double     zeroRatio;

        KRateReset *liborRateReset = NULL;

        double     QMid, bsVol, expT, St;    // Estimating the offset

try {

        if (debugLevel > DEBUG_LEVEL_DRIFT) {
                dppLog << format("================%s start================", 
                                  routine) << endl;
        }


        bsDiscCurve = GetCurveIdx(mBSDiscCVName);

        KVector(TDate)::iterator iterResetDate = 
                find(mBSDelayDates.begin(), mBSDelayDates.end(), 
                     delayResetDate);

        // Check if currDate is a reset date.
        if (iterResetDate == mBSDelayDates.end()) {
                throw KFailure("%s: Calibration date (tpidx = %ld, %s)"
                                " is NOT a model reset date.\n",
                                routine,
                                tpIdx, GtoFormatDate(delayResetDate));        
        }

        liborRateReset = GetLiborRate(delayResetDate);


        nDim = mIRDim + mBSDim;
        nIRDim = mIRDim;

        // Check dimension valid
        if (nDim   < 1 || nDim   > 3 ||
            mBSDim < 1 || mBSDim > 3)
                throw KFailure("%s: invalid factor numbers.\n"
                "Total dimension is %d, IR dimension is %d, "
                "basis dimension is %d.\n",
                routine, nDim, mIRDim, mBSDim);

        // Allocate tmp working space
        //
        statePr  = sliceNew(nDim);
        spread   = sliceNew(nDim);

        liborRateTS   = sliceNew(nDim);

        
        //
        // Calibrate 
        //

        // initial guess of drift
        if (delayResetDate == mBSDelayDates.back())
        {
                // Estimate from closed-form 2Q const
                //
                QMid = (mBSQLo + mBSQHi) / 2.;

                bsVol = mBSVolDiag.VolInterp(delayResetDate);
                expT  = ((double)delayResetDate - (double)TPToday()) / 365.;

                // Spot vol input 
                if (IS_ALMOST_ZERO(bsVol+1.0))
                {
                    // 
                    // !!! 1-Factor basis ONLY !!!
                    //
                    double beta    = SHIFT_ZERO(mBeta[nDim-1]);
                    double spotVol = mAlpha[nDim-1];
                    St = spotVol * sqrt(0.5*(1e0 - exp(-2.* beta * expT))/beta) 
                       * (1. + mBSFSh) / (1. + QMid * mBSFSh);
                }
                else
                    St  = bsVol * sqrt(expT) 
                       * (1. + mBSFSh) / (1. + QMid * mBSFSh);

                ConvexityC_BS2Q(St, mBSQLo, mBSQHi, mBSFSh, &Zt);
                
                Zt /= GetBasisSpreadVolBbq(t);
        }
        else {
                // Use the previous one as initial estimate
                //
                nextResetDate = *(iterResetDate+1);
                Zt = GetBasisSpreadCenter(nextResetDate);
        }


        KTSlice liborTS(*this,
                        format("Libor slice on %s.", 
                                GtoFormatDate(delayResetDate)),
                        mLiborCVName);
        KTSlice discountTS(*this,
                           format("Discount slice on %s.",
                                       GtoFormatDate(delayResetDate)),
                           mBSDiscCVName);
        KTSlice statePrTS(*this,
                          format("State price slice on %s.",
                                 GtoFormatDate(delayResetDate)),
                          mBasisCVName);
        KTSlice spreadTS(*this,
                         format("Basis spread slice on %s.",
                                 GtoFormatDate(delayResetDate)),
                         mBasisCVName);

        // Calibrate basis spreads on the roll back loop on 
        // delayResetDate.
        //
                
        // Discount zero at basis accural end date
        AEDate = GetBasisAccrualEndDate(delayResetDate);
        
        zeroPrice = GetZeroPrice(bsDiscCurve)[TPIdx(AEDate)];

        //
        // Use the previous drift as initial guess.
        //
        // CalcSpread(t, Zt, spread);

        // Get the state price at t,
        sliceUnaryOper(statePr,
                       nDim,
                       GetStatePrice(delayResetDate),
                       COPY);

        // Convert to state price of discount curve
        if (bsDiscCurve != KV_DIFF) {
                zeroRatio = GetZeroPrice(bsDiscCurve)[t]
                            /GetZeroPrice(KV_DIFF)[t];

                sliceScalarOper(statePr,
                                nDim,
                                zeroRatio,
                                MULT);
        }
                

        // Compute the underlying libor rate
        //
        KPirTree::Get(liborTS, *liborRateReset);

        liborRateIRTS = (double*)liborTS.mData;
                                
        if (debugLevel > DEBUG_LEVEL_DRIFT) {
                dppLog << endl;
                dppLog << liborTS << endl;
        }

        // Compute the discount factor between t and AEDate 
        // from zero bank.

        AEDate = GetBasisAccrualEndDate(delayResetDate);
        
        KZeroReset discZero(mBSDiscCVName,
                            delayResetDate,
                            AEDate);

        KTMXirTree::Get(discountTS, discZero);
        discountIR = (double*)discountTS.mData;

        fwdBasisRate = GetBasisFwdRate(delayResetDate);

        //
        // Use drift from previous date as initial guess
        // and solve for the new drift
        //
        SolveBasisOffsetNew(
                        t,
                        fwdBasisRate,
                        liborRateIRTS,
                        statePr,
                        discountIR,
                        zeroPrice,
                        &Zt);


        this->mBSTpCenters.insert(TDate_Double(delayResetDate, Zt));

        // Store the basis rate slice
        // Compute the spread using the new drift offset.
        CalcSpread(t, Zt, spread);

        sliceExpand(t, nDim, nIRDim, liborRateIRTS, liborRateTS);

        sliceUnaryOper((double*)(mBSTmpRate->mData),
                        nDim,
                        liborRateTS,
                        COPY);

        if (mBSType == SUB_SPREAD)
                sliceUnaryOper( (double*)mBSTmpRate->mData,
                                nDim,
                                spread,
                                SUB);        
        else if (mBSType == ADD_SPREAD)
                sliceUnaryOper( (double*)mBSTmpRate->mData,
                                nDim,
                                spread,
                                ADD);        
        else
                sliceUnaryOper( (double*)mBSTmpRate->mData,
                                nDim,
                                spread,
                                MULT);        
        
        // Set the time point for mBSTmpRate
        //
        mBSTmpRate->SetTpIdx(t);

        //
        // Add the discounted basis rate slice in mBSBank
        // to be used in spread calculation
        if (mBSSprdOn)
        {
                // Check if current date is a EvDate in the basis bank
                //
                KVector(TDate)::iterator iterEvDate =
                        find(mBSBank.mEvDates.begin(), mBSBank.mEvDates.end(), 
                             delayResetDate);

                if (iterEvDate != mBSBank.mEvDates.end())
                {
                        KTSlice *discBSSlice 
                                = new KTSlice(*this,
                                      format("Discounted basis rate "
                                      "slice delayed reset on %s",
                                      GtoFormatDate(delayResetDate)),
                                      mBasisCVName);
                        TSliceExpand(t,
                                          discountTS,
                                          *discBSSlice);

                        TSliceUnaryOper(*discBSSlice,
                                        *mBSTmpRate,
                                        MULT);

                        mBSBank.InsertSlice(delayResetDate, discBSSlice);
                    }
        }



        if (debugLevel > DEBUG_LEVEL_DRIFT) {
                dppLog << endl;
                dppLog << format("Spread at TPDate = %s:",
                          GtoFormatDate(delayResetDate)) << endl;
                slicePrint(
                           spread,
                           nDim,
                           t,
                           FALSE,
                           dppLog);

                dppLog << format("Libor at TPDate = %s:",
                          GtoFormatDate(delayResetDate)) << endl;
                slicePrint(
                           liborRateTS,
                           nDim,
                           t,
                           FALSE,
                           dppLog);
                dppLog << endl;
        }
                


        // Free tmp memory
        sliceDelete(statePr);
        sliceDelete(spread);

        sliceDelete(liborRateTS);

    }
    catch (KFailure) {
        // Free tmp memory
        sliceDelete(statePr);
        sliceDelete(spread);

        sliceDelete(liborRateTS);

        throw KFailure("%s: failed.\n", routine);
    }
}





//--------------------------------------------------------------
// Calculate the basis spread at t using the 2-Q mapping function.
//

void
KBirTree::CalcSpread(
        int    tpIdx,                // (I) time point index
        double Zt,              // (I) forward shift
        double *spread)         // (O) basis spreads
{
static        char        routine[] = "KBirTree::CalcSpread";

        double        du;                        // grid step

        double  pJump00,
                pJump10, pJump11,
                pJump20, pJump21, pJump22;

        int     t = tpIdx;               // more convenient ...
        
        TDate   currDate = TPDate(t);    // Current date

        int     nIRDim = this->mIRDim;   // IR dimensions
        int     nBSDim = this->mBSDim;   // IR dimensions
        int     nDim = mIRDim + mBSDim;  // tree dimensions

        double  QSh;                     // Used to calculate initial Zt
        double  baseSh;                  // Base shift with distribution
                                         // at X=0 corresponds to forward spread

        double  xSwitch;                 // Actual 2q switch point in X space

        double  VolBbq;                  // Sigma used in bone mapping
        double  fwdSpread,               // forward spread 
                FwdSpreadA,              // rescaled fwd spread / Fa
                //InsSpread,             // Instantaneous spread
                rJumpHi,                 // rate jump size high and low
                rJumpLo, 
                rJumpHi2,                // rate jump size high and low
                rJumpLo2, 
                rJumpHi5,                // rate jump size high and low
                rJumpLo5, 

                MLo, MHi,                // constants used in grid calc
                SLo, SHi,
                Zidx,                    // xSwitch index j0, j1, j2 adjusted
                Grid,                    // grid step
                *spreadL;                //


        int     nFact = mNbFactor,       // Used in some macros.
                i, k, j0, j1, j2,        // factor indices
                Mid;                     // switch bet lo and hi Q


        // Previous period jump size.
        //
        du = sqrt(JUMPCOEFF * LengthJ[t-1]);

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

        // Forward spread
        // Linear interpolate from mBSTpFSprds
        //
        fwdSpread = GetBasisFwdSpread(currDate, true);

        // Avoid zero spread, which is allowed in normal case but
        // canceled out after the mapping.  Set the minimum to 0.01bp
        SHIFT_ZERO(fwdSpread);


        FwdSpreadA  = fwdSpread / (1. + mBSFSh);


        // Vol backbone factor
        //
        VolBbq = GetBasisSpreadVolBbq(t);
 

        // If nBSDim = 0, deterministic spread, return FwdSpreadA.
        if (nBSDim == 0){
                sliceScalarOper(spread,
                                nDim,
                                FwdSpreadA,
                                COPY);
                return;
        }


        // The center shift Zt is relative to a base shift in 
        // the spread distribution, which is defined in such a way
        // that spread distribution at X=0 ALWAYS corresponds to the
        // forward spread. 
        //
        QSh = (mBSFSh > 0) ? mBSQHi : mBSQLo;
        if (fabs(QSh) > QCUTOFF) 
        {
                baseSh = log(1. + QSh * mBSFSh) / (QSh * VolBbq);
        }
        else
        {
                baseSh = mBSFSh / VolBbq;
        }

        // Actual 2q switch point in the X space
        //
        xSwitch = Zt + baseSh;

        
        switch (nDim) {
        case 1:
            // Set up left q mapping parameters
            if (fabs(mBSQLo) > QCUTOFF) 
            {
                MLo     = FwdSpreadA / mBSQLo;
                SLo     = FwdSpreadA - FwdSpreadA / mBSQLo;
                rJumpLo = exp(mBSQLo * VolBbq * pJump00);
            }
            else
            {
                MLo     = FwdSpreadA;
                SLo     = FwdSpreadA;
                rJumpLo = FwdSpreadA * VolBbq * pJump00;
            }
            
            // Set up right q mapping parameters
            if (fabs(mBSQHi) > QCUTOFF) 
            {
                MHi     = FwdSpreadA / mBSQHi;
                SHi     = FwdSpreadA - FwdSpreadA / mBSQHi;
                rJumpHi = exp(mBSQHi * VolBbq * pJump00);
            }
            else
            {
                MHi     = FwdSpreadA;
                SHi     = 1 + FwdSpreadA;
                rJumpHi = FwdSpreadA * VolBbq * pJump00;
            }
            
        
            switch (nBSDim) {
            case 1:          // nIRDim = 0, nBSDim = 1;

                spreadL = spread + NodeOffset(1, 0, 0, t);

                // Find switch point index

                 Zidx  = xSwitch; 
                 Mid   = (int) ceil(-Zidx / pJump00) - 1;
                 Mid   = Min(Max(Mid, mBottom1[t] - 1), mTop1[t]);


                // LEFT part of distribution
                 Zidx += (pJump00) * mBottom1[t];


                if (fabs(mBSQLo) > QCUTOFF) 
                {
                        Grid = MLo * exp (mBSQLo * VolBbq * Zidx);

                        for (i = mBottom1[t]; i <= Mid; i++) 
                        {
                                spreadL[i] = SLo + Grid;
                                Grid *= rJumpLo;
                        }           
                } 
                else 
                {
                        Grid = MLo * VolBbq * Zidx;

                        for (i = mBottom1[t]; i <= Mid; i++) 
                        {
                                spreadL[i] = SLo + Grid;
                                Grid += rJumpLo;
                        }
                }

                // RIGHT part of distribution
                Zidx = xSwitch + (pJump00) * (Mid + 1);
                
                if (fabs(mBSQHi) > QCUTOFF) {
                        Grid = MHi * exp (mBSQHi * VolBbq * Zidx);

                        for (i = Mid + 1; i <= mTop1[t]; i++) 
                        {
                                spreadL[i] = SHi + Grid;
                                Grid *= rJumpHi;
                        }
                } 
                else 
                {
                        Grid = MHi * VolBbq * Zidx;

                        for (i = Mid +1; i <= mTop1[t]; i++) 
                        {
                                spreadL[i] = SHi + Grid;
                                Grid += rJumpHi;
                        }
                }

                break;

            default:
                throw KFailure("%s: invalid factor dimensions."
                               " nDim=%d, nBSDim=%d.\n", routine, nDim, nBSDim);
            }

            break;

        case 2: 

            // Set up left q mapping parameters
            if (fabs(mBSQLo) > QCUTOFF) 
            {
                MLo      = FwdSpreadA / mBSQLo;
                SLo      = FwdSpreadA - FwdSpreadA / mBSQLo;      
                rJumpLo2 = exp (mBSQLo * VolBbq * pJump11);
            }
            else
            {
                MLo      = FwdSpreadA;
                SLo      = FwdSpreadA;      
                rJumpLo2 = FwdSpreadA * VolBbq * pJump11;
            }

            // Set-up right q mapping parameters
            if (fabs(mBSQHi) > QCUTOFF)
            {
                MHi      = FwdSpreadA / mBSQHi;
                SHi      = FwdSpreadA - FwdSpreadA / mBSQHi;      
                rJumpHi2 = exp (mBSQHi * VolBbq * pJump11);
            }
            else
            {
                MHi      = FwdSpreadA;
                SHi      = FwdSpreadA;      
                rJumpHi2 = FwdSpreadA * VolBbq * pJump11;
            }


            switch (nBSDim) {
            case 1:          // nIRDim = 1, nBSDim = 1;
                for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                {
                        // Find switch point index for spread curve
                        Zidx  = xSwitch + j0 * pJump10;
                        Mid   = (int) ceil(-Zidx / pJump11) - 1;
                        Mid   = Min(Max(Mid, mBottom2[t][j0] - 1),
                                        mTop2[t][j0]);
                        Zidx += (pJump11) * mBottom2[t][j0];

                                spreadL = spread + NodeOffset(2, j0, 0, t);
                        
                        // LEFT part of distribution
                        if (fabs(mBSQLo) > QCUTOFF) 
                        {
                                Grid = MLo * exp (mBSQLo * VolBbq * Zidx);

                                for (i = mBottom2[t][j0]; i <= Mid; i++) 
                                {
                                        spreadL[i] = SLo + Grid;
                                        Grid *= rJumpLo2;
                                }           
                        } 
                        else 
                        {
                                Grid = MLo * VolBbq * Zidx;

                                for (i = mBottom2[t][j0]; i <= Mid; i++) 
                                {
                                        spreadL[i] = SLo + Grid;
                                        Grid += rJumpLo2;
                                }
                        }
 
                        // RIGHT part of distribution
                        Zidx  = xSwitch + j0 * pJump10
                            + (Mid + 1) * pJump11;
 
                        if (fabs(mBSQHi) > QCUTOFF) 
                        {
                                Grid = MHi * exp (mBSQHi * VolBbq * Zidx);

                                for (i = Mid + 1; i <= mTop2[t][j0]; i++) 
                                {
                                        spreadL[i] = SHi + Grid;
                                        Grid *= rJumpHi2;
                                }
                        } 
                        else 
                        {
                                Grid = MHi * VolBbq * Zidx;

                                for (i = Mid +1; i <= mTop2[t][j0]; i++) 
                                {
                                        spreadL[i] = SHi + Grid;
                                        Grid += rJumpHi2;
                                }
                        }
                }   // for j0

                break;

            case 2:        // nIRDim = 0, nBSDim = 2;
                for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                {
                        // Find switch point index for spread curve
                        Zidx  = xSwitch 
                                + j0 * (pJump00 + pJump10);
                        Mid   = (int) ceil(-Zidx / pJump11) - 1;
                        Mid   = Min(Max(Mid, mBottom2[t][j0] - 1),
                                        mTop2[t][j0]);
                        Zidx += (pJump11) * mBottom2[t][j0];

                        spreadL = spread + NodeOffset(2, j0, 0, t);
                        
                        // LEFT part of distribution
                        if (fabs(mBSQLo) > QCUTOFF) 
                        {
                                Grid = MLo * exp (mBSQLo * VolBbq * Zidx);

                                for (i = mBottom2[t][j0]; i <= Mid; i++) 
                                {
                                        spreadL[i] = SLo + Grid;
                                        Grid *= rJumpLo2;
                                }           
                        } 
                        else 
                        {
                                Grid = MLo * VolBbq * Zidx;

                                for (i = mBottom2[t][j0]; i <= Mid; i++) 
                                {
                                        spreadL[i] = SLo + Grid;
                                        Grid += rJumpLo2;
                                }
                        }
 
                        // RIGHT part of distribution
                        Zidx  = xSwitch 
                            + j0 * (pJump00 + pJump10)
                            + (Mid + 1) * pJump11;
 
                        if (fabs(mBSQHi) > QCUTOFF) 
                        {
                                Grid = MHi * exp (mBSQHi * VolBbq * Zidx);

                                for (i = Mid + 1; i <= mTop2[t][j0]; i++) 
                                {
                                        spreadL[i] = SHi + Grid;
                                        Grid *= rJumpHi2;
                                }
                        } 
                        else 
                        {
                                Grid = MHi * VolBbq * Zidx;

                                for (i = Mid +1; i <= mTop2[t][j0]; i++) 
                                {
                                        spreadL[i] = SHi + Grid;
                                        Grid += rJumpHi2;
                                }
                        }
                }   // for j0

                break;

            default:
                throw KFailure("%s: invalid factor dimensions."
                               " nDim=%d, nBSDim=%d.\n", routine, nDim, nBSDim);
            }

            break;

        case 3:        // nDim = 3

            // Set up left q mapping parameters
            if (fabs(mBSQLo) > QCUTOFF) 
            {
                MLo      = FwdSpreadA / mBSQLo;
                SLo      = FwdSpreadA - FwdSpreadA / mBSQLo;      
                rJumpLo5 = exp (mBSQLo * VolBbq * pJump22);
            }
            else
            {
                MLo      = FwdSpreadA;
                SLo      = FwdSpreadA;      
                rJumpLo5 = FwdSpreadA * VolBbq * pJump22;
            }

            // Set-up right q mapping parameters
            if (fabs(mBSQHi) > QCUTOFF)
            {
                MHi      = FwdSpreadA / mBSQHi;
                SHi      = FwdSpreadA - FwdSpreadA / mBSQHi;      
                rJumpHi5 = exp (mBSQHi * VolBbq * pJump22);
            }
            else
            {
                MHi      = FwdSpreadA;
                SHi      = FwdSpreadA;      
                rJumpHi5 = FwdSpreadA * VolBbq * pJump22;
            }


            // Compute the spreads at t
            //
            switch (nBSDim) {
            case 1:          // nIRDim = 2, nBSDim = 1;
                for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
                {
                        // LEFT part of distribution
                        Zidx  = xSwitch 
                            + j0 * (pJump20)
                            + j1 * (pJump21);
                        Mid   = (int) ceil(-Zidx / pJump22) - 1;
                        Mid   = Min(Max(Mid, mBottom3[t][j0][j1] - 1),
                                        mTop3[t][j0][j1]);
                        Zidx += (pJump22) * mBottom3[t][j0][j1];
            
                        spreadL = spread + NodeOffset(3, j0, j1, t);


                        if (fabs(mBSQLo) > QCUTOFF) 
                        {
                            Grid = MLo * exp(mBSQLo * VolBbq * Zidx);

                            for (j2 = mBottom3[t][j0][j1]; j2 <= Mid; j2++) 
                            {
                                spreadL[j2] = SLo + Grid;
                                Grid *= rJumpLo5;
                            }
                        } 
                        else 
                        {
                            Grid = MLo * VolBbq * Zidx;

                            for (j2 = mBottom3[t][j0][j1]; j2 <= Mid; j2++) 
                            {
                                spreadL[j2] = SLo + Grid;
                                Grid += rJumpLo5;
                            }
                        }

                        // RIGHT part of distribution 
                        Zidx = xSwitch 
                            + j0 * (pJump20)
                            + j1 * (pJump21)
                            + (Mid+1) * (pJump22);


                        if (fabs(mBSQHi) > QCUTOFF) 
                        {
                            Grid = MHi * exp (mBSQHi * VolBbq * Zidx);

                            for (j2 = Mid + 1; j2 <= mTop3[t][j0][j1]; j2++) 
                            {
                                spreadL[j2] = SHi + Grid;
                                Grid *= rJumpHi5;
                            }
                        } 
                        else 
                        {
                            Grid = MHi * VolBbq * Zidx;

                            for (j2 = Mid + 1; j2 <= mTop3[t][j0][j1]; j2++)                                     {
                                spreadL[j2] = SHi + Grid;
                                Grid += rJumpHi5;
                            }
                        }
                }  // for j0, j1

                break;

            case 2:          // nIRDim = 1, nBSDim = 2;
                for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
                {
                        // LEFT part of distribution
                        Zidx  = xSwitch 
                            + j0 * (pJump10 + pJump20)
                            + j1 * (pJump11 + pJump21);
                        Mid   = (int) ceil(-Zidx / pJump22) - 1;
                        Mid   = Min(Max(Mid, mBottom3[t][j0][j1] - 1),
                                        mTop3[t][j0][j1]);
                        Zidx += (pJump22) * mBottom3[t][j0][j1];
            
                        spreadL = spread + NodeOffset(3, j0, j1, t);


                        if (fabs(mBSQLo) > QCUTOFF) 
                        {
                            Grid = MLo * exp(mBSQLo * VolBbq * Zidx);

                            for (j2 = mBottom3[t][j0][j1]; j2 <= Mid; j2++) 
                            {
                                spreadL[j2] = SLo + Grid;
                                Grid *= rJumpLo5;
                            }
                        } 
                        else 
                        {
                            Grid = MLo * VolBbq * Zidx;

                            for (j2 = mBottom3[t][j0][j1]; j2 <= Mid; j2++) 
                            {
                                spreadL[j2] = SLo + Grid;
                                Grid += rJumpLo5;
                            }
                        }

                        // RIGHT part of distribution 
                        Zidx = xSwitch 
                            + j0 * (pJump10 + pJump20)
                            + j1 * (pJump11 + pJump21)
                            + (Mid+1) * (pJump22);


                        if (fabs(mBSQHi) > QCUTOFF) 
                        {
                            Grid = MHi * exp (mBSQHi * VolBbq * Zidx);
                            for (j2 = Mid + 1; j2 <= mTop3[t][j0][j1]; k++) 
                            {
                                spreadL[j2] = SHi + Grid;
                                Grid *= rJumpHi5;
                            }
                        } 
                        else 
                        {
                            Grid = MHi * VolBbq * Zidx;

                            for (j2 = Mid + 1; j2 <= mTop3[t][j0][j1]; k++) 
                            {
                                spreadL[j2] = SHi + Grid;
                                Grid += rJumpHi5;
                            }
                        }
                }  // for j0, j1

                break;

            case 3:          // nIRDim = 0, nBSDim = 3;
                for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
                {
                        // LEFT part of distribution
                        Zidx  = xSwitch 
                            + j0 * (pJump00 + pJump10 + pJump20)
                            + j1 * (pJump11 + pJump21);
                        Mid   = (int) ceil(-Zidx / pJump22) - 1;
                        Mid   = Min(Max(Mid, mBottom3[t][j0][j1] - 1),
                                        mTop3[t][j0][j1]);
                        Zidx += (pJump22) * mBottom3[t][j0][j1];
            
                        spreadL = spread + NodeOffset(3, j0, j1, t);


                        if (fabs(mBSQLo) > QCUTOFF) 
                        {
                            Grid = MLo * exp(mBSQLo * VolBbq * Zidx);

                            for (j2 = mBottom3[t][j0][j1]; j2 <= Mid; j2++) 
                            {
                                spreadL[j2] = SLo + Grid;
                                Grid *= rJumpLo5;
                            }
                        } 
                        else 
                        {
                            Grid = MLo * VolBbq * Zidx;

                            for (j2 = mBottom3[t][j0][j1]; j2 <= Mid; j2++) 
                            {
                                spreadL[j2] = SLo + Grid;
                                Grid += rJumpLo5;
                            }
                        }

                        // RIGHT part of distribution 
                        Zidx = xSwitch 
                            + j0 * (pJump00 + pJump10 + pJump20)
                            + j1 * (pJump11 + pJump21)
                            + (Mid+1) * (pJump22);


                        if (fabs(mBSQHi) > QCUTOFF) 
                        {
                            Grid = MHi * exp (mBSQHi * VolBbq * Zidx);

                            for (j2 = Mid + 1; j2 <= mTop3[t][j0][j1]; j2++) 
                            {
                                spreadL[j2] = SHi + Grid;
                                Grid *= rJumpHi5;
                            }
                        } 
                        else 
                        {
                            Grid = MHi * VolBbq * Zidx;

                            for (j2 = Mid + 1; j2 <= mTop3[t][j0][j1]; j2++) 
                            {
                                spreadL[j2] = SHi + Grid;
                                Grid += rJumpHi5;
                            }
                        }
                }  // for j0, j1

                break;

            default:
                throw KFailure("%s: invalid factor dimensions."
                               " nDim=%d, nBSDim=%d.\n", routine, nDim, mBSDim);
            }

            break;

        default:
            throw KFailure("%s: N/A for nDim=%d.\n", routine, nDim);
        }

}




//--------------------------------------------------------------
// Calibrate the basis offset at idxDev so that 1-period basis rate
// discounted at the IR curve equals the 1-period forward basis rate.
// mBSType = SPREAD.

void
KBirTree::SolveBasisOffset(
        int    tpIdx,           // (I) time point index
        double fwdBasisRate,    // (I) basis forward rate
        double *spread,         // (I) spread on the time slice
        double *liborRateIR,    // (I) libor rates on the time slice
        double *statePr,        // (I) state prices on the time slice
        double *discountIR,     // (I) discount between t and T
        double zeroPrice,       // (I) discount zero at T
        double Zt,              // (I) base shift in Taylor expansion
        double *Eta)            // (O) incremental offset
{
static        char        routine[] = "KBirTree::SolveBasisOffset";

        int        t = tpIdx;                // more convenient ...

        int        nIRDim = this->mIRDim;        // IR dimensions
        int        nDim = mIRDim + mBSDim;        // tree dimensions

        if (debugLevel > DEBUG_LEVEL_DRIFT) {
                dppLog << endl;
 
                dppLog << format("=========== %s: ============", routine);
                dppLog << endl;
 
                dppLog << "Input Spread:" << endl;
 
                slicePrint(
                        spread,
                        nDim,
                        t,
                        FALSE,
                        dppLog);
 
                dppLog << endl;

                dppLog << "Input Libor Rate:" << endl;

                slicePrint(
                        liborRateIR,
                        nIRDim,
                        t,
                        FALSE,
                        dppLog);
 
                dppLog << endl;

                dppLog << "Input State Price:" << endl;
 
                slicePrint(
                        statePr,
                        nDim,
                        t,
                        FALSE,
                        dppLog);
 
                dppLog << endl;
 
                dppLog << "Input Discount:" << endl;
 
                slicePrint(
                        discountIR,
                        nIRDim,
                        t,
                        FALSE,
                        dppLog);

                dppLog << endl;
 
                dppLog << "Basis Forward Rate:\t" << setprecision(12)
                       <<  fwdBasisRate << endl;
 
                dppLog << "Zero Price:\t" << setprecision(12)
                       <<  zeroPrice << endl;
        }


        switch (mBSType){
        case SUB_SPREAD:
        case ADD_SPREAD:
                SolveBasisOffsetSpread( tpIdx,
                                          fwdBasisRate,
                                         spread,        
                                         liborRateIR,  
                                         statePr,     
                                        discountIR, 
                                        zeroPrice, 
                                        Zt,
                                        Eta);     
                break;
                
        case PER_SPREAD:
                SolveBasisOffsetPercent(tpIdx,
                                          fwdBasisRate,
                                         spread,        
                                         liborRateIR,  
                                         statePr,     
                                        discountIR, 
                                        zeroPrice, 
                                        Zt,
                                        Eta);     
                break;
         
        default:
                throw KFailure("%s: N/A for mBSType=%d.\n", routine, mBSType);
        } 

}



//--------------------------------------------------------------
// Calibrate the basis offset at idxDev so that 1-period basis rate
// discounted at the IR curve equals the 1-period forward basis rate.
// mBSType = SPREAD.

void
KBirTree::SolveBasisOffsetNew(
        int    tpIdx,           // (I) time point index
        double fwdBasisRate,    // (I) basis forward rate
        double *liborRateIR,    // (I) libor rates on the time slice
        double *statePr,        // (I) state prices on the time slice
        double *discountIR,     // (I) discount between t and T
        double zeroPrice,       // (I) discount zero at T
        double *Zt)             // (I/O) spread shift offset
{
static        char        routine[] = "KBirTree::SolveBasisOffsetNew";

        int        t = tpIdx;                // more convenient ...

        int        nIRDim = this->mIRDim;        // IR dimensions
        int        nDim = mIRDim + mBSDim;        // tree dimensions

        if (debugLevel > DEBUG_LEVEL_DRIFT) {
                dppLog << endl;
 
                dppLog << format("=========== %s: ============", routine);
                dppLog << endl;
 

                dppLog << "Input Libor Rate:" << endl;

                slicePrint(
                        liborRateIR,
                        nIRDim,
                        t,
                        FALSE,
                        dppLog);
 
                dppLog << endl;

                dppLog << "Input State Price:" << endl;
 
                slicePrint(
                        statePr,
                        nDim,
                        t,
                        FALSE,
                        dppLog);
 
                dppLog << endl;
 
                dppLog << "Input Discount:" << endl;
 
                slicePrint(
                        discountIR,
                        nIRDim,
                        t,
                        FALSE,
                        dppLog);

                dppLog << endl;
 
                dppLog << "Basis Forward Rate:\t" << setprecision(12)
                       <<  fwdBasisRate << endl;
 
                dppLog << "Zero Price:\t" << setprecision(12)
                       <<  zeroPrice << endl;
        }


        switch (mBSType){
        case SUB_SPREAD:
        case ADD_SPREAD:
                SolveBasisOffsetSpreadNew( 
                                         tpIdx,
                                         fwdBasisRate,
                                         liborRateIR,  
                                         statePr,     
                                         discountIR, 
                                         zeroPrice, 
                                         Zt);     
                break;
                
        case PER_SPREAD:
                SolveBasisOffsetPercentNew(
                                         tpIdx,
                                         fwdBasisRate,
                                         liborRateIR,  
                                         statePr,     
                                         discountIR, 
                                         zeroPrice, 
                                         Zt);     
                break;
         
        default:
                throw KFailure("%s: N/A for mBSType=%d.\n", routine, mBSType);
        } 

}



//--------------------------------------------------------------
// Calibrate the basis offset at idxDev so that 1-period basis rate
// discounted at the IR curve equals the 1-period forward basis rate.
// mBSType = SUB_SPREAD or ADD_SPREAD.
// Solving drift shift Zt using Taylor expansion and NR.
//

void
KBirTree::SolveBasisOffsetSpread(
        int    tpIdx,           // (I) time point index
        double fwdBasisRate,    // (I) basis forward rate
        double *spread,         // (I) spread on the time slice
        double *liborRateIR,    // (I) libor rates on the time slice
        double *statePr,        // (I) state prices on the time slice
        double *discountIR,     // (I) discount between t and T
        double zeroPrice,       // (I) discount zero at T
        double Zt,              // (I) base shift in Taylor expansion
        double *Eta)            // (O) incremental offset
{
static        char        routine[] = "KBirTree::SolveBasisOffsetSpread";

        double  du;                     // grid step

        double  pJump00,
                pJump10, pJump11,
                pJump20, pJump21, pJump22;

        int     t = tpIdx;              // more convenient ...

        TDate   currDate = TPDate(t);   // Current date
        
        int     nIRDim = this->mIRDim;  // IR dimensions
        int     nDim = mIRDim + mBSDim; // tree dimensions

        double  QSh;                    // Used to calculate initial Zt
        double  baseSh;                 // Base shift with distribution
                                        // at X=0 corresponds to forward spread

        double  xSwitch;                // Actual 2q switch point in X space

        double  VolBbq;                 // Sigma used in bone mapping

        double  fwdSpread,              // fwd spread
                FwdSpreadA,             // rescaled fwd spread / Fa
                //InsSpread,            // Instantaneous spread                
                Zidx,                   // Zt index j0, j1, j2 adjusted
                eIdxRate,               // expected value of index rate 
                P[NPOLY+1],             // expansion
                CQLo[NPOLY+1],          // Low coefficients
                CQHi[NPOLY+1],          // Hight coefficients
                alpha,                  // Coeff of Eta-polynomial
                x,
                uLo,
                uHi;


        double  *statePrIR = NULL,      // variables in nIRDim dimensions
                *statePrIRL,
                *liborRateIRL,
                *discountIRL;

        double  *spreadL,               // variables in nDim dimensions
                *statePrL,
                *discountL;

        double  *discount = NULL;

        int     nFact = mNbFactor,      // Used in some macros.
                i, j0, j1,              // factor indices
                n,
                Mid;                    // switch bet lo and hi Q


    try {

        // Degenerate the state price slice from nDim to nIRDim, which 
        // is the sum over the higher nDim-nIRDim dimensions
        statePrIR = sliceNew(nIRDim);
        sliceDegenerate(t, nDim, nIRDim, statePr, statePrIR);

        // Expand the discount slice from nIRDim to nDim, which 
        // are constant in the higher nDim-nIRDim dimensions
        discount = sliceNew(nDim);
        sliceExpand(t, nDim, nIRDim, discountIR, discount);
        
        // Previous period jump size.
        //
        du = sqrt(JUMPCOEFF * LengthJ[t-1]);

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
        


        // Forward spread
        fwdSpread = GetBasisFwdSpread(currDate);

        // Avoid zero spread, which is allowed in normal case but
        // canceled out after the mapping.  Set the minimum to 0.01bp
        SHIFT_ZERO(fwdSpread);

        // Set up backbone
        //
        FwdSpreadA  = fwdSpread / (1. + mBSFSh);


        // Vol backbone factor
        //
        VolBbq = GetBasisSpreadVolBbq(t);


        // The base shift in the spread distribution is defined in 
        // such a way that spread distribution at X=0 ALWAYS corresponds 
        // to the forward spread. 
        //
        QSh = (mBSFSh > 0) ? mBSQHi : mBSQLo;
        if (fabs(QSh) > QCUTOFF) 
        {
                baseSh = log(1. + QSh * mBSFSh) / (QSh * VolBbq);
        }
        else
        {
                baseSh = mBSFSh / VolBbq;
        }

        // Actual 2q switch point in X the space
        //
        xSwitch = Zt + baseSh;


        // Calculate the constant coefficients for P[n].
        //
        CQLo[1] = CQHi[1] = 1e0;
        for(n=2; n<=NPOLY; n++)
        {
                CQLo[n] = CQLo[n-1] * VolBbq * mBSQLo/(double)n;
                CQHi[n] = CQHi[n-1] * VolBbq * mBSQHi/(double)n;
        }


        // I. Compute the expected value of index (Libor) rate
        // in nIRDim dimensional slice.

        // Calculate the expected value of index (Libor) rate 
        // in nIRDim dimensional slice
        eIdxRate = 0.0e0;

        switch (nIRDim) {
        case 1:
                discountIRL  = discountIR  + NodeOffset(1, 0, 0, t);
                statePrIRL   = statePrIR   + NodeOffset(1, 0, 0, t);
                liborRateIRL = liborRateIR + NodeOffset(1, 0, 0, t);

                for (j0 = mBottom1[t]; j0 <= mTop1[t]; j0++) {
                        eIdxRate += statePrIRL[j0]  
                                  * liborRateIRL[j0]
                                  * discountIRL[j0];
                }

                break;

        case 2:
                for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                {
                        discountIRL  = discountIR  + NodeOffset(2, j0, 0, t);
                        statePrIRL   = statePrIR   + NodeOffset(2, j0, 0, t);
                        liborRateIRL = liborRateIR + NodeOffset(2, j0, 0, t);

                        for (i = mBottom2[t][j0]; i <= mTop2[t][j0]; i++) {
                                eIdxRate += statePrIRL[i]  
                                           * liborRateIRL[i]
                                           * discountIRL[i];
                        }  // for i

                }  // for j0
        
                break;

        case 3:
                for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
                {
                        discountIRL  = discountIR  + NodeOffset(3, j0, j1, t);
                        statePrIRL   = statePrIR   + NodeOffset(3, j0, j1, t);
                        liborRateIRL = liborRateIR + NodeOffset(3, j0, j1, t);

                        for (i = mBottom3[t][j0][j1]; i <= mTop3[t][j0][j1]; i++) 
                        {
                                eIdxRate += statePrIRL[i]  
                                           * liborRateIRL[i]
                                           * discountIRL[i];
                        }
                }   // for j0, j1

                break;

        default:
                throw KFailure("%s: N/A for nIRDim=%d.\n", routine, nIRDim);
        }


        // II. Compute the coefficients of Eta-polynomials associated
        // with drift of spreads in nDim dimensional slice.

        // Reset polynomial expansion
        //
        for (i = 0; i <= NPOLY; i++) P[i] = 0e0;

        uLo = FwdSpreadA * (1. - mBSQLo); 
        uHi = FwdSpreadA * (1. - mBSQHi); 
    

        // Calculate the coefficients of Eta-polynomials for
        // solving the drift of spreads in nDim dimensional slice
        switch (nDim) {
        case 1:

                // 1-D calibration
                Zidx = xSwitch; 
                Mid  = (int) ceil(-Zidx / pJump00) - 1;
                Mid  = Min( Max(Mid, mBottom1[t] - 1), mTop1[t]);

                // Compute coefficients

                discountL = discount + NodeOffset(1, 0, 0, t);
                statePrL  = statePr  + NodeOffset(1, 0, 0, t);
                spreadL   = spread   + NodeOffset(1, 0, 0, t);

                for (j0 = mBottom1[t]; j0 <= Mid; j0++) {
                        x      = statePrL[j0] * discountL[j0];
                        alpha  = spreadL[j0] * mBSQLo + uLo;

                        // Add backbone 
                        alpha *= VolBbq;
                        
                        P[0] += spreadL[j0] * x;
                        for (n = 1; n <= NPOLY; n++)
                                P[n] += alpha * x * CQLo[n];

                }
                for (j0 = Mid + 1; j0 <= mTop1[t]; j0++) {
                        x      = statePrL[j0] * discountL[j0];
                        alpha  = spreadL[j0] * mBSQHi + uHi;

                        // Add backbone 
                        alpha *= VolBbq;

                        P[0] += spreadL[j0] * x;
                        for (n = 1; n <= NPOLY; n++)
                                P[n] += alpha * x * CQHi[n];
                }

                break;

        case 2:
                for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                {
                        // Find switch point index
                        Zidx  = xSwitch 
                            + j0 * (pJump00 + pJump10);
                        Mid   = (int) ceil(-Zidx / pJump11) - 1;
                        Mid   = Min(Max(Mid, mBottom2[t][j0] - 1),
                                        mTop2[t][j0]);
                        //Zidx += (pJump11) * mBottom2[t][j0];

                        discountL = discount + NodeOffset(2, j0, 0, t);
                        statePrL  = statePr  + NodeOffset(2, j0, 0, t);
                        spreadL   = spread   + NodeOffset(2, j0, 0, t);

                        // Compute coefficients

                        // LEFT part of distribution
                        for (i = mBottom2[t][j0]; i <= Mid; i++) {
                                x      = statePrL[i] * discountL[i];
                                alpha  = spreadL[i] * mBSQLo + uLo;
                        
                                // Add backbone 
                                alpha *= VolBbq;

                                P[0] += spreadL[i] * x;
                                for (n = 1; n <= NPOLY; n++)
                                        P[n] += alpha * x * CQLo[n];
                        }

                        // RIGHT part of distribution
                        for (i = Mid + 1; i <= mTop2[t][j0]; i++) {
                                x      = statePrL[i] * discountL[i];
                                alpha  = spreadL[i] * mBSQHi + uHi;

                                // Add backbone 
                                alpha *= VolBbq;
                        
                                P[0] += spreadL[i] * x;
                                for (n = 1; n <= NPOLY; n++)
                                        P[n] += alpha * x * CQHi[n];
                        }
                }  // for j0
        
                break;

        case 3:
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
                        //Zidx += (pJump22) * mBottom3[t][j0][j1];
            
                        discountL = discount + NodeOffset(3, j0, j1, t);
                        statePrL  = statePr  + NodeOffset(3, j0, j1, t);
                        spreadL   = spread   + NodeOffset(3, j0, j1, t);

                        // Compute coefficients

                        // LEFT part of distribution
                        for (i = mBottom3[t][j0][j1]; i <= Mid; i++) {
                                x      = statePrL[i] * discountL[i];
                                alpha  = spreadL[i] * mBSQLo + uLo;

                                // Add backbone 
                                alpha *= VolBbq;
                        
                                P[0] += spreadL[i] * x;
                                for (n = 1; n <= NPOLY; n++)
                                        P[n] += alpha * x * CQLo[n];
                        }

                        // RIGHT part of distribution
                        for (i = Mid + 1; i <= mTop3[t][j0][j1]; i++) {
                                x      = statePrL[i] * discountL[i];
                                alpha  = spreadL[i] * mBSQHi + uHi;
                        
                                // Add backbone 
                                alpha *= VolBbq;
                        
                                P[0] += spreadL[i] * x;
                                for (n = 1; n <= NPOLY; n++)
                                        P[n] += alpha * x * CQHi[n];
                        }
                }   // for j0, j1

                break;

        default:
                throw KFailure("%s: N/A for nDim=%d.\n", routine, nDim);
        }


        // Solve for zero price at t+1
        if (mBSType == SUB_SPREAD)         /* Basis = Libor - S */
        {
            P[0] -= eIdxRate;
            P[0] += zeroPrice * fwdBasisRate;
        }
        else if (mBSType == ADD_SPREAD)    /* Basis = Libor + S */
        {
            P[0] += eIdxRate;
            P[0] -= zeroPrice * fwdBasisRate;
        }
        else
            throw KFailure("%s: invalid spread type (%d).\n",
                            routine,
                            mBSType);

                                        
        if (debugLevel > DEBUG_LEVEL_DRIFT) {
                dppLog << endl;
 
                dppLog << format("=========== %s: ============", routine);
                dppLog << endl;
 
                for (i=0; i<=NPOLY; ++i)
                        dppLog << "P[" << i << "] = " << setprecision(12)
                               << P[i] << endl;
        }


        IF_FAILED_THROW( DrlNRPoly(
                0e0,
                P, 
                NPOLY,
                Eta));


        if (fabs(*Eta*VolBbq) > 1e0) {
                throw KFailure("%s: problem in calculation "
                        "of drift (Eta = %lf).\n",
                        routine, *Eta);
        }


        if (debugLevel > DEBUG_LEVEL_DRIFT) {
                dppLog << endl;
                dppLog << "Center offset: " <<  *Eta << endl;
                dppLog << endl;
        }


        // Free tmp memory
        sliceDelete(statePrIR);
        sliceDelete(discount);

        statePrIRL   = NULL;
        liborRateIRL = NULL;
        discountIRL  = NULL;

        spreadL   = NULL;
        statePrL  = NULL;
        discountL = NULL;

    }
    catch (KFailure) {
        // Free tmp memory
        sliceDelete(statePrIR);
        sliceDelete(discount);

        statePrIRL   = NULL;
        liborRateIRL = NULL;
        discountIRL  = NULL;

        spreadL   = NULL;
        statePrL  = NULL;
        discountL = NULL;


        throw KFailure("%s: failed.\n", routine);
    }

}



//--------------------------------------------------------------
// Calibrate the basis offset at idxDev so that 1-period basis rate
// discounted at the IR curve equals the 1-period forward basis rate.
// mBSType = SUB_SPREAD or ADD_SPREAD.
// Solving the shift Zt in NR directly by computing F(X) and F'(X)
//

void
KBirTree::SolveBasisOffsetSpreadNew(
        int    tpIdx,           // (I) time point index
        double fwdBasisRate,    // (I) basis forward rate
        double *liborRateIR,    // (I) libor rates on the time slice
        double *statePr,        // (I) state prices on the time slice
        double *discountIR,     // (I) discount between t and T
        double zeroPrice,       // (I) discount zero at T
        double *Zt)             // (I/O) Shift at t
{
static        char        routine[] = "KBirTree::SolveBasisOffsetSpreadNew";

        double  du;                     // grid step

        double  pJump00,
                pJump10, pJump11,
                pJump20, pJump21, pJump22;

        int     t = tpIdx;              // more convenient ...

        TDate   currDate = TPDate(t);   // Current date
        
        int     nIRDim = this->mIRDim;  // IR dimensions
        int     nDim = mIRDim + mBSDim; // tree dimensions

        double  QSh;                    // Used to calculate initial Zt
        double  baseSh;                 // Base shift with distribution
                                        // at X=0 corresponds to forward spread

        double  xSwitch;                // Actual 2q switch point in X space

        double  VolBbq;                 // Sigma used in bone mapping

        double  fwdSpread,              // fwd spread
                FwdSpreadA,             // rescaled fwd spread / Fa
                Zidx,                   // Zt index j0, j1, j2 adjusted
                eIdxRate,               // expected value of index rate 
                x,
                uLo,
                uHi;

        double  MLo, MHi,               // constants used in grid calc
                SLo, SHi,
                Grid,                   // grid step
                dGrid,
                spread, 
                DQ,
                rJumpHi,                // rate jump size high and low
                rJumpLo, 
                rJumpHi2,               // rate jump size high and low
                rJumpLo2, 
                rJumpHi5,               // rate jump size high and low
                rJumpLo5; 


        double  *statePrIR = NULL,      // variables in nIRDim dimensions
                *statePrIRL,
                *liborRateIRL,
                *discountIRL;

        double  *statePrL,              // variables in nDim dimensions
                *discountL;

        double  *discount = NULL;


        int     nFact = mNbFactor,      // Used in some macros.
                i, j0, j1,              // factor indices
                n,
                Mid;                    // switch bet lo and hi Q

        double  Zt0 = *Zt;              // Initial guestt

        double  F, dF;                  // Root finding of F and dF/dx.
        double  Eta;                    // Change of Zt

        int     iterNR;

        const   int     MAXITER = 100;
        const   double  MAXXERR = 1e-7;

    try {

        // Degenerate the state price slice from nDim to nIRDim, which 
        // is the sum over the higher nDim-nIRDim dimensions
        statePrIR = sliceNew(nIRDim);
        sliceDegenerate(t, nDim, nIRDim, statePr, statePrIR);

        // Expand the discount slice from nIRDim to nDim, which 
        // are constant in the higher nDim-nIRDim dimensions
        discount = sliceNew(nDim);
        sliceExpand(t, nDim, nIRDim, discountIR, discount);

        // Previous period jump size.
        //
        du = sqrt(JUMPCOEFF * LengthJ[t-1]);

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
        


        // Forward spread
        fwdSpread = GetBasisFwdSpread(currDate);

        // Avoid zero spread, which is allowed in normal case but
        // canceled out after the mapping.  Set the minimum to 0.01bp
        SHIFT_ZERO(fwdSpread);

        // Avoid special treatment for normal.
        SHIFT_ZERO(mBSQLo);
        SHIFT_ZERO(mBSQHi);

        // Set up backbone
        //
        FwdSpreadA  = fwdSpread / (1. + mBSFSh);


        // Vol backbone factor
        //
        VolBbq = GetBasisSpreadVolBbq(t);


        // The base shift in the spread distribution is defined in 
        // such a way that spread distribution at X=0 ALWAYS corresponds 
        // to the forward spread. 
        //
        QSh = (mBSFSh > 0) ? mBSQHi : mBSQLo;
        if (fabs(QSh) > QCUTOFF) 
        {
                baseSh = log(1. + QSh * mBSFSh) / (QSh * VolBbq);
        }
        else
        {
                baseSh = mBSFSh / VolBbq;
        }

        // Actual 2q switch point in the X space
        // This will be used as initial guess for Zt
        xSwitch = Zt0 + baseSh;



        //========================================================
        //
        //  ONLY 1 + 1 IS SUPPORTED AND HAS BEEN TESTED!!!
        //
        //========================================================


        // I. Compute the expected value of index (Libor) rate
        // in nIRDim dimensional slice.

        // Calculate the expected value of index (Libor) rate 
        // in nIRDim dimensional slice
        eIdxRate = 0.0e0;

        switch (nIRDim) {
        case 1:
                discountIRL  = discountIR  + NodeOffset(1, 0, 0, t);
                statePrIRL   = statePrIR   + NodeOffset(1, 0, 0, t);
                liborRateIRL = liborRateIR + NodeOffset(1, 0, 0, t);

                for (j0 = mBottom1[t]; j0 <= mTop1[t]; j0++) {
                    eIdxRate += statePrIRL[j0]  
                             * liborRateIRL[j0]
                             * discountIRL[j0];
                }

                break;

        case 2:
                for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                {
                    discountIRL  = discountIR  + NodeOffset(2, j0, 0, t);
                    statePrIRL   = statePrIR   + NodeOffset(2, j0, 0, t);
                    liborRateIRL = liborRateIR + NodeOffset(2, j0, 0, t);

                    for (i = mBottom2[t][j0]; i <= mTop2[t][j0]; i++) {
                        eIdxRate += statePrIRL[i]  
                                 * liborRateIRL[i]
                                 * discountIRL[i];
                    }  // for i

                }  // for j0
        
                break;

        case 3:
                for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
                {
                    discountIRL  = discountIR  + NodeOffset(3, j0, j1, t);
                    statePrIRL   = statePrIR   + NodeOffset(3, j0, j1, t);
                    liborRateIRL = liborRateIR + NodeOffset(3, j0, j1, t);

                    for (i = mBottom3[t][j0][j1]; i <= mTop3[t][j0][j1]; i++) 
                    {
                        eIdxRate += statePrIRL[i]  
                                 * liborRateIRL[i]
                                 * discountIRL[i];
                    }
                }   // for j0, j1

                break;

        default:
                throw KFailure("%s: N/A for nIRDim=%d.\n", routine, nIRDim);
        }



        // II. Solve for the offset shift Zt for basis spread
        // Using Newton-Ralphson method
        //
        // F = spread*discount*StatPrice - Libor*discount*StatPrice + B*Z
        //

        for (iterNR = 0; iterNR < MAXITER; iterNR++)
        {

            //
            // Compute F(X) and F'(X) given xSwitch
            //

            // Reset F and F' 
            //
            F = dF = 0e0;

            uLo = FwdSpreadA * (1. - mBSQLo); 
            uHi = FwdSpreadA * (1. - mBSQHi); 
            MLo = FwdSpreadA / mBSQLo;
            MHi = FwdSpreadA / mBSQHi;
            SLo = FwdSpreadA - FwdSpreadA / mBSQLo;
            SHi = FwdSpreadA - FwdSpreadA / mBSQHi;
            DQ   = FwdSpreadA * VolBbq;

            switch (nDim) {
            case 1:

                rJumpLo = exp(mBSQLo * VolBbq * pJump00);
                rJumpHi = exp(mBSQHi * VolBbq * pJump00);

                // 1-D calibration
                Zidx = xSwitch; 
                Mid  = (int) ceil(-Zidx / pJump00) - 1;
                Mid  = Min( Max(Mid, mBottom1[t] - 1), mTop1[t]);

                // Compute coefficients

                discountL = discount + NodeOffset(1, 0, 0, t);
                statePrL  = statePr  + NodeOffset(1, 0, 0, t);

                // LEFT part of distribution
                Zidx += (pJump00) * mBottom1[t];
                Grid  = MLo * exp (mBSQLo * VolBbq * Zidx);
                dGrid = DQ  * exp (mBSQLo * VolBbq * Zidx);

                for (j0 = mBottom1[t]; j0 <= Mid; j0++) {
                        x       = statePrL[j0] * discountL[j0];
                        spread  = SLo + Grid;

                        Grid  *= rJumpLo;
                        dGrid *= rJumpLo;

                        F  += spread * x;
                        dF += dGrid * x;

                }

                // RIGHT part of distribution
                Zidx  = xSwitch + (pJump00) * (Mid + 1);
                Grid  = MHi * exp (mBSQHi * VolBbq * Zidx);
                dGrid = DQ  * exp (mBSQHi * VolBbq * Zidx);

                for (j0 = Mid + 1; j0 <= mTop1[t]; j0++) {
                        x      = statePrL[j0] * discountL[j0];
                        spread  = SHi + Grid;

                        Grid  *= rJumpHi;
                        dGrid *= rJumpHi;

                        F  += spread * x;
                        dF += dGrid * x;
                }

                break;

            case 2:

                rJumpLo2 = exp (mBSQLo * VolBbq * pJump11);
                rJumpHi2 = exp (mBSQHi * VolBbq * pJump11); 

                switch (mBSDim) {
                case 1:     // nIRDim = 1, mBSDim = 1;  THIS IS THE MODE SUPPORTED

                    for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                    {
                        // Find switch point index
                        Zidx  = xSwitch + j0 * pJump10;
                        Mid   = (int) ceil(-Zidx / pJump11) - 1;
                        Mid   = Min(Max(Mid, mBottom2[t][j0] - 1),
                                        mTop2[t][j0]);

                        discountL = discount + NodeOffset(2, j0, 0, t);
                        statePrL  = statePr  + NodeOffset(2, j0, 0, t);

                        // LEFT part of distribution
                        Zidx += (pJump11) * mBottom2[t][j0];
                        Grid  = MLo * exp (mBSQLo * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQLo * VolBbq * Zidx);

                        for (i = mBottom2[t][j0]; i <= Mid; i++) {
                                x      = statePrL[i] * discountL[i];
                                spread = SLo + Grid;
                        
                                Grid  *= rJumpLo2;
                                dGrid *= rJumpLo2;

                                F  += spread * x;
                                dF += dGrid * x;
                        }


                        // RIGHT part of distribution
                        Zidx  = xSwitch + j0 * pJump10 + (Mid + 1) * pJump11;
                        Grid  = MHi * exp (mBSQHi * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQHi * VolBbq * Zidx);

                        for (i = Mid + 1; i <= mTop2[t][j0]; i++) {
                                x      = statePrL[i] * discountL[i];
                                spread = SHi + Grid;

                                Grid  *= rJumpHi2;
                                dGrid *= rJumpHi2;

                                F  += spread * x;
                                dF += dGrid * x;
                        }
                    }  // for j0
        
                    break;

                case 2:     // nIRDim = 0, mBSDim = 2;
                    for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                    {
                        // Find switch point index
                        Zidx  = xSwitch + j0 * (pJump00 + pJump10);
                        Mid   = (int) ceil(-Zidx / pJump11) - 1;
                        Mid   = Min(Max(Mid, mBottom2[t][j0] - 1),
                                        mTop2[t][j0]);

                        discountL = discount + NodeOffset(2, j0, 0, t);
                        statePrL  = statePr  + NodeOffset(2, j0, 0, t);


                        // LEFT part of distribution
                        Zidx += (pJump11) * mBottom2[t][j0];
                        Grid  = MLo * exp (mBSQLo * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQLo * VolBbq * Zidx);

                        for (i = mBottom2[t][j0]; i <= Mid; i++) {
                                x      = statePrL[i] * discountL[i];
                                spread = SLo + Grid;
                        
                                Grid  *= rJumpLo2;
                                dGrid *= rJumpLo2;

                                F  += spread * x;
                                dF += dGrid * x;
                        }


                        // RIGHT part of distribution
                        Zidx  = xSwitch + j0 * (pJump00 + pJump10) 
                                        + (Mid + 1) * pJump11;
                        Grid  = MHi * exp (mBSQHi * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQHi * VolBbq * Zidx);

                        for (i = Mid + 1; i <= mTop2[t][j0]; i++) {
                                x      = statePrL[i] * discountL[i];
                                spread = SHi + Grid;

                                Grid  *= rJumpHi2;
                                dGrid *= rJumpHi2;

                                F  += spread * x;
                                dF += dGrid * x;
                        }
                    }  // for j0
        
                    break;
        
                default:
                    throw KFailure("%s: invalid factor dimensions."
                               " nDim=%d, nBSDim=%d.\n", routine, nDim, mBSDim);
                }

                break;

            case 3:

                rJumpLo5 = exp (mBSQLo * VolBbq * pJump22);
                rJumpHi5 = exp (mBSQHi * VolBbq * pJump22); 

                switch (mBSDim) {
                case 1:     // nIRDim = 2, nBSDim = 1;

                    for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                    for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
                    {
                        Zidx  = xSwitch
                              + j0 * (pJump20)
                              + j1 * (pJump21);
                        Mid   = (int) ceil(-Zidx / pJump22) - 1;
                        Mid   = Min(Max(Mid, mBottom3[t][j0][j1] - 1),
                                        mTop3[t][j0][j1]);

                        discountL = discount + NodeOffset(3, j0, j1, t);
                        statePrL  = statePr  + NodeOffset(3, j0, j1, t);

                        // LEFT part of distribution
                        Zidx += (pJump22) * mBottom3[t][j0][j1];
                        Grid  = MLo * exp (mBSQLo * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQLo * VolBbq * Zidx);

                        for (i = mBottom3[t][j0][j1]; i <= Mid; i++) {
                                x      = statePrL[i] * discountL[i];
                                spread = SLo + Grid;
                        
                                Grid  *= rJumpLo5;
                                dGrid *= rJumpLo5;

                                F  += spread * x;
                                dF += dGrid * x;

                        }

                        // RIGHT part of distribution
                        Zidx = xSwitch
                             + j0 * (pJump20)
                             + j1 * (pJump21)
                             + (Mid+1) * (pJump22);
                        Grid  = MHi * exp (mBSQHi * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQHi * VolBbq * Zidx);

                        for (i = Mid + 1; i <= mTop3[t][j0][j1]; i++) {
                                x      = statePrL[i] * discountL[i];
                                spread = SHi + Grid;
                        
                                Grid  *= rJumpHi5;
                                dGrid *= rJumpHi5;

                                F  += spread * x;
                                dF += dGrid * x;
                        }
                    }

                    break;

                case 2:     // nIRDim = 1, nBSDim = 2;

                    for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                    for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
                    {
                        Zidx  = xSwitch
                              + j0 * (pJump10 + pJump20)
                              + j1 * (pJump11 + pJump21);
                        Mid   = (int) ceil(-Zidx / pJump22) - 1;
                        Mid   = Min(Max(Mid, mBottom3[t][j0][j1] - 1),
                                        mTop3[t][j0][j1]);

                        discountL = discount + NodeOffset(3, j0, j1, t);
                        statePrL  = statePr  + NodeOffset(3, j0, j1, t);

                        // LEFT part of distribution
                        Zidx += (pJump22) * mBottom3[t][j0][j1];
                        Grid  = MLo * exp (mBSQLo * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQLo * VolBbq * Zidx);

                        for (i = mBottom3[t][j0][j1]; i <= Mid; i++) {
                                x      = statePrL[i] * discountL[i];
                                spread = SLo + Grid;
                        
                                Grid  *= rJumpLo5;
                                dGrid *= rJumpLo5;

                                F  += spread * x;
                                dF += dGrid * x;

                        }

                        // RIGHT part of distribution
                        Zidx = xSwitch
                             + j0 * (pJump10 + pJump20)
                             + j1 * (pJump11 + pJump21)
                             + (Mid+1) * (pJump22);
                        Grid  = MHi * exp (mBSQHi * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQHi * VolBbq * Zidx);

                        for (i = Mid + 1; i <= mTop3[t][j0][j1]; i++) {
                                x      = statePrL[i] * discountL[i];
                                spread = SHi + Grid;
                        
                                Grid  *= rJumpHi5;
                                dGrid *= rJumpHi5;

                                F  += spread * x;
                                dF += dGrid * x;
                        }
                    }

                    break;

                case 3:     // nIRDim = 0, nBSDim = 3;
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

                        discountL = discount + NodeOffset(3, j0, j1, t);
                        statePrL  = statePr  + NodeOffset(3, j0, j1, t);

                        // LEFT part of distribution
                        Zidx += (pJump22) * mBottom3[t][j0][j1];
                        Grid  = MLo * exp (mBSQLo * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQLo * VolBbq * Zidx);

                        for (i = mBottom3[t][j0][j1]; i <= Mid; i++) {
                                x      = statePrL[i] * discountL[i];
                                spread = SLo + Grid;
                        
                                Grid  *= rJumpLo5;
                                dGrid *= rJumpLo5;

                                F  += spread * x;
                                dF += dGrid * x;
                        }

                        // RIGHT part of distribution
                        Zidx = xSwitch
                             + j0 * (pJump00 + pJump10 + pJump20)
                             + j1 * (pJump11 + pJump21)
                             + (Mid+1) * (pJump22);
                        Grid  = MHi * exp (mBSQHi * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQHi * VolBbq * Zidx);

                        for (i = Mid + 1; i <= mTop3[t][j0][j1]; i++) {
                                x      = statePrL[i] * discountL[i];
                                spread = SHi + Grid;
                        
                                Grid  *= rJumpHi5;
                                dGrid *= rJumpHi5;

                                F  += spread * x;
                                dF += dGrid * x;
                            }
                        }   // for j0, j1

                    break;

                default:
                    throw KFailure("%s: invalid factor dimensions."
                               " nDim=%d, nBSDim=%d.\n", routine, nDim, mBSDim);
                }

                break;    // nDim = 3

            default:
                throw KFailure("%s: N/A for nDim=%d.\n", routine, nDim);

            }    // end of switch nDim


            // Add constant terms
            //
            if (mBSType == SUB_SPREAD)         /* Basis = Libor - S */
            {
                F -= eIdxRate;
                F += zeroPrice * fwdBasisRate;
            }
            else if (mBSType == ADD_SPREAD)    /* Basis = Libor + S */
            {
                F += eIdxRate;
                F -= zeroPrice * fwdBasisRate;
            }
            else
                throw KFailure("%s: invalid spread type (%d).\n",
                                routine,
                                mBSType);

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
                *Zt = xSwitch - baseSh;
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
                dppLog << "Center offset: " <<  *Zt - Zt0 << endl;
                dppLog << endl;
        }


        // Free tmp memory
        sliceDelete(statePrIR);
        sliceDelete(discount);

        statePrIRL   = NULL;
        liborRateIRL = NULL;
        discountIRL  = NULL;

        statePrL  = NULL;
        discountL = NULL;

    }
    catch (KFailure) {
        // Free tmp memory
        sliceDelete(statePrIR);
        sliceDelete(discount);

        statePrIRL   = NULL;
        liborRateIRL = NULL;
        discountIRL  = NULL;

        statePrL  = NULL;
        discountL = NULL;


        throw KFailure("%s: failed.\n", routine);
    }

}



//--------------------------------------------------------------
// Calibrate the basis offset at idxDev so that 1-period basis rate
// discounted at the IR curve equals the 1-period forward basis rate.
// mBSType = PERCENT.

void
KBirTree::SolveBasisOffsetPercent(
        int    tpIdx,           // (I) time point index
        double fwdBasisRate,    // (I) basis forward rate
        double *spread,         // (I) spread on the time slice
        double *liborRateIR,    // (I) libor rates on the time slice
        double *statePr,        // (I) state prices on the time slice
        double *discountIR,     // (I) discount between t and T
        double zeroPrice,       // (I) discount zero at T
        double Zt,              // (I) base shift in Taylor expansion
        double *Eta)            // (O) incremental offset
{
static        char        routine[] = "KBirTree::SolveBasisOffsetPercent";

        double  du;                     // grid step

        double  pJump00,
                pJump10, pJump11,
                pJump20, pJump21, pJump22;


        int     t = tpIdx;              // more convenient ...
        
        TDate   currDate = TPDate(t);   // Current date

        int     nIRDim = this->mIRDim;  // IR dimensions
        int     nDim = mIRDim + mBSDim; // tree dimensions

        double  QSh;                    // Used to calculate initial Zt
        double  baseSh;                        // Base shift with distribution
                                        // at X=0 corresponds to forward spread

        double  xSwitch;                // Actual 2q switch point in X space

        double  VolBbq;                 // Sigma used in bone mapping

        double  fwdSpread,              // fwd spread
                FwdSpreadA,             // rescaled fwd spread / Fa
                //InsSpread,            // Instantaneous spread
                pJumpM,                 // jump in the last dimension
                Zidx,                   // xSwitch index j0, j1, j2 adjusted
                P[NPOLY+1],             // expansion
                CQLo[NPOLY+1],          // Low coefficients
                CQHi[NPOLY+1],          // Hight coefficients
                alpha,                  // Coeff of Eta-polynomial
                x,
                uLo,
                uHi;


        double  *tmpSliceIR = NULL;     // variables in nIRDim dimensions

        double  *spreadL,               // variables in nDim dimensions
                *statePrL,
                *tmpSliceL;

        double  *tmpSlice = NULL;

        int     nFact = mNbFactor,      // Used in some macros.
                i, j0, j1,              // factor indices
                n,
                Mid;                    // switch bet lo and hi Q


    try {

        // Previous period jump size.
        //
        du = sqrt(JUMPCOEFF * LengthJ[t-1]);

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
        

        // Forward spread
        fwdSpread = GetBasisFwdSpread(currDate);

        // Avoid zero spread, which is allowed in normal case but
        // canceled out after the mapping.  Set the minimum to 0.01bp
        SHIFT_ZERO(fwdSpread);

        // Set up backbone
        //
        FwdSpreadA  = fwdSpread / (1. + mBSFSh);


        // Vol backbone factor
        //
        VolBbq = GetBasisSpreadVolBbq(t);


        // The base shift in the spread distribution is defined in 
        // such a way that spread distribution at X=0 ALWAYS corresponds 
        // to the forward spread. 
        //
        QSh = (mBSFSh > 0) ? mBSQHi : mBSQLo;
        if (fabs(QSh) > QCUTOFF) 
        {
                baseSh = log(1. + QSh * mBSFSh) / (QSh * VolBbq);
        }
        else
        {
                baseSh = mBSFSh / VolBbq;
        }

        // Actual 2q switch point in X the space
        //
        xSwitch = Zt + baseSh;


        //
        // Compute the coefficients of Eta-polynomials associated
        // with drift of spreads in nDim dimensional slice.
        //

        // Reset polynomial expnsion
        //
        for (i = 0; i <= NPOLY; i++) P[i] = 0e0;

        // Calculate the constant coefficients for P[n].
        //
        CQLo[1] = CQHi[1] = 1e0;
        for(n=2; n<=NPOLY; n++)
        {
                CQLo[n] = CQLo[n-1] * VolBbq * mBSQLo/(double)n;
                CQHi[n] = CQHi[n-1] * VolBbq * mBSQHi/(double)n;
        }


        uLo = FwdSpreadA * (1. - mBSQLo); 
        uHi = FwdSpreadA * (1. - mBSQHi); 
    
        // Calculate Libor * Discount in nIRDim dimension
        tmpSliceIR = sliceNew(nIRDim);

        sliceUnaryOper( tmpSliceIR,
                        nIRDim,
                        discountIR,
                        COPY);

        sliceUnaryOper( tmpSliceIR,
                        nIRDim,
                        liborRateIR,
                        MULT);
        
                        
        // Expand tmpSlice from nIRDim to nDim, which 
        // are constant in the higher nDim-nIRDim dimensions
        tmpSlice = sliceNew(nDim);
        sliceExpand(t, nDim, nIRDim, tmpSliceIR, tmpSlice);
        

        // Calculate the coefficients of Eta-polynomials for
        // solving the drift of spreads in nDim dimensional slice
        switch (nDim) {
        case 1:
                // jump in highest dimension 
                pJumpM = pJump00;

                // 1-D calibration
                Zidx = xSwitch; 
                Mid  = (int) ceil(-Zidx / pJumpM) - 1;
                Mid  = Min( Max(Mid, mBottom1[t] - 1), mTop1[t]);

                // Compute coefficients

                tmpSliceL  = tmpSlice  + NodeOffset(1, 0, 0, t);
                statePrL   = statePr   + NodeOffset(1, 0, 0, t);
                spreadL    = spread    + NodeOffset(1, 0, 0, t);

                // LEFT part of distribution
                for (j0 = mBottom1[t]; j0 <= Mid; j0++) {
                        x      = statePrL[j0] * tmpSliceL[j0];
                        alpha  = spreadL[j0] * mBSQLo + uLo;

                        // Add backbone
                        alpha *= VolBbq;
                        
                        P[0] += spreadL[j0] * x;
                        for (n = 1; n <= NPOLY; n++)
                                P[n] += alpha * x * CQLo[n];
                }

                // RIGHT part of distribution
                for (j0 = Mid + 1; j0 <= mTop1[t]; j0++) {
                        x      = statePrL[j0] * tmpSliceL[j0];
                        alpha  = spreadL[j0] * mBSQHi + uHi;

                        // Add backbone
                        alpha *= VolBbq;

                        P[0] += spreadL[j0] * x;
                        for (n = 1; n <= NPOLY; n++)
                                P[n] += alpha * x * CQHi[n];
                }

                break;

        case 2:
                for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                {
                        // Find switch point index
                        Zidx  = xSwitch 
                            + j0 * (pJump00 + pJump10);
                        Mid   = (int) ceil(-Zidx / pJump11) - 1;
                        Mid   = Min(Max(Mid, mBottom2[t][j0] - 1),
                                        mTop2[t][j0]);
                        //Zidx += (pJump11) * mBottom2[t][j0];

                        tmpSliceL  = tmpSlice  + NodeOffset(2, j0, 0, t);
                        statePrL   = statePr   + NodeOffset(2, j0, 0, t);
                        spreadL    = spread    + NodeOffset(2, j0, 0, t);

                        // Compute coefficients

                        // LEFT part of distribution
                        for (i = mBottom2[t][j0]; i <= Mid; i++) {
                                x = statePrL[i] * tmpSliceL[i];
                                alpha  = spreadL[i] * mBSQLo + uLo;
                        
                                // Add backbone
                                alpha *= VolBbq;

                                P[0] += spreadL[i] * x;
                                for (n = 1; n <= NPOLY; n++)
                                        P[n] += alpha * x * CQLo[n];
                        }

                        // RIGHT part of distribution
                        for (i = Mid + 1; i <= mTop2[t][j0]; i++) {
                                x      = statePrL[i] * tmpSliceL[i];
                                alpha  = spreadL[i] * mBSQHi + uHi;
                        
                                // Add backbone
                                alpha *= VolBbq;

                                P[0] += spreadL[i] * x;
                                for (n = 1; n <= NPOLY; n++)
                                        P[n] += alpha * x * CQHi[n];
                        }
                }  // for j0
        
                break;

        case 3:
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
                        //Zidx += (pJump22) * mBottom3[t][j0][j1];
            
                        tmpSliceL = tmpSlice + NodeOffset(3, j0, j1, t);
                        statePrL  = statePr  + NodeOffset(3, j0, j1, t);
                        spreadL   = spread   + NodeOffset(3, j0, j1, t);

                        // Compute coefficients

                        // LEFT part of distribution
                        for (i = mBottom3[t][j0][j1]; i <= Mid; i++) {
                                x      = statePrL[i] * tmpSliceL[i];
                                alpha  = spreadL[i] * mBSQLo + uLo;
                        
                                // Add backbone
                                alpha *= VolBbq;

                                P[0] += spreadL[i] * x;
                                for (n = 1; n <= NPOLY; n++)
                                        P[n] += alpha * x * CQLo[n];
                        }

                        // RIGHT part of distribution
                        for (i = Mid + 1; i <= mTop3[t][j0][j1]; i++) {
                                x      = statePrL[i] * tmpSliceL[i];
                                alpha  = spreadL[i] * mBSQHi + uHi;
                        
                                // Add backbone
                                alpha *= VolBbq;

                                P[0] += spreadL[i] * x;
                                for (n = 1; n <= NPOLY; n++)
                                        P[n] += alpha * x * CQHi[n];

                        }
                }   // for j0, j1

                break;

        default:
                throw KFailure("%s: N/A for nDim=%d.\n", routine, nDim);
        }


        // Solve for zero price at t+1
        P[0] -= zeroPrice * fwdBasisRate;
                                        
        if (debugLevel > DEBUG_LEVEL_DRIFT) {
                dppLog << endl;
 
                dppLog << format("=========== %s: ============", routine);
                dppLog << endl;
 
                for (i=0; i<=NPOLY; ++i)
                        dppLog << "P[" << i << "] = " << setprecision(12)
                               << P[i] << endl;
        }

        IF_FAILED_THROW( DrlNRPoly(
                0e0,
                P, 
                NPOLY,
                Eta));


        if (fabs(*Eta*VolBbq) > 1e0) {
                throw KFailure("%s: problem in calculation "
                        "of drift (Eta = %lf).\n",
                        routine, *Eta);
        }




        if (debugLevel > DEBUG_LEVEL_DRIFT) {
                dppLog << endl;
                dppLog << "Center offset: " <<  *Eta << endl;
                dppLog << endl;
        }


        // Free tmp memory
        sliceDelete(tmpSliceIR);
        sliceDelete(tmpSlice);

        spreadL  = NULL;        
        statePrL = NULL;
        tmpSliceL= NULL;

    }
    catch (KFailure) {
        // Free tmp memory
        sliceDelete(tmpSliceIR);
        sliceDelete(tmpSlice);

        spreadL  = NULL;        
        statePrL = NULL;
        tmpSliceL= NULL;

        throw KFailure("%s: failed.\n", routine);
    }

}




//--------------------------------------------------------------
// Calibrate the basis offset at idxDev so that 1-period basis rate
// discounted at the IR curve equals the 1-period forward basis rate.
// Solving the shift Zt in NR directly by computing F(X) and F'(X)
// mBSType = PERCENT.
//


void
KBirTree::SolveBasisOffsetPercentNew(
        int    tpIdx,           // (I) time point index
        double fwdBasisRate,    // (I) basis forward rate
        double *liborRateIR,    // (I) libor rates on the time slice
        double *statePr,        // (I) state prices on the time slice
        double *discountIR,     // (I) discount between t and T
        double zeroPrice,       // (I) discount zero at T
        double *Zt)             // (I/O) Shift at t
{
static        char        routine[] = "KBirTree::SolveBasisOffsetPercentNew";

        double  du;                     // grid step

        double  pJump00,
                pJump10, pJump11,
                pJump20, pJump21, pJump22;

        int     t = tpIdx;              // more convenient ...

        TDate   currDate = TPDate(t);   // Current date
        
        int     nIRDim = this->mIRDim;  // IR dimensions
        int     nDim = mIRDim + mBSDim; // tree dimensions

        double  QSh;                    // Used to calculate initial Zt
        double  baseSh;                 // Base shift with distribution
                                        // at X=0 corresponds to forward spread

        double  xSwitch;                // Actual 2q switch point in X space

        double  VolBbq;                 // Sigma used in bone mapping

        double  fwdSpread,              // fwd spread
                FwdSpreadA,             // rescaled fwd spread / Fa
                Zidx,                   // Zt index j0, j1, j2 adjusted
                eIdxRate,               // expected value of index rate 
                x,
                uLo,
                uHi;

        double  MLo, MHi,               // constants used in grid calc
                SLo, SHi,
                Grid,                   // grid step
                dGrid,
                spread, 
                DQ,
                rJumpHi,                // rate jump size high and low
                rJumpLo, 
                rJumpHi2,               // rate jump size high and low
                rJumpLo2, 
                rJumpHi5,               // rate jump size high and low
                rJumpLo5; 


        double  *tmpSliceIR = NULL,      // variables in nIRDim dimensions
                *tmpSlice   = NULL;

        double  *statePrL,              // variables in nDim dimensions
                *tmpSliceL;


        int     nFact = mNbFactor,      // Used in some macros.
                i, j0, j1,              // factor indices
                n,
                Mid;                    // switch bet lo and hi Q

        double  Zt0 = *Zt;              // Initial guestt

        double  F, dF;                  // Root finding of F and dF/dx.
        double  Eta;                    // Change of Zt

        int     iterNR;

        const   int     MAXITER = 100;
        const   double  MAXXERR = 1e-7;

    try {


        // Previous period jump size.
        //
        du = sqrt(JUMPCOEFF * LengthJ[t-1]);

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
        


        // Forward spread
        fwdSpread = GetBasisFwdSpread(currDate);

        // Avoid zero spread, which is allowed in normal case but
        // canceled out after the mapping.  Set the minimum to 0.01bp
        SHIFT_ZERO(fwdSpread);

        // Avoid special treatment for normal.
        SHIFT_ZERO(mBSQLo);
        SHIFT_ZERO(mBSQHi);

        // Set up backbone
        //
        FwdSpreadA  = fwdSpread / (1. + mBSFSh);


        // Vol backbone factor
        //
        VolBbq = GetBasisSpreadVolBbq(t);


        // The base shift in the spread distribution is defined in 
        // such a way that spread distribution at X=0 ALWAYS corresponds 
        // to the forward spread. 
        //
        QSh = (mBSFSh > 0) ? mBSQHi : mBSQLo;
        if (fabs(QSh) > QCUTOFF) 
        {
                baseSh = log(1. + QSh * mBSFSh) / (QSh * VolBbq);
        }
        else
        {
                baseSh = mBSFSh / VolBbq;
        }

        // Actual 2q switch point in the X space
        // This will be used as initial guess for Zt

        xSwitch = Zt0 + baseSh;



        //========================================================
        //
        //  ONLY 1 + 1 IS SUPPORTED AND HAS BEEN TESTED!!!
        //
        //========================================================


        // Solve for the offset shift Zt for basis spread
        // Using Newton-Ralphson method
        //
        // F = spread * (discountIR*Libor) * StatPrice - B*Z
        //                     ||
        //                  tmpSliceIR
        //

        // Calculate Libor * Discount in nIRDim dimension
        tmpSliceIR = sliceNew(nIRDim);

        sliceUnaryOper( tmpSliceIR,
                        nIRDim,
                        discountIR,
                        COPY);

        sliceUnaryOper( tmpSliceIR,
                        nIRDim,
                        liborRateIR,
                        MULT);


        // Expand tmpSlice from nIRDim to nDim, which
        // are constant in the higher nDim-nIRDim dimensions
        tmpSlice = sliceNew(nDim);
        sliceExpand(t, nDim, nIRDim, tmpSliceIR, tmpSlice);


        for (iterNR = 0; iterNR < MAXITER; iterNR++)
        {

            //
            // Compute F(X) and F'(X) given xSwitch
            //

            // Reset F and F' 
            //
            F = dF = 0e0;

            uLo = FwdSpreadA * (1. - mBSQLo); 
            uHi = FwdSpreadA * (1. - mBSQHi); 
            MLo = FwdSpreadA / mBSQLo;
            MHi = FwdSpreadA / mBSQHi;
            SLo = FwdSpreadA - FwdSpreadA / mBSQLo;
            SHi = FwdSpreadA - FwdSpreadA / mBSQHi;
            DQ   = FwdSpreadA * VolBbq;

            switch (nDim) {
            case 1:

                rJumpLo = exp(mBSQLo * VolBbq * pJump00);
                rJumpHi = exp(mBSQHi * VolBbq * pJump00);

                // 1-D calibration
                Zidx = xSwitch; 
                Mid  = (int) ceil(-Zidx / pJump00) - 1;
                Mid  = Min( Max(Mid, mBottom1[t] - 1), mTop1[t]);

                // Compute coefficients

                tmpSliceL = tmpSlice + NodeOffset(1, 0, 0, t);
                statePrL  = statePr  + NodeOffset(1, 0, 0, t);

                // LEFT part of distribution
                Zidx += (pJump00) * mBottom1[t];
                Grid  = MLo * exp (mBSQLo * VolBbq * Zidx);
                dGrid = DQ  * exp (mBSQLo * VolBbq * Zidx);

                for (j0 = mBottom1[t]; j0 <= Mid; j0++) {
                        x       = statePrL[j0] * tmpSliceL[j0];
                        spread  = SLo + Grid;

                        Grid  *= rJumpLo;
                        dGrid *= rJumpLo;

                        F  += spread * x;
                        dF += dGrid * x;

                }

                // RIGHT part of distribution
                Zidx  = xSwitch + (pJump00) * (Mid + 1);
                Grid  = MHi * exp (mBSQHi * VolBbq * Zidx);
                dGrid = DQ  * exp (mBSQHi * VolBbq * Zidx);

                for (j0 = Mid + 1; j0 <= mTop1[t]; j0++) {
                        x      = statePrL[j0] * tmpSliceL[j0];
                        spread  = SHi + Grid;

                        Grid  *= rJumpHi;
                        dGrid *= rJumpHi;

                        F  += spread * x;
                        dF += dGrid * x;
                }

                break;

            case 2:

                rJumpLo2 = exp (mBSQLo * VolBbq * pJump11);
                rJumpHi2 = exp (mBSQHi * VolBbq * pJump11); 

                switch (mBSDim) {
                case 1:     // nIRDim = 1, mBSDim = 1;  THIS IS THE MODE SUPPORTED

                    for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                    {
                        // Find switch point index
                        Zidx  = xSwitch + j0 * pJump10;
                        Mid   = (int) ceil(-Zidx / pJump11) - 1;
                        Mid   = Min(Max(Mid, mBottom2[t][j0] - 1),
                                        mTop2[t][j0]);

                        tmpSliceL = tmpSlice + NodeOffset(2, j0, 0, t);
                        statePrL  = statePr  + NodeOffset(2, j0, 0, t);

                        // LEFT part of distribution
                        Zidx += (pJump11) * mBottom2[t][j0];
                        Grid  = MLo * exp (mBSQLo * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQLo * VolBbq * Zidx);

                        for (i = mBottom2[t][j0]; i <= Mid; i++) {
                                x      = statePrL[i] * tmpSliceL[i];
                                spread = SLo + Grid;
                        
                                Grid  *= rJumpLo2;
                                dGrid *= rJumpLo2;

                                F  += spread * x;
                                dF += dGrid * x;
                        }


                        // RIGHT part of distribution
                        Zidx  = xSwitch + j0 * pJump10 + (Mid + 1) * pJump11;
                        Grid  = MHi * exp (mBSQHi * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQHi * VolBbq * Zidx);

                        for (i = Mid + 1; i <= mTop2[t][j0]; i++) {
                                x      = statePrL[i] * tmpSliceL[i];
                                spread = SHi + Grid;

                                Grid  *= rJumpHi2;
                                dGrid *= rJumpHi2;

                                F  += spread * x;
                                dF += dGrid * x;
                        }
                    }  // for j0
        
                    break;

                case 2:     // nIRDim = 0, mBSDim = 2;
                    for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                    {
                        // Find switch point index
                        Zidx  = xSwitch + j0 * (pJump00 + pJump10);
                        Mid   = (int) ceil(-Zidx / pJump11) - 1;
                        Mid   = Min(Max(Mid, mBottom2[t][j0] - 1),
                                        mTop2[t][j0]);

                        tmpSliceL = tmpSlice + NodeOffset(2, j0, 0, t);
                        statePrL  = statePr  + NodeOffset(2, j0, 0, t);


                        // LEFT part of distribution
                        Zidx += (pJump11) * mBottom2[t][j0];
                        Grid  = MLo * exp (mBSQLo * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQLo * VolBbq * Zidx);

                        for (i = mBottom2[t][j0]; i <= Mid; i++) {
                                x      = statePrL[i] * tmpSliceL[i];
                                spread = SLo + Grid;
                        
                                Grid  *= rJumpLo2;
                                dGrid *= rJumpLo2;

                                F  += spread * x;
                                dF += dGrid * x;
                        }


                        // RIGHT part of distribution
                        Zidx  = xSwitch + j0 * (pJump00 + pJump10) 
                                        + (Mid + 1) * pJump11;
                        Grid  = MHi * exp (mBSQHi * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQHi * VolBbq * Zidx);

                        for (i = Mid + 1; i <= mTop2[t][j0]; i++) {
                                x      = statePrL[i] * tmpSliceL[i];
                                spread = SHi + Grid;

                                Grid  *= rJumpHi2;
                                dGrid *= rJumpHi2;

                                F  += spread * x;
                                dF += dGrid * x;
                        }
                    }  // for j0
        
                    break;
        
                default:
                    throw KFailure("%s: invalid factor dimensions."
                               " nDim=%d, nBSDim=%d.\n", routine, nDim, mBSDim);
                }

                break;

            case 3:

                rJumpLo5 = exp (mBSQLo * VolBbq * pJump22);
                rJumpHi5 = exp (mBSQHi * VolBbq * pJump22); 

                switch (mBSDim) {
                case 1:     // nIRDim = 2, nBSDim = 1;

                    for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                    for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
                    {
                        Zidx  = xSwitch
                              + j0 * (pJump20)
                              + j1 * (pJump21);
                        Mid   = (int) ceil(-Zidx / pJump22) - 1;
                        Mid   = Min(Max(Mid, mBottom3[t][j0][j1] - 1),
                                        mTop3[t][j0][j1]);

                        tmpSliceL = tmpSlice + NodeOffset(3, j0, j1, t);
                        statePrL  = statePr  + NodeOffset(3, j0, j1, t);

                        // LEFT part of distribution
                        Zidx += (pJump22) * mBottom3[t][j0][j1];
                        Grid  = MLo * exp (mBSQLo * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQLo * VolBbq * Zidx);

                        for (i = mBottom3[t][j0][j1]; i <= Mid; i++) {
                                x      = statePrL[i] * tmpSliceL[i];
                                spread = SLo + Grid;
                        
                                Grid  *= rJumpLo5;
                                dGrid *= rJumpLo5;

                                F  += spread * x;
                                dF += dGrid * x;

                        }

                        // RIGHT part of distribution
                        Zidx = xSwitch
                             + j0 * (pJump20)
                             + j1 * (pJump21)
                             + (Mid+1) * (pJump22);
                        Grid  = MHi * exp (mBSQHi * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQHi * VolBbq * Zidx);

                        for (i = Mid + 1; i <= mTop3[t][j0][j1]; i++) {
                                x      = statePrL[i] * tmpSliceL[i];
                                spread = SHi + Grid;
                        
                                Grid  *= rJumpHi5;
                                dGrid *= rJumpHi5;

                                F  += spread * x;
                                dF += dGrid * x;
                        }
                    }

                    break;

                case 2:     // nIRDim = 1, nBSDim = 2;

                    for (j0 = mBottom1[t];     j0 <= mTop1[t];     j0++)
                    for (j1 = mBottom2[t][j0]; j1 <= mTop2[t][j0]; j1++)
                    {
                        Zidx  = xSwitch
                              + j0 * (pJump10 + pJump20)
                              + j1 * (pJump11 + pJump21);
                        Mid   = (int) ceil(-Zidx / pJump22) - 1;
                        Mid   = Min(Max(Mid, mBottom3[t][j0][j1] - 1),
                                        mTop3[t][j0][j1]);

                        tmpSliceL = tmpSlice + NodeOffset(3, j0, j1, t);
                        statePrL  = statePr  + NodeOffset(3, j0, j1, t);

                        // LEFT part of distribution
                        Zidx += (pJump22) * mBottom3[t][j0][j1];
                        Grid  = MLo * exp (mBSQLo * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQLo * VolBbq * Zidx);

                        for (i = mBottom3[t][j0][j1]; i <= Mid; i++) {
                                x      = statePrL[i] * tmpSliceL[i];
                                spread = SLo + Grid;
                        
                                Grid  *= rJumpLo5;
                                dGrid *= rJumpLo5;

                                F  += spread * x;
                                dF += dGrid * x;

                        }

                        // RIGHT part of distribution
                        Zidx = xSwitch
                             + j0 * (pJump10 + pJump20)
                             + j1 * (pJump11 + pJump21)
                             + (Mid+1) * (pJump22);
                        Grid  = MHi * exp (mBSQHi * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQHi * VolBbq * Zidx);

                        for (i = Mid + 1; i <= mTop3[t][j0][j1]; i++) {
                                x      = statePrL[i] * tmpSliceL[i];
                                spread = SHi + Grid;
                        
                                Grid  *= rJumpHi5;
                                dGrid *= rJumpHi5;

                                F  += spread * x;
                                dF += dGrid * x;
                        }
                    }

                    break;

                case 3:     // nIRDim = 0, nBSDim = 3;
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

                        tmpSliceL = tmpSlice + NodeOffset(3, j0, j1, t);
                        statePrL  = statePr  + NodeOffset(3, j0, j1, t);

                        // LEFT part of distribution
                        Zidx += (pJump22) * mBottom3[t][j0][j1];
                        Grid  = MLo * exp (mBSQLo * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQLo * VolBbq * Zidx);

                        for (i = mBottom3[t][j0][j1]; i <= Mid; i++) {
                                x      = statePrL[i] * tmpSliceL[i];
                                spread = SLo + Grid;
                        
                                Grid  *= rJumpLo5;
                                dGrid *= rJumpLo5;

                                F  += spread * x;
                                dF += dGrid * x;
                        }

                        // RIGHT part of distribution
                        Zidx = xSwitch
                             + j0 * (pJump00 + pJump10 + pJump20)
                             + j1 * (pJump11 + pJump21)
                             + (Mid+1) * (pJump22);
                        Grid  = MHi * exp (mBSQHi * VolBbq * Zidx);
                        dGrid = DQ  * exp (mBSQHi * VolBbq * Zidx);

                        for (i = Mid + 1; i <= mTop3[t][j0][j1]; i++) {
                                x      = statePrL[i] * tmpSliceL[i];
                                spread = SHi + Grid;
                        
                                Grid  *= rJumpHi5;
                                dGrid *= rJumpHi5;

                                F  += spread * x;
                                dF += dGrid * x;
                            }
                        }   // for j0, j1

                    break;

                default:
                    throw KFailure("%s: invalid factor dimensions."
                               " nDim=%d, nBSDim=%d.\n", routine, nDim, mBSDim);
                }

                break;    // nDim = 3

            default:
                throw KFailure("%s: N/A for nDim=%d.\n", routine, nDim);

            }    // end of switch nDim


            // Add constant terms
            //
            F -= zeroPrice * fwdBasisRate;


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
                *Zt = xSwitch - baseSh;
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
                dppLog << "Center offset: " <<  *Zt - Zt0 << endl;
                dppLog << endl;
        }


        // Free tmp memory
        sliceDelete(tmpSliceIR);
        sliceDelete(tmpSlice);

        statePrL  = NULL;
        tmpSliceL = NULL;

    }
    catch (KFailure) {
        // Free tmp memory
        sliceDelete(tmpSliceIR);
        sliceDelete(tmpSlice);

        statePrL  = NULL;
        tmpSliceL = NULL;


        throw KFailure("%s: failed.\n", routine);
    }

}



