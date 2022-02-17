/************************************************************************
 * Module:  PenGuin
 * File:
 * Function:    
 * Author:      D. Liu
 ************************************************************************/
#include "vtlkio2idx.h"

#include "kutilios.h"
#include "kmrntree.h"



extern "C" {
#include "drltime.h"
#include "drlmath.h"
};


//===============================================================
//
//===============================================================


static SharedPointer<KVPKnockIO> DownCastToKnockIO
    (const SharedPointer<KVPKnockIO2Idx> ins)
{
    SharedPointer<KVPKnockIO> insDownCast; 
    SharedPointerDownCastTo(ins, insDownCast);

    return insDownCast;
}


//---------------------------------------------------------------

KVPToolKnockIO2Idx::KVPToolKnockIO2Idx(SharedPointer<KVPKnockIO2Idx> ins, 
                                       KVTree &vt)
    : KVPToolKnockIO(DownCastToKnockIO(ins),vt)
{
static  char    routine[] = "KVPToolKnockIO2Idx::KVPToolKnockIO2Idx";

    int   idx, n;
    TDate today = vt.TPToday();

try {
    
    // Trigger rate index curve name 
    mRateIdx2Name = ins->mRateIndex2[0]->Rate().CurveName(); 

    // Only include events with observation date >= today 
    ins->ValidEvents(today);

    mKnockIO2Idx = ins;

    // Insert the second rate index into the tree
    n = (mKnockIO2Idx->mRateIndex2).size();

    for (idx=0; idx<n; idx++) 
        vt.Insert(*(mKnockIO2Idx->mRateIndex2[idx]));
}
catch (KFailure) {

    throw KFailure("%s: failed on knock IO schedule:\n", routine);
}
}

//---------------------------------------------------------------

KVPToolKnockIO2Idx::~KVPToolKnockIO2Idx()
{
    if(debugLevel)
        dppLog << GetName() << ": deleted." << endl;
}


//---------------------------------------------------------------


void
KVPToolKnockIO2Idx::Update()
{
static  char    routine[] = "KVPToolKnockIO2Idx::Update";

    int idx, n = mKnockIO2Idx->mObsDates.size();
    TDate   curDate = mVTree->TPDateCurrent(),
            settleDate, obsDate, modelObsDate;
    int idxD;
    int depCurveIdx, rateCurveIdx, rateCurveIdx2;

    bool    isKnockOut = false;

    double  idxRate, idxStep, idxRate2, idxStep2;
    double  lowerStep, upperStep, rangeInStep; // Smoothing steps 1st idx
    double  lowerStep2, upperStep2;            // Smoothing steps 2nd idx

    KTSlice     *rateTS    = NULL;   // Knock IO rate 1 index slice
    KTSlice     *rate2TS   = NULL;   // Knock IO rate 2 index slice
    KTSlice     *koIndTS   = NULL;   // Knock IO indicator slice (smooth)
    KTSlice     *koUSIndTS = NULL;   // Knock IO indicator slice (unsmooth)
    KTSliceNode tsNode;         // Slice node

    double      underlyingValue = 0e0;

    double  xt0,        // Probability
            xt1,        // expected ko time
            xt2;        // expected ko time^2 

    const   String& discCurveName = mKnockIO2Idx->mDiscZcName;

    String  curveName    = GetCurveName(); // max dim of underlying slice
    String  rateIndexCV  = mRateIdxName; 
    String  rateIndexCV2 = mRateIdx2Name;


    try {

    if (!NeedsUpdate()) return;

    double  exTime    = mVTree->TPTimeCurrent();
    double  exTimeSqr = exTime * exTime;


    // 
    // The slice dimensions of mValue and mExer should be
    // the maximum of rate slice and underlying slice
    //
    depCurveIdx   = mVTree->GetCurveIdx(curveName);
    rateCurveIdx  = mVTree->GetCurveIdx(mRateIdxName);
    rateCurveIdx2 = mVTree->GetCurveIdx(mRateIdx2Name);

    if (rateCurveIdx > depCurveIdx)
    {
        if (rateCurveIdx2 > rateCurveIdx)
            curveName = mRateIdx2Name;
        else
            curveName = mRateIdxName;
    }
    else
    {
        if (rateCurveIdx2 > depCurveIdx)
            curveName = mRateIdx2Name;
    }


    // Set Knock out flag
    //
    if (mKnockIO2Idx->mIOType == CRX_KNOCK_OUT)
        isKnockOut = true;

    
    //------------------------------------------------------
    //
    // (1) Perform dev on the value (both smoothed and unsmoothed)
    // and on the exercise bank.
    //
    //------------------------------------------------------
    if (mValue != NULL) {
        mValue->Dev(discCurveName);
        if (mRunStat) {
            mXT0->Ev();
            mXT1->Ev();
            mXT2->Ev();
        }
    }


    if (mValueUS != NULL) {
        mValueUS->Dev(discCurveName);
    }

    
    // Dev the exercise value
    //
    for(KMap(TDate,KTSlice*)::iterator itEx=mExer.begin();
         itEx != mExer.end();
         ++itEx) {
        ((*itEx).second)->Dev(discCurveName);
    }


    // Dev the rebate. Apply only to knock-out
    //
    if (isKnockOut)
    {
        for(KMap(TDate,KTSlice*)::iterator itRebate=mRebate.begin();
                itRebate != mRebate.end(); ++itRebate) {
        ((*itRebate).second)->Dev(discCurveName);
        }


        //
        // Ensure mRebate slice exists
        //
        if (mRebateValue != NULL) 
        {
            mRebateValue->Dev(discCurveName);
        }

        if (mRebateValueUS != NULL) {
            mRebateValueUS->Dev(discCurveName);
        }
    }
    

    //------------------------------------------------------
    //
    // (2) Add settlements
    //
    //------------------------------------------------------
    for (idx=0; idx<n; idx++) {
        settleDate = mKnockIO2Idx->mSettleDates[idx];

        if (settleDate == curDate) {

        if(debugLevel)
            dppLog << GetName() << format(": Processing settlement %d "
                  " (settlement %s)",
                  idx,
                  DrlTDatePrint(NULL, mKnockIO2Idx->mSettleDates[idx]))
                  << endl;

        //
        // Create new settlement value and exercise slice
        //
        mExer[settleDate] = new KTSlice(*mVTree, "Exercise", curveName);
        (*mExer[settleDate]) = 0e0;

        // Calculate value of the underlying
        // by summimg all dependencies.
        //
        for (idxD=0; idxD<this->NumDep(); idxD++) {
            if (!(Dep(idxD)->IsType("KVPToolRate")))
            (*mExer[settleDate]) += Dep(idxD)->GetValue();
        }


        //
        // Create new rebate value and slice
        //
        if (isKnockOut)
        {
            mRebate[settleDate] 
                = new KTSlice(*mVTree, "Rebate", curveName);
            (*mRebate[settleDate]) = mKnockIO2Idx->mRebates[idx];
        }


        if (debugLevel == DEBUG_LEVEL_KIO)
        {
                    mVTree->TSliceSpecialOper(*mExer[settleDate],
                                          "TSPRT_MINMAX",
                                          (ostream*) &dppLog);
                    dppLog << endl;

            if (isKnockOut)
            {
                dppLog << routine << ": Rebate" << endl;
                        mVTree->TSliceSpecialOper(*mRebate[settleDate],
                                          "TSPRT_MINMAX",
                                          (ostream*) &dppLog);
                        dppLog << endl;
            }
        }

        }   // if settleDate
    }   // for idx


    //------------------------------------------------------
    //
    // (3) Perform knock-in exercises.
    //
    //------------------------------------------------------

    for (idx=0; idx<n; idx++) {
        obsDate      = mKnockIO2Idx->mObsDates[idx];
        modelObsDate = mModelObsDates[idx];

        if (modelObsDate == curDate) {

        if(debugLevel)
            dppLog << GetName() << format(": Processing observation "
              "%d on date %s (settlement %s)", 
              idx,
              DrlTDatePrint(NULL, curDate), 
              DrlTDatePrint(NULL, mKnockIO2Idx->mSettleDates[idx]))
              << endl;

        //
        // Ensure value slice exists
        //
        if (mValue == NULL) {
            mValue = new KTSlice(*mVTree, "mValue", curveName);
            *mValue = 0e0;
            if (mRunStat) {
                mXT0 = new KTSlice(*mVTree, "mXT0", curveName);
                *mXT0 = 0e0;
                mXT1 = new KTSlice(*mVTree, "mXT1", curveName);
                *mXT1 = 0e0;
                mXT2 = new KTSlice(*mVTree, "mXT2", curveName);
                *mXT2 = 0e0;
            }
        }

        if (mKnockIO2Idx->mSmooth == DOUBLE_SMOOTH && mValueUS == NULL) {
            mValueUS = 
                new KTSlice(*mVTree, "mValueUS", curveName);
            *mValueUS = 0e0;
        }

        if (isKnockOut)
        {
            if (mRebateValue == NULL) 
            {
            mRebateValue 
                = new KTSlice(*mVTree, "mRebateValue", curveName);
            *mRebateValue = 0e0;
            }

            if (mRebateValueUS == NULL && 
                mKnockIO2Idx->mSmooth == DOUBLE_SMOOTH)
            {
            mRebateValueUS = 
                new KTSlice(*mVTree, "mRebateValueUS", curveName);
            *mRebateValueUS = 0e0;
            }
        }


        if (debugLevel == DEBUG_LEVEL_KIO)
        {
            dppLog << "Knock IO value slice BEFORE exercise: " 
                   << endl;
            dppLog << *mValue << endl;

            if (mKnockIO2Idx->mSmooth == DOUBLE_SMOOTH) {
                dppLog << *mValueUS << endl;
            }

            if (isKnockOut)
                dppLog << *mRebateValue << endl;
        }


        //
        // Retrieve settlement value time slice
        //
        settleDate = mKnockIO2Idx->mSettleDates[idx];

        KMap(TDate,KTSlice*)::iterator itEx = mExer.find(settleDate);

        if (itEx == mExer.end()) {
            throw KFailure("%s: Current observation date is %s. "
                       "Can't find exercise slice on "
                       "settlement date %s.\n",
                       routine, 
                       GtoFormatDate(curDate),
                       GtoFormatDate(settleDate));
        }

        KTSlice& exerTS = *((*itEx).second);


        //
        // Retrieve rebate value time slice
        //
        KTSlice *rebateTS = NULL;
        KMap(TDate,KTSlice*)::iterator itRebate 
                    = mRebate.find(settleDate);

        if (isKnockOut)
        {
            if (itRebate == mRebate.end()) {
                throw KFailure("%s: Current observation date is "
                       "%s.  Can't find rebate slice on "
                       "settlement date %s.\n",
                       routine, 
                       GtoFormatDate(curDate),
                       GtoFormatDate(settleDate));
                }
            else
                rebateTS = (*itRebate).second;
        }


        //
        // Perform exercise
        //

        //
        // Create new observation rate slice
        //
        rateIndexCV = (mKnockIO2Idx->mRateIndex[idx]) ->Rate().CurveName(); 
        rateIndexCV2= (mKnockIO2Idx->mRateIndex2[idx])->Rate().CurveName(); 

        rateTS    = new KTSlice(*mVTree, "Observ", rateIndexCV);
        rate2TS   = new KTSlice(*mVTree, "Observ", rateIndexCV2);
        koIndTS   = new KTSlice(*mVTree, "Observ", curveName);
        koUSIndTS = new KTSlice(*mVTree, "Observ", curveName);

        mVTree->Get(*rateTS,  *(mKnockIO2Idx->mRateIndex[idx]));
        mVTree->Get(*rate2TS, *(mKnockIO2Idx->mRateIndex2[idx]));
        
            
        
        if (debugLevel == DEBUG_LEVEL_KIO)
        {
            dppLog << "Knock IO rate index slice" << endl;
            dppLog << *rateTS  << endl;
            dppLog << *rate2TS << endl;

            dppLog << "Underlying slice" << endl;
            dppLog << exerTS << endl;
        }



        // Create smooth and unsmooth KO indicator   
        if(mKnockIO2Idx->mIOWindow == CRX_KNOCK_IN)
        {
            if(mKnockIO2Idx->mIOWindow2 == CRX_KNOCK_IN)
            {
              for(tsNode.begin(*mValue); !tsNode.end(); ++tsNode)
              {
                idxRate  = (*rateTS)[tsNode];
                idxStep  = rateTS->GetNodeStepMax(tsNode);

                idxRate2 = (*rate2TS)[tsNode];
                idxStep2 = rate2TS->GetNodeStepMax(tsNode);

                // Generate unsmoothed KO indicator 
                if(idxRate  > mKnockIO2Idx->mBarrierLo[idx] *(1.-BARRIER_TOL) && 
                   idxRate  < mKnockIO2Idx->mBarrierHi[idx] *(1.+BARRIER_TOL) &&
                   idxRate2 > mKnockIO2Idx->mBarrierLo2[idx]*(1.-BARRIER_TOL) &&
                   idxRate2 < mKnockIO2Idx->mBarrierHi2[idx]*(1.+BARRIER_TOL))
                    rangeInStep = 1.;
                else
                    rangeInStep = 0.;

                (*koUSIndTS)[tsNode] = rangeInStep;     
                    

                // Generate smooth KO indicator 

                // First condition 
                lowerStep = DrlSmoothStepFcn(
                                idxRate - mKnockIO2Idx->mBarrierLo[idx],
                                idxStep);
                upperStep = DrlSmoothStepFcn(
                                idxRate - mKnockIO2Idx->mBarrierHi[idx],
                                idxStep);

                // Second condition 
                lowerStep2 = DrlSmoothStepFcn(
                                idxRate2 - mKnockIO2Idx->mBarrierLo2[idx],
                                idxStep2);
                upperStep2 = DrlSmoothStepFcn(
                                idxRate2 - mKnockIO2Idx->mBarrierHi2[idx],
                                idxStep2);

                rangeInStep = (1. - upperStep)  * lowerStep *
                              (1. - upperStep2) * lowerStep2;

                (*koIndTS)[tsNode] = rangeInStep;
              }
            }
            else if(mKnockIO2Idx->mIOWindow2 == CRX_KNOCK_OUT)
            {
              for(tsNode.begin(*mValue); !tsNode.end(); ++tsNode)
              {
                idxRate  = (*rateTS)[tsNode];
                idxStep  = rateTS->GetNodeStepMax(tsNode);

                idxRate2 = (*rate2TS)[tsNode];
                idxStep2 = rate2TS->GetNodeStepMax(tsNode);

                // Generate unsmoothed KO indicator 
                if(idxRate  > mKnockIO2Idx->mBarrierLo[idx] *(1.-BARRIER_TOL) && 
                   idxRate  < mKnockIO2Idx->mBarrierHi[idx] *(1.+BARRIER_TOL) &&
                  (idxRate2 > mKnockIO2Idx->mBarrierHi2[idx]*(1.-BARRIER_TOL) ||
                   idxRate2 < mKnockIO2Idx->mBarrierLo2[idx]*(1.+BARRIER_TOL)))
                    rangeInStep = 1.;
                else
                    rangeInStep = 0.;

                (*koUSIndTS)[tsNode] = rangeInStep;     
                    

                // Generate smooth KO indicator 

                // First condition 
                lowerStep = DrlSmoothStepFcn(
                                idxRate - mKnockIO2Idx->mBarrierLo[idx],
                                idxStep);
                upperStep = DrlSmoothStepFcn(
                                idxRate - mKnockIO2Idx->mBarrierHi[idx],
                                idxStep);

                // Second condition 
                lowerStep2 = DrlSmoothStepFcn(
                                idxRate2 - mKnockIO2Idx->mBarrierLo2[idx],
                                idxStep2);
                upperStep2 = DrlSmoothStepFcn(
                                idxRate2 - mKnockIO2Idx->mBarrierHi2[idx],
                                idxStep2);

                rangeInStep = (1. - upperStep)  * lowerStep *
                              (1. - (1. - upperStep2) * lowerStep2);

                (*koIndTS)[tsNode] = rangeInStep;
              }
            }
            else
            {
                throw KFailure("%s: invalid knock window type (%d).\n",
                               routine, mKnockIO2Idx->mIOWindow2);
            }
        }
        else if(mKnockIO2Idx->mIOWindow == CRX_KNOCK_OUT)
        {
            if(mKnockIO2Idx->mIOWindow2 == CRX_KNOCK_IN)
            {
              for(tsNode.begin(*mValue); !tsNode.end(); ++tsNode)
              {
                idxRate  = (*rateTS)[tsNode];
                idxStep  = rateTS->GetNodeStepMax(tsNode);

                idxRate2 = (*rate2TS)[tsNode];
                idxStep2 = rate2TS->GetNodeStepMax(tsNode);

                // Generate unsmoothed KO indicator 
                if((idxRate > mKnockIO2Idx->mBarrierHi[idx] *(1.-BARRIER_TOL) ||
                   idxRate  < mKnockIO2Idx->mBarrierLo[idx] *(1.+BARRIER_TOL))&&
                   idxRate2 > mKnockIO2Idx->mBarrierLo2[idx]*(1.-BARRIER_TOL) &&
                   idxRate2 < mKnockIO2Idx->mBarrierHi2[idx]*(1.+BARRIER_TOL))
                    rangeInStep = 1.;
                else
                    rangeInStep = 0.;

                (*koUSIndTS)[tsNode] = rangeInStep;     
                    

                // Generate smooth KO indicator 

                // First condition 
                lowerStep = DrlSmoothStepFcn(
                                idxRate - mKnockIO2Idx->mBarrierLo[idx],
                                idxStep);
                upperStep = DrlSmoothStepFcn(
                                idxRate - mKnockIO2Idx->mBarrierHi[idx],
                                idxStep);

                // Second condition 
                lowerStep2 = DrlSmoothStepFcn(
                                idxRate2 - mKnockIO2Idx->mBarrierLo2[idx],
                                idxStep2);
                upperStep2 = DrlSmoothStepFcn(
                                idxRate2 - mKnockIO2Idx->mBarrierHi2[idx],
                                idxStep2);

                rangeInStep = (1. - (1. - upperStep)  * lowerStep) *
                              (1. - upperStep2) * lowerStep2;

                (*koIndTS)[tsNode] = rangeInStep;
              }
            }
            else if(mKnockIO2Idx->mIOWindow2 == CRX_KNOCK_OUT)
            {
              for(tsNode.begin(*mValue); !tsNode.end(); ++tsNode)
              {
                idxRate  = (*rateTS)[tsNode];
                idxStep  = rateTS->GetNodeStepMax(tsNode);

                idxRate2 = (*rate2TS)[tsNode];
                idxStep2 = rate2TS->GetNodeStepMax(tsNode);

                // Generate unsmoothed KO indicator 
                if((idxRate > mKnockIO2Idx->mBarrierHi[idx] *(1.-BARRIER_TOL) || 
                   idxRate  < mKnockIO2Idx->mBarrierLo[idx] *(1.+BARRIER_TOL))&&
                  (idxRate2 > mKnockIO2Idx->mBarrierHi2[idx]*(1.-BARRIER_TOL) ||
                   idxRate2 < mKnockIO2Idx->mBarrierLo2[idx]*(1.+BARRIER_TOL)))
                    rangeInStep = 1.;
                else
                    rangeInStep = 0.;

                (*koUSIndTS)[tsNode] = rangeInStep;     
                    

                // Generate smooth KO indicator 

                // First condition 
                lowerStep = DrlSmoothStepFcn(
                                idxRate - mKnockIO2Idx->mBarrierLo[idx],
                                idxStep);
                upperStep = DrlSmoothStepFcn(
                                idxRate - mKnockIO2Idx->mBarrierHi[idx],
                                idxStep);

                // Second condition 
                lowerStep2 = DrlSmoothStepFcn(
                                idxRate2 - mKnockIO2Idx->mBarrierLo2[idx],
                                idxStep2);
                upperStep2 = DrlSmoothStepFcn(
                                idxRate2 - mKnockIO2Idx->mBarrierHi2[idx],
                                idxStep2);

                rangeInStep = (1. - (1. - upperStep)  * lowerStep) *
                              (1. - (1. - upperStep2) * lowerStep2);

                (*koIndTS)[tsNode] = rangeInStep;
              }
            }
            else
            {
                throw KFailure("%s: invalid knock window type (%d).\n",
                               routine, mKnockIO2Idx->mIOWindow2);
            }
        }
        else
        {
            throw KFailure("%s: invalid knock window type (%d).\n",
                           routine, mKnockIO->mIOWindow);
        }

        
        
        // Perform Knock IO 
        if (mKnockIO2Idx->mSmooth == NO_SMOOTH)
        {
          for(tsNode.begin(*mValue); !tsNode.end(); ++tsNode)
          {
            // Use unsmoothed indicator
            rangeInStep = (*koUSIndTS)[tsNode];         

            (*mValue)[tsNode] += 
                (exerTS[tsNode] - (*mValue)[tsNode]) * rangeInStep;

            if (isKnockOut)
            {
                (*mRebateValue)[tsNode] += rangeInStep *
                    ((*rebateTS)[tsNode] - (*mRebateValue)[tsNode]);
            }

            if (mRunStat) 
            {
                (*mXT0)[tsNode] += (1e0       - (*mXT0)[tsNode]) * rangeInStep;
                (*mXT1)[tsNode] += (exTime    - (*mXT1)[tsNode]) * rangeInStep;
                (*mXT2)[tsNode] += (exTimeSqr - (*mXT2)[tsNode]) * rangeInStep;
            }

            if (debugLevel == DEBUG_LEVEL_KIO)
            {
                dppLog << "KIO: underlying value: " 
                       << exerTS[tsNode];
                dppLog << "KIO after exercise: "
                       << (*mValue)[tsNode] << endl;
            }
          }
        }
        else if (mKnockIO2Idx->mSmooth == SINGLE_SMOOTH) 
        {
          for(tsNode.begin(*mValue); !tsNode.end(); ++tsNode)
          {
            // Use smooth indicator
            rangeInStep = (*koIndTS)[tsNode];         

            (*mValue)[tsNode] += 
                (exerTS[tsNode] - (*mValue)[tsNode]) * rangeInStep;

            if (isKnockOut)
            {
                (*mRebateValue)[tsNode] += rangeInStep *
                    ((*rebateTS)[tsNode] - (*mRebateValue)[tsNode]);
            }

            if (mRunStat) 
            {
                (*mXT0)[tsNode] += (1e0       - (*mXT0)[tsNode]) * rangeInStep;
                (*mXT1)[tsNode] += (exTime    - (*mXT1)[tsNode]) * rangeInStep;
                (*mXT2)[tsNode] += (exTimeSqr - (*mXT2)[tsNode]) * rangeInStep;
            }

            if (debugLevel == DEBUG_LEVEL_KIO)
            {
                dppLog << "KIO: underlying value: " 
                       << exerTS[tsNode];
                dppLog << "KIO after exercise: "
                       << (*mValue)[tsNode] << endl;
            }
          }
        }
        else if (mKnockIO2Idx->mSmooth == DOUBLE_SMOOTH) 
        {
          for(tsNode.begin(*mValue); !tsNode.end(); ++tsNode)
          {
            // Use smooth indicator
            rangeInStep = (*koIndTS)[tsNode];         

            // Use unsmoothed value inside the smoothing window
            if (rangeInStep > 0.0 && rangeInStep < 1.0) 
            {
                (*mValue)[tsNode] =   (*mValueUS)[tsNode] +
                    (exerTS[tsNode] - (*mValueUS)[tsNode]) * rangeInStep;

                if (isKnockOut)
                {
                    (*mRebateValue)[tsNode] =  (*mRebateValueUS)[tsNode]  +
                        ((*rebateTS)[tsNode] - (*mRebateValueUS)[tsNode]) *
                        rangeInStep;
                }
            }
            else // Use smooth value
            {
                (*mValue)[tsNode] =   (*mValue)[tsNode] +
                    (exerTS[tsNode] - (*mValue)[tsNode]) * rangeInStep;

                if (isKnockOut)
                {
                    (*mRebateValue)[tsNode] =  (*mRebateValue)[tsNode]  +
                        ((*rebateTS)[tsNode] - (*mRebateValue)[tsNode]) *
                        rangeInStep;
                }
            }

            if (mRunStat) 
            {
                (*mXT0)[tsNode] += (1e0       - (*mXT0)[tsNode]) * rangeInStep;
                (*mXT1)[tsNode] += (exTime    - (*mXT1)[tsNode]) * rangeInStep;
                (*mXT2)[tsNode] += (exTimeSqr - (*mXT2)[tsNode]) * rangeInStep;
            }

            // Update unsmoothed slice
            rangeInStep = (*koUSIndTS)[tsNode];

            (*mValueUS)[tsNode] = (*mValueUS)[tsNode] +
                (exerTS[tsNode] - (*mValueUS)[tsNode]) * rangeInStep;

            if (isKnockOut)
            {
                (*mRebateValueUS)[tsNode] = (*mRebateValueUS)[tsNode]  +
                    ((*rebateTS)[tsNode]  - (*mRebateValueUS)[tsNode]) *
                    rangeInStep;
            }

            if (debugLevel == DEBUG_LEVEL_KIO)
            {
                dppLog << "KIO: underlying value: " 
                       << exerTS[tsNode];
                dppLog << "KIO after exercise: "
                       << (*mValue)[tsNode] << endl;
            }
          }            
        }

        if (debugLevel == DEBUG_LEVEL_KIO)
        {
            dppLog << "Knock IO value slice AFTER exercise: " 
                   << endl;
            dppLog << *mValue << endl;

            if (isKnockOut)
                dppLog << *mRebateValue << endl;
        }


        //
        // Erase exerTs and rebateTS
        //
        delete (*itEx).second;
        mExer.erase(itEx);

        if (isKnockOut)
        {
            delete (*itRebate).second;
            mRebate.erase(itRebate);
            
            rebateTS = NULL;
        }


        }   // if modelObsDate

    }     // for idx


    //
    // Last TP: store the result.
    //
    if (mVTree->TPIdxCurrent() == 0) {
        if (mValue == NULL) {
            throw KFailure("%s: no observation event "
                "after today date.\n", routine);
        }


        // Store PV
        mResults.insert(
            String_Double("PV", mValue->GetCenter())
        );

        // Calculate PV value of the underlying
        // by summimg all dependencies.
        //
        for (idxD=0; idxD<this->NumDep(); idxD++) {
            if (!(Dep(idxD)->IsType("KVPToolRate")))
            underlyingValue += Dep(idxD)->GetResults()["PV"];
        }
        

        mResults.insert(
            String_Double("UNDER", underlyingValue)
        );

        // Expected time to hit
        //
        if (mRunStat) {
            xt0 = mXT0->GetCenter();
            xt1 = mXT1->GetCenter();
            xt2 = mXT2->GetCenter();

            double  pe,         // prob exercise
                te,         // expected exer time
                se,         // var exer time
                tZ;         // last exer time
            TDate   lastDate = mKnockIO->mObsDates.back();

            tZ = (lastDate - mVTree->TPToday()) / 365e0;

            pe = xt0;
            te = xt1 + (1-pe)*tZ;

            mResults["XTPROB"] = pe;
            if (isKnockOut)
                mResults["XTPROB"] = 1e0 - pe;
            mResults["XTEXP"]  = xt1;
            mResults["XTFUG"]  = te;

            se = (xt2 + (1-pe)*tZ*tZ) - te * te;
            if (se < 0e0) {
                mResults["XTSDV"] = sqrt(-se);
            } else {
                mResults["XTSDV"] = sqrt(se);
            }
        } else {
            mResults["XTPROB"] = -9999.99;
            mResults["XTEXP"]  = -9999.99;
            mResults["XTFUG"]  = -9999.99;
            mResults["XTSDV"]  = -9999.99;
        }

        // Correct for knock-out using parity, plus the rebate.
        //
        if(isKnockOut)
        {
            double  koPV   = mResults["UNDER"]
                       - mResults["PV"]
                       + mRebateValue->GetCenter();

            mResults["PV"] = koPV;
        }

    }


    //
    // Done updating, set flag
    //
    UpdateDone();


    delete rateTS;
    delete rate2TS;
    delete koIndTS;
    delete koUSIndTS;

    }
    catch (KFailure) {

    delete rateTS;
    delete rate2TS;
    delete koIndTS;
    delete koUSIndTS;

    // Don't need to be here anymore.
    // Default destructor should take care of it
    //
    for (idx=0; idx<n; idx++) {
        //
        // Retrieve settlement value time slice
        //
        settleDate = mKnockIO->mSettleDates[idx];

        KMap(TDate,KTSlice*)::iterator itEx = mExer.find(settleDate);

        if (itEx != mExer.end()) {
            delete (*itEx).second;
            mExer.erase(itEx);
        }

        KMap(TDate,KTSlice*)::iterator itRebate 
                    = mRebate.find(settleDate);

        if (itRebate != mRebate.end()) {
            delete (*itRebate).second;
            mRebate.erase(itRebate);
        }

    }

    throw KFailure("%s: failed.\n", routine);
    }
}

