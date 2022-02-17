/************************************************************************
 * Module:        PenGuin
 * File:
 * Function:    
 * Author:  Robert Mattingly, 
 * Modified by David Liu for risky
 ************************************************************************/
#include "vtlfleg3Idx.h"
#include "kvtspar.h"        // Slice parser
#include "kmktcrv.h"

#include "kutilios.h"


extern        KMap(String,double)        globVpConstTable;        // Constant tables


extern "C" {
#include "drltime.h"
};

#define DEBUG_CURVE
#undef        DEBUG_CURVE

//===============================================================
//
//===============================================================


//---------------------------------------------------------------

KVPToolFloatLeg3Idx::KVPToolFloatLeg3Idx
(
 SharedPointer<KVPFloatLeg3Idx> ins,        // (I) instrument to value
 KVTree &vt                             // (I) virtual tree
 ) : KVPToolFloatLeg(vt)
{
  static        char        routine[] = "KVPToolFloatLeg3Idx::KVPToolFloatLeg3Idx";
  int        idx, n;
  int         cIdx, cIdxMax;
    
  TDate        todayDate = vt.TPToday();
  TDate        valueDate2, valueDate3;
  TDate        resetDate, modelResetDate;
    
  KRateReset        *stubRate = NULL;
    
  String        curveName, curveName2, curveName3;
    
  const        String&        discCurveName = ins->mDiscZcName;

  try 
    {
  
      SharedPointer<KVPFloatLeg> vpflt;
      SharedPointerConvertTo(ins, vpflt);
      KVPToolFloatLeg::CreateVPToolFloatLeg(vpflt, vt);        
      
      // Check valid
      
      mFloatLeg3Idx = ins;
      
      n = mFloatLeg3Idx->mResetDates.size();
      
      mFloatLeg3Idx->CheckValid();
      
      // enough space to store the model reset dates
      mModelResetDates2.resize(n);
      mModelResetDates3.resize(n);
      
      //
      // The value date of the leg will be determined from 
      // the curve index of
      // 1. First reset period if reset in the future.
      // OR
      // 2. Stub period if reset in the past.
      //
      KVector(TDate)::iterator firstResetDate = 
        min_element(
                    mFloatLeg3Idx->mResetDates.begin(),
                    mFloatLeg3Idx->mResetDates.end()
                    );
      
      valueDate2 = -1;
      valueDate3 = -1;

      //
      // Add 2nd set of index rates in the tree
      //
      
      // The largest rate curve index from the 
      // 1st set of index rates
      //
      cIdxMax = mVTree->GetCurveIdx(mCurveName);
      for (idx=0; idx<n; idx++) {
        
        resetDate = mFloatLeg3Idx->mResetDates[idx];
        
        if (mFloatLeg3Idx->mRateResets2[idx]->Rate().IsFloating())
          curveName2 = mFloatLeg3Idx->mRateResets2[idx]->Rate().CurveName();
        else
          curveName2 = mFloatLeg3Idx->mDiscZcName;
        
        // Value date is curve dependent
        valueDate2 = vt.GetValueDate(vt.GetCurveIdx(curveName2));
        
        // Value date is the min of the two indices.
        // The value date of the 1st index leg is known already.
        // 
        // First reset could be in the past or future.
        // * if in the future, then this would be THE value date.
        // * if in the past, then it would be replaced by the first
        //   stub period treated below.
        //
        if (*firstResetDate == mFloatLeg->mResetDates[idx] &&
            *firstResetDate > todayDate)
          mValueDate = MIN(mValueDate, valueDate2);
        
        //
        // Identify largest geometry to allocate slice
        // The geometry is given by the curve index associated to
        // a name by the virtual tree.
        // THE CONVENTION IS THAT INCREASING INDICES REPRESENT
        // INCREASING GEOMETRIES.
        //
        cIdx = mVTree->GetCurveIdx(curveName2);
        if (cIdx > cIdxMax)
          cIdxMax = cIdx;
        
        if (mFloatLeg3Idx->mAccEndDates[idx] <= valueDate2) // obsolete events
          continue;
        else         // insert valid events
          {
            //
            // Add events to tree.
            // Rate reset: we collect the model reset date
            // anywhere between reset date and accural end date.
            // This will allow us to compute rate anywhere
            // between reset date and accural end date, 
            // and ensure that rate stub or American type
            // of exercise can be evaluated before reaching
            // the reset date.
            // 
            // (i.e. where the rate is known in the tree)
            // Pay date: we need the zero from pay to reset.
            // Create (and store) rate and zero resets
            // for each coupon payment
            //
            //
            if (cIdx < KV_BASIS ||
                resetDate >= todayDate)        
              {
                mModelResetDates2[idx] 
                  = vt.Insert(*(mFloatLeg3Idx->mRateResets2[idx]),
                              mFloatLeg3Idx->mAccEndDates[idx]);
              }
            
            // 
            // Determine if front stub exists, and insert
            // appropriate stub rate to the zero bank for calculation
            // later on.
            //
            if (resetDate < todayDate &&
                mFloatLeg3Idx->mRateResets2[idx]->Rate().IsFloating())
              // Stub period
              {
                // The value date of the leg is determined from
                // the first valid period.
                mValueDate = MIN(valueDate2, mValueDate);
                
                switch ((long)(mFloatLeg3Idx->mStubConv)) {
                case GTO_STUB_NONE:
                case GTO_STUB_SIMPLE:
                case GTO_STUB_BOND:
                  if (cIdx >= KV_BASIS)
                    {
                      mModelResetDates2[idx] 
                        = vt.Insert(*(mFloatLeg3Idx->mRateResets2[idx]),
                                    mFloatLeg3Idx->mAccEndDates[idx]);
                    }
                  // Else get the rate reset value from reset bank
                  
                  break;
                  
                case KV_STUB_SPOT_SIMPLE:
                  // Same index rate spot reset
                  stubRate = new KRateReset
                    (
                     todayDate,
                     valueDate2,
                     mFloatLeg3Idx->mRateResets2[idx]->Rate()
                     );
                  
                  modelResetDate 
                    = vt.Insert(*stubRate,
                                mFloatLeg3Idx->mAccEndDates[idx]);
                  delete stubRate;
                  
                  if (cIdx >= KV_BASIS)
                    mModelResetDates[idx] = modelResetDate;
                  
                  break;
                  
                case KV_STUB_PAR:        // use short rate between VDate and AE
                  
                  // Simple stub rate between value date
                  // and accEndDate, with DCC and spotOffset
                  // the same as RateReset2.
                  //
                  stubRate = new KRateReset
                    (
                     todayDate,
                     valueDate2,
                     mFloatLeg3Idx->mAccEndDates[idx],
                     mFloatLeg3Idx->mRateResets2[idx]->Rate()
                     );
                  
                  
                  modelResetDate 
                    = vt.Insert(*stubRate);
                  
                  delete stubRate;
                  
                  if (cIdx >= KV_BASIS)
                    mModelResetDates[idx] = modelResetDate;
                  
                  break;
                  
                default:
                  throw KFailure("%s: Invalid stub convention (%d).\n",
                                 routine, (long)(mFloatLeg->mStubConv));
                } // switch case
              }

            // In order for the multi-index rate payment formula to work
            // the modelResetDate must be the same
            if (mModelResetDates2[idx] !=
                mModelResetDates[idx])
              throw KFailure("%s: %d-th model reset date (%s) for "
                             "rate 1 != model reset date (%s) for rate "
                             "2.\n",
                             routine,
                             idx,
                             GtoFormatDate(mModelResetDates[idx]),
                             GtoFormatDate(mModelResetDates2[idx]));
            
          }
      }

      //
      // Add 3rd set of index rates in the tree
      //
      
      // The largest rate curve index from the 
      // 1st set of index rates
      //
      for (idx=0; idx<n; idx++) {
        
        resetDate = mFloatLeg3Idx->mResetDates[idx];
        
        if (mFloatLeg3Idx->mRateResets3[idx]->Rate().IsFloating())
          curveName3 = mFloatLeg3Idx->mRateResets3[idx]->Rate().CurveName();
        else
          curveName3 = mFloatLeg3Idx->mDiscZcName;
        
        // Value date is curve dependent
        valueDate3 = vt.GetValueDate(vt.GetCurveIdx(curveName3));
        
        // Value date is the min of the two indices.
        // The value date of the 1st index leg is known already.
        // 
        // First reset could be in the past or future.
        // * if in the future, then this would be THE value date.
        // * if in the past, then it would be replaced by the first
        //   stub period treated below.
        //
        if (*firstResetDate == mFloatLeg->mResetDates[idx] &&
            *firstResetDate > todayDate)
          mValueDate = MIN(mValueDate, valueDate3);
        
        //
        // Identify largest geometry to allocate slice
        // The geometry is given by the curve index associated to
        // a name by the virtual tree.
        // THE CONVENTION IS THAT INCREASING INDICES REPRESENT
        // INCREASING GEOMETRIES.
        //
        cIdx = mVTree->GetCurveIdx(curveName3);
        if (cIdx > cIdxMax)
          cIdxMax = cIdx;
        
        
        if (mFloatLeg3Idx->mAccEndDates[idx] <= valueDate3) // obsolete events
          continue;
        else         // insert valid events
          {
            //
            // Add events to tree.
            // Rate reset: we collect the model reset date
            // anywhere between reset date and accural end date.
            // This will allow us to compute rate anywhere
            // between reset date and accural end date, 
            // and ensure that rate stub or American type
            // of exercise can be evaluated before reaching
            // the reset date.
            // 
            // (i.e. where the rate is known in the tree)
            // Pay date: we need the zero from pay to reset.
            // Create (and store) rate and zero resets
            // for each coupon payment
            //
            //
            if (cIdx < KV_BASIS ||
                resetDate >= todayDate)        
              {
                mModelResetDates3[idx] 
                  = vt.Insert(*(mFloatLeg3Idx->mRateResets3[idx]),
                              mFloatLeg3Idx->mAccEndDates[idx]);
              }
            
            // 
            // Determine if front stub exists, and insert
            // appropriate stub rate to the zero bank for calculation
            // later on.
            //
            if (resetDate < todayDate &&
                mFloatLeg3Idx->mRateResets3[idx]->Rate().IsFloating())
              // Stub period
              {
                // The value date of the leg is determined from
                // the first valid period.
                mValueDate = MIN(valueDate3, mValueDate);
                
                switch ((long)(mFloatLeg3Idx->mStubConv)) {
                case GTO_STUB_NONE:
                case GTO_STUB_SIMPLE:
                case GTO_STUB_BOND:
                  if (cIdx >= KV_BASIS)
                    {
                      mModelResetDates3[idx] 
                        = vt.Insert(*(mFloatLeg3Idx->mRateResets3[idx]),
                                    mFloatLeg3Idx->mAccEndDates[idx]);
                    }
                  // Else get the rate reset value from reset bank
                  
                  break;
                  
                case KV_STUB_SPOT_SIMPLE:
                  // Same index rate spot reset
                  stubRate = new KRateReset
                    (
                     todayDate,
                     valueDate3,
                     mFloatLeg3Idx->mRateResets3[idx]->Rate()
                     );
                  
                  modelResetDate 
                    = vt.Insert(*stubRate,
                                mFloatLeg3Idx->mAccEndDates[idx]);
                  delete stubRate;
                  
                  if (cIdx >= KV_BASIS)
                    mModelResetDates[idx] = modelResetDate;
                  
                  break;
                  
                case KV_STUB_PAR:        // use short rate between VDate and AE
                  
                  // Simple stub rate between value date
                  // and accEndDate, with DCC and spotOffset
                  // the same as RateReset2.
                  //
                  stubRate = new KRateReset
                    (
                     todayDate,
                     valueDate3,
                     mFloatLeg3Idx->mAccEndDates[idx],
                     mFloatLeg3Idx->mRateResets3[idx]->Rate()
                     );
                  
                  
                  modelResetDate 
                    = vt.Insert(*stubRate);
                  
                  delete stubRate;
                  
                  if (cIdx >= KV_BASIS)
                    mModelResetDates[idx] = modelResetDate;
                  
                  break;
                  
                default:
                  throw KFailure("%s: Invalid stub convention (%d).\n",
                                 routine, (long)(mFloatLeg->mStubConv));
                } // switch case
              }

            // In order for the multi-index rate payment formula to work
            // the modelResetDate must be the same
            if (mModelResetDates3[idx] !=
                mModelResetDates[idx])
              throw KFailure("%s: %d-th model reset date (%s) for "
                             "rate 1 != model reset date (%s) for rate "
                             "3.\n",
                             routine,
                             idx,
                             GtoFormatDate(mModelResetDates[idx]),
                             GtoFormatDate(mModelResetDates3[idx]));
            
          }
      }
      
      //
      // Insert one extra period at the end of the floating leg to 
      // ensure that spread drift can be interpolated for any back 
      // stub in basis rate. 
      //
      KVector(TDate)::iterator lastDate = 
        max_element
        (
         mFloatLeg3Idx->mAccEndDates.begin(),
         mFloatLeg3Idx->mAccEndDates.end()
         );
      
      KRateReset extraRateReset2
        (
         (*lastDate),
         (*lastDate) + 
         mFloatLeg3Idx->mRateResets2[n-1]->Rate().SpotOffset(),
         mFloatLeg3Idx->mRateResets2[n-1]->Rate()
         );
      
      vt.Insert(extraRateReset2);
      
      KRateReset extraRateReset3
        (
         (*lastDate),
         (*lastDate) + 
         mFloatLeg3Idx->mRateResets3[n-1]->Rate().SpotOffset(),
         mFloatLeg3Idx->mRateResets3[n-1]->Rate()
         );
      
      vt.Insert(extraRateReset3);

      //
      // Curve name of the leg is the one with largest slice geometry
      // among the 3 sets of index rates
      // need to include discount curve in the case of risky
      //
      curveName = mFloatLeg3Idx->mDiscZcName;
      cIdx = mVTree->GetCurveIdx(curveName);
      if (cIdx > cIdxMax)
          cIdxMax = cIdx;

      //
      mCurveName = mVTree->GetCurveName(cIdxMax);
      
    }
  catch (KFailure) 
    {
      if((long)(mFloatLeg3Idx->mStubConv) == KV_STUB_SPOT_SIMPLE) 
        delete stubRate;
    
      throw KFailure("%s: failed on instrument:\n", routine);
    }
  
 
}





//---------------------------------------------------------------

KVPToolFloatLeg3Idx::~KVPToolFloatLeg3Idx()
{
}


//---------------------------------------------------------------
// Updates floating cash flow price by discounting 
// the price from the next timestep to this timestep, and then adding 
// in any payments that occur on this timestep.

void
KVPToolFloatLeg3Idx::Update()
{
  static        char        routine[] = "KVPToolFloatLeg3Idx::Update";
  int        idx, n = mFloatLeg->mResetDates.size();
  TDate        currDate = mVTree->TPDateCurrent(),
    resetDate,
    modelResetDate,
    accStartDate,
    accEndDate,
    payDate;
  double        notional,
    dayCnt;
        
  const        String& discCurveName =  mFloatLeg->mDiscZcName;
  const        String& curveName = GetCurveName();

#ifdef  DEBUG_CURVE
  dppLog << routine << ": Curve : " << curveName << endl;
  dppLog << routine << ": Discount: " << discCurveName << endl;
#endif

  KTSlice        *newPymtTs = NULL;
  KTSlice        *rateArrayTs = NULL;

  KRateReset         *rateReset1 = NULL;
  KRateReset      *rateReset2 = NULL;
  KRateReset      *rateReset3 = NULL;

  try {
    if (!NeedsUpdate()) return;

    //
    // Perform dev on the value and on the exercise bank.
    // We use the disc curve name of the product
    //
    if (mUnValue != NULL) {

      mUnValue->Dev(discCurveName);

      // Dev pending resets
      for(KMap(int,KTSlice*)::iterator p=mPResets.begin();
          p != mPResets.end();
          ++p) {
        if(debugLevel)
          dppLog << GetName() <<
            format(": Updating (%d pending resets).",
                   (int) mPResets.size())
                 << endl;

        mPResets[(*p).first]->Dev(discCurveName);
      }
    }


    //
    // Process
    //
    for (idx=0; idx<n; idx++) 
      {
        if (mFloatLeg->mAccEndDates[idx] <= mValueDate) // obsolete events
          continue;
        else         // Valid events
          {
            resetDate      = mFloatLeg3Idx->mResetDates[idx];
            modelResetDate = mModelResetDates[idx];
            accStartDate   = mFloatLeg3Idx->mAccStartDates[idx];
            accEndDate     = mFloatLeg3Idx->mAccEndDates[idx];
            payDate        = mFloatLeg3Idx->mPayDates[idx];
            notional       = mFloatLeg3Idx->mNotionals[idx];
            rateReset1     = mFloatLeg3Idx->mRateResets[idx];
            rateReset2     = mFloatLeg3Idx->mRateResets2[idx];
            rateReset3     = mFloatLeg3Idx->mRateResets3[idx];


            //
            // Ensure slices allocated
            //
                
            if ((mUnValue == NULL) && (payDate >= currDate)) 
              {
                mUnValue = new KTSlice(*mVTree, "mUnValue", curveName);
                mClValue = new KTSlice(*mVTree, "mClValue", curveName);
                mTmpTs   = new KTSlice(*mVTree, "mTmpTs",   curveName);
                mFRateTs = new KTSlice(*mVTree, "mFRateTs", curveName);
                    
                (*mUnValue) = 0e0;
                (*mClValue) = 0e0;
                (*mTmpTs)   = 0e0;
              }



            //
            // if Reset date: the floating rate can be computed.
            // Because we may need to stub (bond), we do not add the 
            // cash flow to the total value, but put it
            // in the pending payments
            // WARNING: the cash-flows in pending payments
            //          are NOT multiplied by either dcc or notional.
            //

            if (currDate == modelResetDate) 
              {
                newPymtTs = new KTSlice
                  (
                   *mVTree,
                   format(
                          "mPReset[%d]-%s-%s", idx,
                          DrlTDatePrint(NULL, resetDate),
                          DrlTDatePrint(NULL, modelResetDate)),
                   curveName);

                //
                // Allocate 3 rate slices
                //
                rateArrayTs = KTSliceNewVector
                  (
                   3,
                   *mVTree,
                   curveName,
                   format("rate reset slices at %s-%s",
                          DrlTDatePrint(NULL, resetDate),
                          DrlTDatePrint(NULL, modelResetDate))
                   );

                if (debugLevel)
                  dppLog << GetName() 
                         << format(": Proc Res  %d "
                                   "(R:%s,MR:%s:P:%s,AS:%s,AE:%s)",
                                   idx,
                                   DrlTDatePrint(NULL, resetDate),
                                   DrlTDatePrint(NULL, modelResetDate),
                                   DrlTDatePrint(NULL, payDate),
                                   DrlTDatePrint(NULL, accStartDate),
                                   DrlTDatePrint(NULL, accEndDate))
                         << endl;
                    
                //
                // Get rate reset
                // Since we have modified the rate reset for
                // different stub conventions in the initialiation,
                // this should always work, even in the first period.
                //
                mVTree->Get(rateArrayTs[0], *rateReset1);
                mVTree->Get(rateArrayTs[1], *rateReset2);
                mVTree->Get(rateArrayTs[2], *rateReset3);
                    
                //
                // If yacction (payment formula) applicable
                //
                if (mFloatLeg->mFormulas[idx].size() != 0) 
                  {
                    KTSliceParEval
                      (
                       *mVTree,
                       3,
                       rateArrayTs,
                       mFloatLeg->mFormulas[idx].c_str(),
                       globVpConstTable,        // Constant tables
                       newPymtTs
                       );
                    (*mTmpTs) = (*newPymtTs);
                  }
                else        // use the 1st rate index for payment calc
                  (*mTmpTs) = rateArrayTs[0];
                    

                delete [] rateArrayTs;
                rateArrayTs = NULL;

                // if paid in future, get zero
                //
                if (payDate != modelResetDate) 
                  {
                    mVTree->Get(*newPymtTs, *mIRZeroResets[idx]);
                    *newPymtTs *= *mTmpTs;
                  } 
                else 
                  {
                    *newPymtTs  = *mTmpTs;
                  }
                
                
                //
                // WARNING: DO NOT Mutliply by dcc or notional at this stage !
                // because we may have to do stub calculations.
                //


                //
                // Add to pending resets bank
                //
                mPResets[idx] = newPymtTs;
                newPymtTs = NULL;        // so we now don't free it !
              }

            //
            // if reset and accrual in future, the paymnt is
            // fully accrued. Look for it in the bank, if it is
            // there, remove it by adding it to the fully accrued pymts
            // TS value.
            // WARNING: need to multiply by dcc and notional at this stage
            //
            if (currDate <= mFullAccDates[idx])
            {
                KMap(int,KTSlice*)::iterator p = mPResets.find(idx);
                if (p != mPResets.end()) {

                        //
                        // Day count fraction
                        //
                        dayCnt = DayCountFraction(
                                                accStartDate,
                                                accEndDate,
                                                mFloatLeg->mDayCc);
                        //
                        // Get pymt (rate reset)
                        //
                        KTSlice *pymtTs = mPResets[(*p).first];


                        // Default accrual interest adjustment
                        // !!! This is a hack. Should be implemented properly 
                        // !!! as part of risky zero bank management, and used 
                        // !!! by CMCDS index calculation as well !!!
                        if (mDefAccrual)
                        {
                                //
                                // Default prob between [accSt, accEnd]
                                // Negative because Z-1
                                //
                                mVTree->Get(*mTmpTs,   *mDefZeroAccEnds[idx]);
                                mVTree->Get(*mFRateTs, *mDefZeroAccStarts[idx]);

                                *mTmpTs /= (*mFRateTs);
                                *mTmpTs -= 1.0;

                                //
                                // Average default accrual is 1/2. 
                                // "-" sign to change Z-1 to 1-Z 
                                //
                                *mTmpTs *= (-0.5);


                                // Add no default payment
                                *mTmpTs += 1.0;

                                // Effective DCF
                                *mTmpTs *= dayCnt;

                                // Multiply by effective DCF
                                *(pymtTs) *= (*mTmpTs);
                        }
                        else
                                *(pymtTs) *= dayCnt;

                        //
                        // Multiply by notional
                        //
                        *(pymtTs) *= notional;


                        //
                        // Add to the fully accrued values
                        //
                        *mUnValue += *(pymtTs);


                        //
                        // Erase form bank
                        //
                        delete mPResets[(*p).first];
                        mPResets.erase(p);
                                }
              }
                }
      }
        
        
    //
    // Done updating, set flag
    //
    UpdateDone();
        
        
  }
  catch (KFailure) 
    {
      delete newPymtTs;
      delete [] rateArrayTs;
        
      rateReset1 = NULL;
      rateReset2 = NULL;
      rateReset3 = NULL;
        
      throw KFailure("%s: failed.\n", routine);
    }
}


//---------------------------------------------------------------
// For two set of index rates, we only handled the case here when 
// both are floating rates.  If one of index rates is fixed,
// the payment formula can always be broken down as floating leg with
// single set of index rates with adjusted payment formula.
// 

KTSlice&
KVPToolFloatLeg3Idx::GetValue()
{
  static        char        routine[] = "KVPToolFloatLeg3Idx::GetValue";
  int        idx, n = mFloatLeg->mResetDates.size();
  TDate        currDate = mVTree->TPDateCurrent(),
    todayDate = mVTree->TPToday(),
    currEffDate,
    resetDate,
    resetEffDate,
    modelResetDate,
    accStartDate,
    accEndDate,
    payDate;
        
  KRateReset         *rateReset1         = NULL;
  KRateReset        *rateReset2         = NULL;
  KRateReset        *rateReset3         = NULL;
  KRateReset        *stubRateReset1 = NULL;
  KRateReset        *stubRateReset2 = NULL;
  KRateReset        *stubRateReset3 = NULL;

  double notional, dayCnt;

  const String&        discCurveName = mFloatLeg->mDiscZcName;

  KTSlice *newPymtTs   = NULL;
  KTSlice *zTs         = NULL;
  KTSlice *rateArrayTs = NULL;

  const        String& curveName = GetCurveName();
  int          curveIdx = mVTree->GetCurveIdx(curveName);

  try {

    // Before any reset date, mValue = NULL,
    // set the default value to 0
    if (mDefValue == NULL) 
      {
        mDefValue = new KTSlice(*mVTree, "mDefValue", GetCurveName());
        ASSERT_OR_THROW(mDefValue != NULL);
        
        *mDefValue = 0e0;
      }

    //
    // Check exists, if not return default value 0
    //
    if (mUnValue == NULL) 
      {
        mDefValue->SetTpIdx(mVTree->TPIdxCurrent());
        return (*mDefValue);
      }

    //
    // The total clean value will consist of 
    // (1) the sum of all already fully accrued (i.e. future) cashflows
    // (2) the cash flows that require either a stub or an estimation
    //     if unknown in the tree (not past the modelResetDate)
    //

    //
    // (1) Add all fully accrued coupons
    //
    (*mClValue) = (*mUnValue);
    (*mTmpTs)   = 0e0;

    //
    // Allocate 3 rate slices
    //
    rateArrayTs = KTSliceNewVector
      (
       3,
       *mVTree,
       curveName,
       format("Stub rate slices at %s\n.",
              DrlTDatePrint(NULL, currDate))
       );
    
    //
    // (2) 
    //
    for (idx=0; idx<n; idx++) {

      if (mFloatLeg->mAccEndDates[idx] <= mValueDate) // obsolete events
        {
          continue;
        }
      else         // Valid events
        {
          resetDate      = mFloatLeg->mResetDates[idx];
          resetEffDate   = mFloatLeg->mResetEffDates[idx];
          modelResetDate = mModelResetDates[idx];
          accStartDate   = mFloatLeg->mAccStartDates[idx];
          accEndDate     = mFloatLeg->mAccEndDates[idx];
          payDate        = mFloatLeg->mPayDates[idx];
          notional       = mFloatLeg->mNotionals[idx];
          rateReset1     = mFloatLeg->mRateResets[idx];
          rateReset2     = mFloatLeg3Idx->mRateResets2[idx];
          rateReset3     = mFloatLeg3Idx->mRateResets3[idx];

          //
          //
          //                     |<---------- f1 ----------->|
          // ------------C-2-----C-----C+2-------------------M-----M+2-----
          //            Reset                  |<---------- f2 ----------->| 
          //
          // Current date in GetValue usually refers to the settlement date 
          // of a swap, equivalent of value date.  The reset date 
          // should be prior to the settlement date by SpotOffset,
          // which is not precisely known due to lack of holiday adj info.
          // 
          // The reset date is only used in constructing stub rate
          // in SIMPLE/PAR conventions.  If we use the correct reset C-2
          // to construct the stub rate f1 (i.e. reset at C-2, and effective
          // at C), then what actually being computed in the tree is 
          // rate f2 (reset at C and effective at C+2, adjusted for forward
          // ratio f1/f2) because of past reset.
          //
          // This can be avoided if we set reset date = current date, ignoring
          // the spotOffset.  This way the rate being computed IS f1 
          // using all the correct zeros.
          //
          // For NONE/BOND stubs, the approx of rate with past reset
          // is unavoidable where the rate reset is clearly defined by
          // the product.  The current date is entered only in the dcf 
          // calculation.  There is no need to define the current reset 
          // date.
          //

          if (currDate <= accStartDate  &&
              currDate <= modelResetDate) {
            //
            // ----C<=AS--------------------
            // ----C<=R---------------------
            // If accrued and reset occured at or after currDate,
            // payment has already been included, as it should be.
            //
            //

            continue;

          } else if (currDate >= accEndDate  &&
                     currDate >  modelResetDate) {
            //
            // ----------AE<=C---------------
            // -----R--------C---------------
            // If accrued and reset in the past, payment not included,
            // as it should not be. Note that since no accrual takes
            // place on the accrueEndDate, currDate = accrueEndDate
            // means accrual occurs in the past.
            //


            continue;

          } else if (currDate < accEndDate  &&
                     currDate > modelResetDate) { 
            //
            // ----R--------C----------AE--------
            //
            // NEED STUB RATE !
            // If currDate falls somewhere in or before the accrual
            // period, and the reset is in the past, the payment
            // has not been included, but it needs to be.
            // This could be one of the following cases with respect to
            // the location of accrual start date:

            // Get zero for pay date and multiply by dcc*rate
            mVTree->Get(*mTmpTs, 
                        KZeroReset(discCurveName, 
                                   currDate,
                                   payDate));
            
            //        
            // Only both index rates being floating
            // are handled here explicitly.
            //
            if ((rateReset1->Rate().IsFloating()) &&
                (rateReset2->Rate().IsFloating()) &&
                (rateReset3->Rate().IsFloating()))
              {
                //
                // --- FLOATING RATE ---
                //
                // The modelResetDate is strictly in the past
                // so that the floating cash-flow has not been
                // added at all. We need to either compute a 
                // simple stub rate (SIMPLE/PAR) or estimate 
                // the past rate reset (BOND/NONE) (the estimation 
                // method is handled by the virtual tree) and 
                // compute the proper stubs.
                // NONE:        + RATE*dcf(S,E)*Z
                // SIMPLE:        + RATE*dcf(C,E)*Z
                // SPOT_SIMPLE:        + spotRateReset*dcf(C,E)*Z
                // BOND:        + RATE*(dcf(S,E)*Z - dcf(S,C)*1)
                // PAR:         + stubRate*dcf(C,E)*Z, 
                // (for pure IR MM rate,  = (Z0 - Zn)*Z/Zn, 
                //                where Z0 and Zn are computed from rate curve,
                //                and Z is computed from discount curve.)

                //
                // 1) ----R---AS---C----------AE------- 
                //
                // 2) ----R-C-AS--------------AE-------
                //
                // 3) ---AS-----R----C--------AE------- 
                //
                // Case 1) Simple stub.  For PAR/SPOT_SIMPLE, 
                // construct the stub rate with 
                //        resetDate    = currDate 
                //      resetEffDate = MAX(currDate, valueDate).
                //
                // Case 2) Fully accured, No stub. But rate unknown. 
                // Constuct the rate with
                //        resetDate    = currDate
                //        resetEffDate = MAX(currDate, mResetEffDate);
                // This would also ensure that GetValue works when
                // currDate < valueDate.
                //
                // Case 3) Basis with 1/3 delay shift.  However,
                // for basis with SIMPLE/PAR stub, since the stub
                // rate will be simple rate reset at spot (current date), 
                // only 0 delay shift is allowed.  Otherwise,
                // treated the same as case 1).
                //
                // Combined, we have
                //         resetDate    = currDate
                //        resetEffDate = MAX(currDate, valueDate, mResetEffDate) 
                //

                // Compute the stub multiplier

        //
        // !!! Default accrual is ignored for stub !!!!
        //
        //

                switch ((long)(mFloatLeg->mStubConv)) {
                case GTO_STUB_NONE:           // Rate reset at beginning of period
                  dayCnt = DayCountFraction(
                                            accStartDate,
                                            accEndDate,
                                            mFloatLeg->mDayCc);
                  *mTmpTs *= dayCnt;
                        
                  stubRateReset1 = rateReset1;
                  stubRateReset2 = rateReset2;
                  stubRateReset3 = rateReset3;

                  break;

                case GTO_STUB_BOND:  // Rate reset at beginning of period
                  dayCnt = DayCountFraction(
                                            accStartDate,
                                            accEndDate,
                                            mFloatLeg->mDayCc);
                  *mTmpTs *= dayCnt;

                  dayCnt = DayCountFraction(
                                            accStartDate,
                                            MAX(currDate, accStartDate),
                                            mFloatLeg->mDayCc);
                  *mTmpTs -= dayCnt;

                  stubRateReset1 = rateReset1;
                  stubRateReset2 = rateReset2;
                  stubRateReset3 = rateReset3;

                  break;

                case GTO_STUB_SIMPLE:    // Rate reset at beginning of prd
                  dayCnt = DayCountFraction(
                                            MAX(currDate, accStartDate),
                                            accEndDate,
                                            mFloatLeg->mDayCc);
                  *mTmpTs *= dayCnt;

                  stubRateReset1 = rateReset1;
                  stubRateReset2 = rateReset2;
                  stubRateReset3 = rateReset3;

                  break;

                case KV_STUB_SPOT_SIMPLE:    // Same index rate spot reset
                  currEffDate = MAX(currDate, 
                                    MAX(mValueDate, resetEffDate));
                  dayCnt = DayCountFraction(
                                            MAX(currDate, accStartDate),
                                            accEndDate,
                                            mFloatLeg->mDayCc);
                  *mTmpTs *= dayCnt;

                  stubRateReset1 = new KRateReset(currDate, 
                                                  currEffDate,
                                                  rateReset1->Rate());

                  stubRateReset2 = new KRateReset(currDate, 
                                                  currEffDate,
                                                  rateReset2->Rate());

                  stubRateReset3 = new KRateReset(currDate, 
                                                  currEffDate,
                                                  rateReset3->Rate());

                  break;

                case KV_STUB_PAR:    // short rate betwn currEffDate and AE
                  currEffDate = MAX(currDate, 
                                    MAX(mValueDate, resetEffDate));

                  // PAR stub requires accEndDate = payDate
                  if (accEndDate != payDate)
                    throw KFailure("%s: accEndDate (%s) != payDate (%s)"
                                   " for par stub.\n",
                                   routine,
                                   GtoFormatDate(accEndDate),
                                   GtoFormatDate(payDate));

                  dayCnt = DayCountFraction(
                                            MAX(currDate, accStartDate),
                                            accEndDate,
                                            mFloatLeg->mDayCc);
                  *mTmpTs *= dayCnt;
                        
                  // Simple rate between [currDate, payDate]
                  // with the same DayCc, SpotOffset, and Spread
                  // as rateReset.
                  //
                  stubRateReset1 = new KRateReset(currDate,
                                                  currEffDate,
                                                  payDate,
                                                  rateReset1->Rate());

                  stubRateReset2 = new KRateReset(currDate,
                                                  currEffDate,
                                                  payDate,
                                                  rateReset2->Rate());

                  stubRateReset3 = new KRateReset(currDate,
                                                  currEffDate,
                                                  payDate,
                                                  rateReset3->Rate());

                  break;

                default:
                  throw KFailure("%s: Error getting reset %d "
                                 " at date %s (stub method not supported).",
                                 routine,
                                 idx, DrlTDatePrint(NULL, currDate));
                }


                
                // Compute the rate (estimated) reset
                mVTree->Get(rateArrayTs[0], *stubRateReset1);
                mVTree->Get(rateArrayTs[1], *stubRateReset2);
                mVTree->Get(rateArrayTs[2], *stubRateReset3);

                if ((long)(mFloatLeg->mStubConv) == KV_STUB_SPOT_SIMPLE ||
                    (long)(mFloatLeg->mStubConv) == KV_STUB_PAR)
                  {
                    delete stubRateReset1;
                    delete stubRateReset2;
                    delete stubRateReset3;
                  }


                //
                // If yacction (payment formula) applicable
                //
                if (mFloatLeg->mFormulas[idx].size() != 0) {
                  newPymtTs = new KTSlice(
                                          *mVTree,
                                          format(
                                                 "mPReset[%d]-%s-%s", idx,
                                                 DrlTDatePrint(NULL, resetDate),
                                                 DrlTDatePrint(NULL, modelResetDate)),
                                          mCurveName);

                  KTSliceParEval(
                                 *mVTree,
                                 3,
                                 rateArrayTs,
                                 mFloatLeg->mFormulas[idx].c_str(),
                                 globVpConstTable,        // Constant tables
                                 newPymtTs);
                  (*mFRateTs) = (*newPymtTs);

                  delete newPymtTs;
                }
                else        // default: use 1st index for payment calc
                  (*mFRateTs) = rateArrayTs[0];

                *mTmpTs *= *mFRateTs;

              }
            else
              throw KFailure("%s: only allowed when both set of "
                             "index rates are floating rates\n."
                             "If one or both rates are fixed, it should be "
                             "possible to decompose it into floating leg "
                             "with single set of index rates with payment "
                             "formula.\n",
                             routine);

            //
            // Multiply by notional
            //
            *(mTmpTs) *= notional;

          } else if (currDate <= modelResetDate &&
                     currDate > accStartDate) {
            // 
            // ------AS-------C<=R-----------
            //
            // Reset is in the future, so we already added in the coupon.
            // But some or all of the accrual period is in the past.
            // 1) NONE: some of the coupon needs to be removed
            // 2) BOND: invalid case for BOND. CheckValid earlier.
            // 3) PAR:  invalid case for reset in arrear.
            // 4) SIMPLE: recompute stub rate and payment from tree.
            //
            // This could happen in the following cases:
            // 1. The reset occurs in arrears.
            // 2. Front stub in basis, where 
            //    currDate = modelReset = todayDate.
            // 3. In basis where reset is at 1/3 of the accrual 
            //    period, and current date is within the first 1/3 
            //    period.
            //
            // For basis with SIMPLE/PAR stub, since the stub
            // rate will be simple rate reset at current date, 
            // only 0 delay shift is allowed.  Therefore case 3
            // actually does not apply. 
            //


            if ((rateReset1->Rate().IsFloating()) &&
                (rateReset2->Rate().IsFloating()) &&
                (rateReset3->Rate().IsFloating()))
              {
                //
                // --- FLOATING RATE ---
                //
                // Retreive the coupon ts
                KMap(int,KTSlice*)::iterator p = mPResets.find(idx);
                if (p == mPResets.end()) {
                  throw KFailure("%s: Error getting reset %d "
                                 " at date %s.",
                                 routine,
                                 idx, DrlTDatePrint(NULL, currDate));
                }


                //
                // Compute the day count adjustment
                //
                switch ((long)(mFloatLeg->mStubConv)) {
                case GTO_STUB_NONE:
                  // Get the full payment
                  dayCnt = DayCountFraction(
                                            accStartDate,
                                            accEndDate,
                                            mFloatLeg->mDayCc);

                  *mTmpTs = *((*p).second);
                  *mTmpTs *= dayCnt;

                  break;

                case GTO_STUB_BOND:
                  // 
                  // Only valid when currDate = modelReset = todayDate.
                  //
                  if (currDate < modelResetDate &&
                      currDate != todayDate)
                    throw KFailure("%s: bond stub is not allowed at "
                                   "currDate (%s) < modelResetDate (%s) and "
                                   "currDate (%s) != todayDate (%s).\n",
                                   routine,
                                   GtoFormatDate(currDate),
                                   GtoFormatDate(modelResetDate),
                                   GtoFormatDate(currDate),
                                   GtoFormatDate(todayDate));
                  else {
                    dayCnt = DayCountFraction(
                                              accStartDate,
                                              accEndDate,
                                              mFloatLeg->mDayCc);

                    *mTmpTs = *((*p).second);
                    *mTmpTs *= dayCnt;

                    // Remove the past accrual 
                    dayCnt = DayCountFraction(
                                              accStartDate,
                                              MIN(currDate, accEndDate),
                                              mFloatLeg->mDayCc);


                    // Compute the rate (estimated or past) reset
                    mVTree->Get(rateArrayTs[0], *rateReset1);
                    mVTree->Get(rateArrayTs[1], *rateReset2);
                    mVTree->Get(rateArrayTs[2], *rateReset3);

                    //
                    // If yacction (payment formula) applicable
                    //
                    if (mFloatLeg->mFormulas[idx].size() != 0) {
                      newPymtTs = new KTSlice
                        (
                         *mVTree,
                         format(
                                "mPReset[%d]-%s-%s", idx,
                                DrlTDatePrint(NULL, resetDate),
                                DrlTDatePrint(NULL, 
                                              modelResetDate)),
                         mCurveName);

                      KTSliceParEval(
                                     *mVTree,
                                     3,
                                     rateArrayTs,
                                     mFloatLeg->mFormulas[idx].c_str(),
                                     globVpConstTable,  // Constant tables
                                     newPymtTs);
                      (*mFRateTs) = (*newPymtTs);

                      delete newPymtTs;
                    }
                    else  // default: use 1st index for payment calc
                      (*mFRateTs) = rateArrayTs[0];

                    (*mFRateTs) *= dayCnt;
                    (*mTmpTs) -= (*mFRateTs);
                  }

                  break;

                case GTO_STUB_SIMPLE:
                  dayCnt = DayCountFraction(
                                            currDate,
                                            MAX(currDate, accEndDate),
                                            mFloatLeg->mDayCc);

                  *mTmpTs = *((*p).second);
                  *mTmpTs *= dayCnt;

                  break;

                case KV_STUB_SPOT_SIMPLE:
                  // Only apply to Basis rate
                  //
                  if (curveIdx < KV_BASIS ||
                      curveIdx > KV_PAR_SPREAD)
                  {
                       throw KFailure("%s: Stub method not supported "
                                "(Only NONE/SIMPLE stubs are allowed "
                                "for reset in arrears)\n.",
                                routine);
                  }
                  else
                  {
                      // Get zero for pay date and multiply by dcc*rate
                      zTs   = new KTSlice(*mVTree, "zTs", discCurveName);
                      mVTree->Get(*zTs,
                                   KZeroReset(discCurveName,
                                              currDate,
                                              payDate));

                      dayCnt = DayCountFraction(
                                            currDate,
                                            MAX(currDate, accEndDate),
                                            mFloatLeg->mDayCc);
                      (*zTs) *= dayCnt;

                      stubRateReset1 = new KRateReset(currDate,
                                                  currDate,
                                                  rateReset1->Rate());

                      stubRateReset2 = new KRateReset(currDate,
                                                  currDate,
                                                  rateReset2->Rate());

                      stubRateReset3 = new KRateReset(currDate,
                                                  currDate,
                                                  rateReset3->Rate());


                      // Compute the spot rate reset
                      mVTree->Get(rateArrayTs[0], *stubRateReset1);
                      mVTree->Get(rateArrayTs[1], *stubRateReset2);
                      mVTree->Get(rateArrayTs[2], *stubRateReset3);

                      delete stubRateReset1;
                      delete stubRateReset2;
                      delete stubRateReset3;

                      // If yacction (payment formula) applicable
                      //
                      if (mFloatLeg->mFormulas[idx].size() != 0) {
                        newPymtTs = new KTSlice
                          (
                           *mVTree,
                           format(
                              "mPReset[%d]-%s-%s", idx,
                              DrlTDatePrint(NULL, resetDate),
                              DrlTDatePrint(NULL,
                                            modelResetDate)),
                                            mCurveName);

                        KTSliceParEval(
                                   *mVTree,
                                   3,
                                   rateArrayTs,
                                   mFloatLeg->mFormulas[idx].c_str(),
                                   globVpConstTable,  // Constant tables
                                   newPymtTs);
                        (*mFRateTs) = (*newPymtTs);

                        delete newPymtTs;
                      }
                      else  // default: use 1st index for payment calc
                        (*mFRateTs) = rateArrayTs[0];

                      *mTmpTs  = (*zTs);
                      *mTmpTs *= (*mFRateTs);
                  }

                  break;

                case KV_STUB_PAR: 
                  // Only apply to Basis rate
                  //
                  if (curveIdx < KV_BASIS ||
                      curveIdx > KV_PAR_SPREAD)
                  {
                       throw KFailure("%s: Stub method not supported "
                                "(Only NONE/SIMPLE stubs are allowed "
                                "for reset in arrears)\n.",
                                routine);
                  }
                  else
                  {
                      // Get zero for pay date and multiply by dcc*rate
                      zTs   = new KTSlice(*mVTree, "zTs", discCurveName);
                      mVTree->Get(*zTs,
                                   KZeroReset(discCurveName,
                                              currDate,
                                              payDate));

                      dayCnt = DayCountFraction(
                                            currDate,
                                            MAX(currDate, accEndDate),
                                            mFloatLeg->mDayCc);
                      (*zTs) *= dayCnt;

                      // Simple stub rate
                      stubRateReset1 = new KRateReset(currDate,
                                                  currDate,
                                                  payDate,
                                                  rateReset1->Rate());

                      stubRateReset2 = new KRateReset(currDate,
                                                  currDate,
                                                  payDate,
                                                  rateReset2->Rate());

                      stubRateReset3 = new KRateReset(currDate,
                                                  currDate,
                                                  payDate,
                                                  rateReset3->Rate());


                      // Compute the spot rate reset
                      mVTree->Get(rateArrayTs[0], *stubRateReset1);
                      mVTree->Get(rateArrayTs[1], *stubRateReset2);
                      mVTree->Get(rateArrayTs[2], *stubRateReset3);

                      delete stubRateReset1;
                      delete stubRateReset2;
                      delete stubRateReset3;


                      // If yacction (payment formula) applicable
                      //
                      if (mFloatLeg->mFormulas[idx].size() != 0) {
                        newPymtTs = new KTSlice
                          (
                           *mVTree,
                           format(
                              "mPReset[%d]-%s-%s", idx,
                              DrlTDatePrint(NULL, resetDate),
                              DrlTDatePrint(NULL,
                                            modelResetDate)),
                           mCurveName
                           );

                        KTSliceParEval(
                                   *mVTree,
                                   3,
                                   rateArrayTs,
                                   mFloatLeg->mFormulas[idx].c_str(),
                                   globVpConstTable,  // Constant tables
                                   newPymtTs);
                        (*mFRateTs) = (*newPymtTs);

                        delete newPymtTs;
                      }
                      else  // default: use 1st index for payment calc
                        (*mFRateTs) = rateArrayTs[0];

                      *mTmpTs  = (*zTs);
                      *mTmpTs *= (*mFRateTs);
                  }
                  break;

                default:
                  throw KFailure("%s: Error getting reset %d "
                                 " at date %s (stub method not supported).",
                                 routine,
                                 idx, DrlTDatePrint(NULL, currDate));
                }  // end of case

              }  

            //
            // Multiply by notional
            //
            *(mTmpTs) *= notional;

            delete zTs;

          } else {
            //
            // Case not handled
            //

            throw KFailure("%s: %s not handled "
                           "Current: %s Reset: %s ModelReset: %s "
                           "AccStart: %s AccEnd: %s\n",
                           routine, GetName(),
                           DrlTDatePrint(NULL, currDate),
                           DrlTDatePrint(NULL, resetDate),
                           DrlTDatePrint(NULL, modelResetDate),
                           DrlTDatePrint(NULL, accStartDate),
                           DrlTDatePrint(NULL, accEndDate));
          }

          //
          // Add payment to clean value
          //
          (*mClValue) += (*mTmpTs);
        }
    }

    delete [] rateArrayTs;
    rateArrayTs = NULL;

    return(*mClValue);

  }
  catch(KFailure) {
    delete zTs;

    if ((long)(mFloatLeg->mStubConv) == KV_STUB_SPOT_SIMPLE ||
        (long)(mFloatLeg->mStubConv) == KV_STUB_PAR)
      {
        delete stubRateReset1;
        delete stubRateReset2;
      }


    delete [] rateArrayTs;
    delete newPymtTs;

    rateReset1 = NULL;
    rateReset2 = NULL;
    rateReset3 = NULL;

    throw KFailure("%s: failed on %s at date %s.\n",
                   routine, GetName(), DrlTDatePrint(NULL, currDate));
  }
}

