/************************************************************************
 * Module:        PenGuin
 * File:
 * Function:    
 * Author:      David Liu
 ************************************************************************/
#include "vtldefprotect.h"
#include "kmktcrv.h"

#include "kutilios.h"

extern "C" {
#include "drltime.h"
#include "drlsort.h"
};

//===============================================================

//---------------------------------------------------------------

KVPToolDefProtect::KVPToolDefProtect(SharedPointer<KVPDefProtect> ins, 
                                     KVTree &vt)
        : KVPToolAtom(vt)
{
static        char        routine[] = "KVPToolDefProtect::KVPToolDefProtect";

        int        idx, n, discIdx;
        TDate      startDate, endDate, settleDate;

    try {

        mDefProtect = ins;

        // Check valid
        mDefProtect->CheckIsValid();

        // Discount curve has to be risky
        //
        discIdx   = mVTree->GetCurveIdx(mDefProtect->GetDiscName());
        if (discIdx < KV_CREDIT_RISKY)
            throw KFailure("%s: discount curve '%s' for protection component "
                           "has to be risky.\n",
                           routine,
                           mDefProtect->mDiscZcName.c_str());

        mValueDate =  vt.GetValueDate(discIdx);

        n = mDefProtect->mStartDates.size();
        mDefZeros.resize(n);

        for (idx=0; idx<n; idx++) {
            startDate  = mDefProtect->mStartDates[idx];
            endDate    = mDefProtect->mEndDates[idx];
            settleDate = mDefProtect->mSettleDates[idx];
            mDefZeros[idx] = NULL;         

            if(endDate > mValueDate)
            {
                mDefZeros[idx] = new KZeroReset(
                           mVTree->GetCurveName(KV_CREDIT_DEFPROB),
                           MAX(mValueDate, startDate),
                           endDate);

                vt.Insert(*mDefZeros[idx]);
                vt.Insert(settleDate);
            }
        }

        mValue    = NULL;
        mDefValue = NULL;
        mTmpTs    = NULL;

    }
    catch (KFailure) {
        throw KFailure("%s: failed on protection initialization:\n", routine);
    }
}

//---------------------------------------------------------------

KVPToolDefProtect::~KVPToolDefProtect()
{
        int idx, n; 

        n = mDefProtect->mStartDates.size();

        delete mValue;
        delete mDefValue;
        delete mTmpTs;

        for (KMap(int,KTSlice*)::iterator itUnder = mUnder.begin();
                        itUnder != mUnder.end(); ++itUnder)
        {
                delete (*itUnder).second;
        }
        mUnder.clear();

        for (idx=0; idx<n; idx++) {
                if (mDefZeros[idx])  delete mDefZeros[idx];
        } 

        if(debugLevel)
                dppLog << GetName() << ": deleted." << endl;
}


//---------------------------------------------------------------


inline const String&
KVPToolDefProtect::GetCurveName()
{
static        char        routine[] = "KVPToolDefProtect::GetCurveName";
        if (mCurveName.empty()) {
          try {
            int        idx, idxMax, cIdxMax, cIdx;
            int        discIdx;

            //
            // Identify largest geometry to allocate slice
            // The geometry is given by the curve index associated to
            // a name by the virtual tree.
            // THE CONVENTION IS THAT INCREASING INDICES REPRESENT
            // INCREASING GEOMETRIES.
            //
            cIdxMax = -100000;
            idxMax = -1;
            for(idx=0; idx<NumDep(); idx++) {
                cIdx = mVTree->GetCurveIdx(Dep(idx)->GetCurveName());

                //
                // In default protection valuation, the underlying
                // should be "riskless" in the current single credit name 
                // framework.
                //
                if (cIdx >= KV_CREDIT_RISKY) 
                        throw KFailure("%s: underlying '%s' must be riskless. "
                              "(curve '%s' is risky).\n",
                              routine,
                              mDefProtect->Dep(idx)->GetName(),
                              Dep(idx)->GetCurveName().c_str());

                if (cIdx > cIdxMax) {
                        idxMax = idx;
                        cIdxMax = cIdx;
                }
            }

            // include discount curve
            discIdx = mVTree->GetCurveIdx(mDefProtect->GetDiscName());
            if (discIdx > cIdxMax)
                cIdxMax = discIdx;

            mCurveName = mVTree->GetCurveName(cIdxMax);


        }
          catch (KFailure) {
            throw KFailure("KVPToolDefProtect::GetCurveName: "
                "failed on protection `%s'.\n", mDefProtect->GetName());
          }
        }
        return mCurveName;
}




//---------------------------------------------------------------


void
KVPToolDefProtect::Update()
{
static        char        routine[] = "KVPToolDefProtect::Update";

        int        idx, n = mDefProtect->mStartDates.size();
        TDate      curDate = mVTree->TPDateCurrent();
        TDate      startDate, endDate, settleDate;
        int        idxD;


const        String&   discCurveName = mDefProtect->mDiscZcName;
const        String&   curveName = GetCurveName();
const        String&   discIRCurveName = mVTree->GetIRDiscCurveName(discCurveName);


    try {

        if (!NeedsUpdate()) return;


        //
        // Perform dev on the value and on the exercise bank.
        //
        if (mValue != NULL) {
                mValue->Dev(discCurveName);
        }


        // In default protection valuation, the underlying
        // should be "riskless", because the default probability
        // is taken into account explicitly.
        //
        for(KMap(int,KTSlice*)::iterator p=mUnder.begin();
             p != mUnder.end();
             ++p) {
                mUnder[(*p).first]->Dev(discIRCurveName);
        }



        //
        // Add settlements
        //
        for (idx=0; idx<n; idx++) {
            settleDate = mDefProtect->mSettleDates[idx];

            if (settleDate == curDate) {

                if (debugLevel) 
                        dppLog << GetName() << ": Processing setllement "
                                << idx <<endl;

                //
                // Create new settlement value slice
                //
                mUnder[idx] = new KTSlice(*mVTree, "Settle", curveName);

                if (mDefProtect->GetDefType() == DEF_KNOCKIN)
                {
                    // Add rebate
                    *(mUnder[idx]) = (mDefProtect->mRebates[idx]);

                    //
                    // Calculate value of the underlying
                    // by summimg all dependencies.
                    //
                    for (idxD=0; idxD<this->NumDep(); idxD++) {
                        (*mUnder[idx]) += this->Dep(idxD)->GetValue(); 
                    } 
                }
                else  // DEF_EXPOSURE
                {

                    //
                    // Only one underlying is allowed
                    // positive value only for exposure.
                    //
                    if (this->NumDep() != 1)
                    {
                        throw KFailure("%s: only one underlying is allowed for "
                                       "default exposure calc (got %d). \n",
                                       routine,
                                       this->NumDep());
                    }

                    *(mUnder[idx]) = 0e0;
                    mUnder[idx]->max(this->Dep(0)->GetValue());

                    // Multiply by 1-R 
                    (*mUnder[idx]) *= (1e0 - mDefProtect->mRecovery);
                }

            }
        }


        //
        // Multiply by default probability 
        //
        for (idx=0; idx<n; idx++) {
            if (mValueDate >= mDefProtect->mEndDates[idx])
                continue;

            startDate  = MAX(mValueDate, mDefProtect->mStartDates[idx]);
            endDate    = mDefProtect->mEndDates[idx];
            settleDate = mDefProtect->mSettleDates[idx];

            //
            // If is start of the protection period
            //
            if (startDate == curDate) {

                if (debugLevel) 
                    dppLog << GetName() << format(
                        ": Processing protection %d (settlement %s)",
                        idx, 
                        DrlTDatePrint(NULL, settleDate)) << endl;

                //
                // Ensure value slice exists
                //
                if (mValue == NULL) {
                        mValue = new KTSlice(*mVTree, "mValue", curveName);
                        mTmpTs = new KTSlice(*mVTree, "mTmpTs", curveName);
                        *mValue = 0e0;
                }


                //
                // Retrieve settlement value time slice
                //
                KMap(int,KTSlice*)::iterator p = mUnder.find(idx);

                if (p == mUnder.end()) {
                        throw KFailure("%s: can't find settle date %s.\n",
                                routine, DrlTDatePrint(NULL, settleDate));
                }

                KTSlice& underTs = *(mUnder[(*p).first]);


                //
                // Calculate default probability for the period 
                // [startDate, endDate].
                //
                mVTree->Get(*mTmpTs, *mDefZeros[idx]);
                *mTmpTs -= 1.0;
                *mTmpTs *= (-1.0);

                // Multiply underlying by default prob for the period
                //
                *mTmpTs *= underTs; 

                //
                // Add to the mValue
                //
                *mValue += *(mTmpTs);

                //
                // Erase underTs
                //
                delete mUnder[(*p).first];
                mUnder.erase(p);

            }
        }



        //
        // Done updating, set flag
        //
        UpdateDone();


    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}


//---------------------------------------------------------------

KTSlice&
KVPToolDefProtect::GetValue()
{
static        char        routine[] = "KVPToolDefProtect::GetValue";

        TDate   currDate = mVTree->TPDateCurrent();
        int     idx;

        // Before any protection period occurs, mValue=NULL,
        // set the protection to default value 0.
        //
        if (mDefValue == NULL) {
            mDefValue = new KTSlice(*mVTree, "mDefValue", GetCurveName());
            ASSERT_OR_THROW(mDefValue != NULL);
            *mDefValue = 0e0;
        }

        //
        // Compute stub if needed
        //
        if (currDate >= mDefProtect->mEndDates.back())
        {
            mDefValue->SetTpIdx(mVTree->TPIdxCurrent());        
            return(*mDefValue);

        }
        else if (currDate > MAX(mValueDate, mDefProtect->mStartDates.front()))
        {

            //
            // Check if we need stub - in the middle of protection period
            //
            KVector(TDate)::iterator iterDate = find(
                                    mDefProtect->mStartDates.begin(),
                                    mDefProtect->mStartDates.end(),
                                    currDate);

            if (iterDate != mDefProtect->mStartDates.end())
                return(*mValue);
            else
            {
                *mDefValue = *mValue;

                // find the current protection period
                //
                DrlTDateArrayCeilIdx(&(mDefProtect->mEndDates[0]),
                                     mDefProtect->mEndDates.size(),
                                     currDate,
                                     &idx);
                                  

                // Stubbed default zero
                KZeroReset defZero = KZeroReset(
                                       mVTree->GetCurveName(KV_CREDIT_DEFPROB),
                                       currDate,
                                       mDefProtect->mEndDates[idx]);

                // Compute the stub default probabiltiy 
                mVTree->Get(*mTmpTs, defZero);
                *mTmpTs -= 1.0;
                *mTmpTs *= (-1.0);

                //
                // Retrieve settlement value time slice
                // for the current period
                //
                KMap(int,KTSlice*)::iterator 
                               p = mUnder.find(idx);

                if (p == mUnder.end()) {
                    throw KFailure("%s: can't find settle date %s.\n",
                           routine, 
                           DrlTDatePrint(NULL, mDefProtect->mSettleDates[idx]));
                }

                KTSlice& underTs = *(mUnder[(*p).first]);

                // Multiply underlying by default prob for the period
                //
                *mTmpTs *= underTs; 

                // Add to mValue 
                *mDefValue += (*mTmpTs);

                return(*mDefValue);

            }  // end of stub calculation
        }
        else   // done with the forward starting protection.
            return(*mValue);
}


