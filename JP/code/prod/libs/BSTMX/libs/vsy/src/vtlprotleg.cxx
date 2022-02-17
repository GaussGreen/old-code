/************************************************************************
 * Module:        PenGuin
 * File:
 * Function:    
 * Author:      D. Liu
 *              Amended by Charles Morcom to include a notional schedule.
 ************************************************************************/
#include "vtlprotleg.h"
#include "kvtspar.h"        // Slice parser
#include "kmktcrv.h"

#include "kutilios.h"


extern "C" {
#include "drltime.h"
};


//===============================================================
//
//===============================================================


//---------------------------------------------------------------
KVPToolProtLeg::KVPToolProtLeg(
        SharedPointer<KVPProtLeg> ins,     // (I) inetrument to value
        KVTree &vt)                        // (I) virtual tree
    : KVPToolAtom(vt), mProtLeg(ins), mDefValue(NULL), mValue(NULL),
      mTmpSlice(NULL) 
{
static        char        routine[] = "KVPToolProtLeg::KVPToolProtLeg";

        TDate        valueDate;
        int          discIdx;
        int n, idx;

    try {


        // Check valid
        mProtLeg->CheckValid();
        n = mProtLeg->mNtlAmts.size();

        // Discount curve has to be risky
        //
        discIdx   = mVTree->GetCurveIdx(mProtLeg->mDiscZcName);
        if (discIdx < KV_CREDIT_RISKY)
            throw KFailure("%s: discount curve '%s' for protection leg "
                           "has to be risky.\n",
                           routine,
                           mProtLeg->mDiscZcName.c_str());
        
        valueDate =  vt.GetValueDate(discIdx);

        // Protection curve name
        // Currently only one credit name allowed to avoid ambiguity.
        // 
        if (mProtLeg->mIsDefRecovery)
            mProtCurveName = mVTree->GetCurveName(KV_CREDIT_PROT_DEFRECOV);
        else
            mProtCurveName = mVTree->GetCurveName(KV_CREDIT_PROT_BINRECOV);

        /* Insert all the zero resets for the protection sections */
        // mProtZero = KVector(KZeroReset*)(n, NULL);
        mProtZero.resize(n);

        for (idx=0; idx<n; idx++) {
            mProtZero[idx] = NULL;

            TDate startDate = mProtLeg->mNtlStartDts[idx];
            TDate endDate = mProtLeg->mNtlEndDts[idx];

            // Skip the past protection period
            if (endDate<valueDate)  continue;

            startDate = MAX(startDate, valueDate);
            // create protection zero
            mProtZero[idx] = new KZeroReset(mProtCurveName,
                                   startDate, 
                                   endDate);
            // register it
            vt.Insert(*mProtZero[idx]);
        }


    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}





//---------------------------------------------------------------

KVPToolProtLeg::~KVPToolProtLeg()
{
        int idx;
        if (mDefValue) delete mDefValue;
        if (mValue) delete mValue;
        for (idx=0; idx<mProtZero.size(); idx++) {
            if (mProtZero[idx]) delete mProtZero[idx];
        }

        if(debugLevel)
                dppLog << GetName() << ": deleted." << endl;
}



//---------------------------------------------------------------
// Updates protection leg price by special default discounting 

void
KVPToolProtLeg::Update()
{
static        char        routine[] = "KVPToolProtLeg::Update";

        TDate   currDate = mVTree->TPDateCurrent();
        int n = mProtLeg->mNtlEndDts.size();
        int idx;
        
    try {

        if (!NeedsUpdate()) return;

        // Always perform risky discount of value so far
        if (mValue != NULL)
            mValue->Dev(mProtLeg->GetDiscName());



        for (idx=0; idx<n; idx++) {
            // create value and temp slices when first pass protection end date
            // leave inside loop in case end dates are ever not in idx order!
            if (mValue==NULL && mProtLeg->mNtlEndDts[idx]>=currDate) {
                // initialize slices
                mValue = new KTSlice(*mVTree, "mValue", mProtCurveName);
                mTmpSlice = new KTSlice(*mVTree, "mTmpSlice", mProtCurveName);
            }
            if (mProtLeg->mNtlStartDts[idx]==currDate /*&& mProtZero[idx]!=NULL*/) {

                /* Get protection value from bank and do the (1-R)*ntl thing */
                mVTree->Get(*mTmpSlice, *mProtZero[idx]);
                if (!mProtLeg->mIsDefRecovery) (*mTmpSlice) *= (1. - mProtLeg->mRecovery);
                (*mTmpSlice) *= mProtLeg->mNtlAmts[idx];

                /* Add to mValue */
                (*mValue) += (*mTmpSlice);
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
// Returns current PV of floating cashflows after tpIdx.
// This routine adjusts the PV depending on the stub convention supplied.

KTSlice&
KVPToolProtLeg::GetValue()
{
static        char        routine[] = "KVPToolProtLeg::GetValue";

        TDate   currDate = mVTree->TPDateCurrent();
        int n = mProtLeg->mNtlEndDts.size();
        int idx;

        // Before the end date,  mValue=NULL,
        // set the default value to 0.
        if (mDefValue == NULL) {
            mDefValue = new KTSlice(*mVTree, "mDefValue", mProtCurveName);
            ASSERT_OR_THROW(mDefValue != NULL);
            *mDefValue = 0e0;
        }

        if (mValue==NULL) {
            mDefValue->SetTpIdx(mVTree->TPIdxCurrent());
            return *mDefValue;
        }

        /* Iterate through different protection periods and add value to mDefValue if current,
           else will have been done already by update */
        for (idx=0; idx<n; idx++) {

            if (mProtLeg->mNtlEndDts[idx]<=currDate) {
                // this is in the past - ignore it
                continue;

            } else if (currDate<=mProtLeg->mNtlStartDts[idx]) {
                // this value has already been included in mValue, so nothing more to do.
                continue;

            } else {
                // this protection section is in mid-calculation. You should extract its
                // current value and return that plus mValue.

                // Protection zero
                KZeroReset protZero = KZeroReset(mProtCurveName,
                                                currDate, 
                                                mProtLeg->mNtlEndDts[idx]);
                mVTree->Get(*mTmpSlice, protZero); 

                // Need to mulpliply by (1-R) when overwrite
                if (!mProtLeg->mIsDefRecovery)
                    *mTmpSlice *= (1. - mProtLeg->mRecovery);
                // Multiply by notional
                *mTmpSlice *= mProtLeg->mNtlAmts[idx];

                // add current value of rest of component and return.
                (*mTmpSlice) += (*mValue);
                return (*mTmpSlice);
            }
        }
        
        // if you get here you are not in the middle of one of the
        // protection periods, so mValue is the value, as no protection
        // is accruing.
        return (*mValue);
}

