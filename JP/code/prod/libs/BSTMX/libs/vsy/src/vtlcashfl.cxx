/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      C. Daher, David Liu
 * Revision:	$Header$
 ************************************************************************/
#include "vtlcashfl.h"

#include "kutilios.h"


extern "C" {
#include "drltime.h"
};


//===============================================================
//
//===============================================================


//---------------------------------------------------------------

KVPToolCashFlows::KVPToolCashFlows(SharedPointer<KVPCashFlows> ins, KVTree &vt)
	: KVPToolAtom(vt)
{
	int	idx, numCf = 0;

	mCashFlows = ins;

	ASSERT_OR_THROW(mCashFlows->mPayDates.size()
				== mCashFlows->mPayAmounts.size());
	for (idx=0; idx<mCashFlows->mPayDates.size(); idx++) {
		if (mCashFlows->mPayDates[idx] <= vt.TPToday()) {
			continue;
		}
		vt.Insert(mCashFlows->mPayDates[idx]);
		numCf++;
	}
	if (numCf == 0) {
		throw KFailure("KVPToolCashFlows: no live cash flows found "
				"for %s.\n",
				ins->GetName());
	}

	mValue = NULL;
	mDefValue = NULL;
	
}

//---------------------------------------------------------------

KVPToolCashFlows::~KVPToolCashFlows()
{
	delete mValue;
	delete mDefValue;
	if (debugLevel) 
		dppLog << GetName() << ": deleted." << endl;
}


//---------------------------------------------------------------


void
KVPToolCashFlows::Update()
{
	int	idx;
	TDate	curDate = mVTree->TPDateCurrent();

	if (!NeedsUpdate()) return;

	//
	// We discount using the curve name of the product.
	//
	if (mValue != NULL)
		mValue->Dev(mCashFlows->GetDiscName());



	// Add cash flows.
	for (idx=0; idx<mCashFlows->mPayDates.size(); idx++) {
	    if (mCashFlows->mPayDates[idx] == curDate) {

		if (debugLevel) 
		    dppLog << GetName() << format(
			": Adding cash flow %d (%16.8f)",
			idx, mCashFlows->mPayAmounts[idx]) << endl;

		if (mValue == NULL) {
			mValue = new KTSlice(*mVTree, "Cfl dev", 
					     mCashFlows->GetDiscName());

			*mValue = 0e0;
		}

		*mValue += mCashFlows->mPayAmounts[idx];
	    }
	}

	UpdateDone();
}


//---------------------------------------------------------------

KTSlice&
KVPToolCashFlows::GetValue()
{
static	char	routine[] = "KVPToolCashFlows::GetValue";

	// Before any cash pay dates, mValue=NULL,
	// set the default value to 0.
        if (mDefValue == NULL) {
                mDefValue = new KTSlice(*mVTree, 
					"mDefValue", 
					mCashFlows->GetDiscName());
                ASSERT_OR_THROW(mDefValue != NULL);
                *mDefValue = 0e0;
        }

	if (mValue == NULL)
	{
		mDefValue->SetTpIdx(mVTree->TPIdxCurrent());
		return (*mDefValue);
	}
	else
		return(*mValue);
}

