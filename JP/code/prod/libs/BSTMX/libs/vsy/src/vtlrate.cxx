/************************************************************************
 * Module:	dvprod
 * File:
 * Function:    vtlrate.cxx
 * Author:      David Liu
 ************************************************************************/
#include "vtlrate.h"

#include "kutilios.h"


//---------------------------------------------------------------

KVPToolRate::KVPToolRate(SharedPointer<KRate> ins, KVTree &vt)
	: KVPToolAtom(vt)
{
	mRate = ins;
	mValue = NULL;
}

//---------------------------------------------------------------

KVPToolRate::~KVPToolRate()
{
	delete mValue;

	if(debugLevel)
		dppLog << GetName() << ": deleted." << endl;
}



//---------------------------------------------------------------


inline const String&
KVPToolRate::GetCurveName()
{
	if (mCurveName.empty()) {
	  try {
	    mCurveName = mRate->CurveName();
	  }
	  catch (KFailure) {
	    throw KFailure("KVPToolRate::GetCurveName: "
		"failed on `%s'.\n", mRate->GetName());
	  }
	}
	return mCurveName;
}



//---------------------------------------------------------------

void
KVPToolRate::Update()
{

	if (NeedsUpdate()) {

		UpdateDone();
	}
}


//---------------------------------------------------------------

KTSlice&
KVPToolRate::GetValue()
{
static	char	routine[] = "KVPToolRate::GetValue";

	TDate	currDate = mVTree->TPDateCurrent();

const	String& curveName = GetCurveName();

	// If the storage slice has not been allocated,
	// allocate it now.
	if (mValue == NULL) {
	    mValue = new KTSlice(*mVTree,
			GetName() + String("::mValue"),
			curveName);
	}

	*mValue = 0e0;


	KRateReset rateReset(currDate, *mRate);

	mVTree->Get(*mValue, rateReset);

	return(*mValue);
}


