/************************************************************************
 * Module:	dvprod
 * File:
 * Function:    vtlcplxrate.cxx
 * Author:      Changhong He
 ************************************************************************/
#include "vtlcplxrate.h"

#include "kutilios.h"


//---------------------------------------------------------------

KVPToolCplxRate::KVPToolCplxRate(SharedPointer<KCplxRate> ins, KVTree &vt)
    : KVPToolAtom(vt)
{
    mCplxRate = ins;
    mValue    = NULL;
}

//---------------------------------------------------------------

KVPToolCplxRate::~KVPToolCplxRate()
{
    delete mValue;

    if(debugLevel)
        dppLog << GetName() << ": deleted." << endl;
}



//---------------------------------------------------------------


inline const String&
KVPToolCplxRate::GetCurveName()
{
    if (mCurveName.empty()) 
    {
 	  try {
	    mCurveName = mCplxRate->CurveName();
        cout << "KVPToolCplxRate " << mCplxRate->GetName() << '\t';
        cout << "curve name " << mCplxRate->CurveName() << endl;
	  }
	  catch (KFailure) {
	    throw KFailure("KVPToolRate::GetCurveName: "
		"failed on `%s'.\n", mCplxRate->GetName());
	  }
//      mCurveName = K_DEFAULT_NAME;
    }
    return mCurveName;
}



//---------------------------------------------------------------

void
KVPToolCplxRate::Update()
{

    if (NeedsUpdate()) 
    {
        UpdateDone();
    }
}


//---------------------------------------------------------------

KTSlice&
KVPToolCplxRate::GetValue()
{
static  char    routine[] = "KVPToolCplxRate::GetValue";

    TDate   currDate = mVTree->TPDateCurrent();

const   String& curveName = GetCurveName();

    // If the storage slice has not been allocated,
    // allocate it now.
    if (mValue == NULL) {
        mValue = new KTSlice(*mVTree,
            GetName() + String("::mValue"),
            curveName);
    }

    *mValue = 0e0;


    KCplxRateReset cplxRateReset(currDate, *mCplxRate);

    mVTree->Get(*mValue, cplxRateReset);

    return(*mValue);
}

