/************************************************************************
 * Module:	dvprod
 * File:
 * Function:    vtlindicator.cxx
 * Author:      Changhong He
 ************************************************************************/
#include "vtlindicator.h"
#include "ktypes.h"
#include "kutilios.h"
#include "kvtspar.h"

extern        KMap(String,double)        globVpConstTable;        // Constant tables

//---------------------------------------------------------------

KVPToolIndicator::KVPToolIndicator(SharedPointer<KIndicator> ins, KVTree &vt)
    : KVPToolAtom(vt)
{
    mIndicator = ins;
    mValue     = NULL;
//    mFormula = "STEP(x0,-999.0,999.0,INSIDE_RANGE,SMOOTH_SINGLE))";
    setFormula(ins->getLBarrier(),
               ins->getHBarrier(),
               ins->getIoWindow(),
               ins->getSmooth());
//    cout << "formula " << mFormula << endl;
}

//---------------------------------------------------------------

KVPToolIndicator::~KVPToolIndicator()
{
    delete mValue;

    if(debugLevel)
        dppLog << GetName() << ": deleted." << endl;
}



//---------------------------------------------------------------


inline const String&
KVPToolIndicator::GetCurveName()
{
    if (mCurveName.empty()) 
    {
 	  try {
	    mCurveName = (mIndicator->CplxRate()).CurveName();
	  }
	  catch (KFailure) {
	    throw KFailure("KVPToolIndicator::GetCurveName: "
		"failed on `%s'.\n", mIndicator->GetName());
	  }
    }
    return mCurveName;
}



//---------------------------------------------------------------

void
KVPToolIndicator::Update()
{

    if (NeedsUpdate()) 
    {
        UpdateDone();
    }
}


//---------------------------------------------------------------

KTSlice&
KVPToolIndicator::GetValue()
{
static  char    routine[] = "KVPToolIndicator::GetValue";

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

    KTSlice *rateTS = new KTSlice(*mVTree, "rateTS", curveName);
    KCplxRateReset cplxRateReset(currDate, mIndicator->CplxRate() );
        

    mVTree->Get(*rateTS, cplxRateReset);

    KTSliceParEval(
                   *mVTree,
                    1,
                    rateTS,
                    mFormula.c_str(),
                    globVpConstTable,        // Constant tables
                    mValue);

/*    cout << routine << " formula " << mFormula << endl;
    dppLog << "**********************" << endl;;
    dppLog << "parser indicator formula " << mFormula << endl;
    dppLog << "cplx rate slice" << endl;
    dppLog << *rateTS;
    dppLog << "mvalue" << endl;
    dppLog << *mValue;
    dppLog << "*********************" << endl;*/
        
    return(*mValue);
}

void
KVPToolIndicator::setFormula(double    LBarrier,
                             double    HBarrier,
                             KKnockIO  IoWindow,
                             KSmooth   Smooth)
{
static  char    routine[] = "KVPToolIndicator::SetFormula";

    try
    {
        mFormula = "STEP(X0,";
        mFormula = mFormula + LBarrier;
        mFormula = mFormula + ",";
        mFormula = mFormula + HBarrier;
        mFormula = mFormula + ",";

        if (IoWindow == CRX_KNOCK_IN)
            mFormula = mFormula + INRANGE;
        else if (IoWindow == CRX_KNOCK_OUT)
            mFormula = mFormula + OUTRANGE;
        else
            throw KFailure ("%s failed: Invalid IoWindow.\n", routine);

        mFormula = mFormula + ",";
        
        if (Smooth == NO_SMOOTH)
            mFormula = mFormula + SMOOTHNO;
        else if (Smooth == SINGLE_SMOOTH)
            mFormula = mFormula + SMOOTHS;
        else if (Smooth == DOUBLE_SMOOTH)
            mFormula = mFormula + SMOOTHD;
        else
            throw KFailure ("%s failed: Invalid smooth.\n", routine);

        mFormula = mFormula + ")";
    }
    catch (KFailure)
    {
	    throw KFailure ("%s failed\n", routine);
    }
        
    return;
}
