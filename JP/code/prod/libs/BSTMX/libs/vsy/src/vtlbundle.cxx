/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      C. Daher, David Liu
 * Revision:	$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/crxvprod/src/vtlbundle.cxx,v 1.1.1.1 2005/06/27 19:23:40 dliu Exp $
 ************************************************************************/
#include "vtlbundle.h"

#include "kutilios.h"


//---------------------------------------------------------------

KVPToolWBundle::KVPToolWBundle(SharedPointer<KVPWBundle> ins, KVTree &vt)
	: KVPToolAtom(vt)
{
	mWBundle = ins;
	mValue = NULL;
}

//---------------------------------------------------------------

KVPToolWBundle::~KVPToolWBundle()
{
	delete mValue;

	if(debugLevel)
		dppLog << GetName() << ": deleted." << endl;
}



//---------------------------------------------------------------


inline const String&
KVPToolWBundle::GetCurveName()
{
	if (mCurveName.empty()) {
	  try {
	    int	idx, idxMax, cIdxMax, cIdx;
            int discIdx;

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
		if (cIdx > cIdxMax) {
			idxMax = idx;
			cIdxMax = cIdx;
		}
	    }

            // include discount curve
            discIdx = mVTree->GetCurveIdx(mWBundle->GetDiscName());
            if (discIdx > cIdxMax)
                cIdxMax = discIdx;

	    mCurveName = mVTree->GetCurveName(cIdxMax);
	  }
	  catch (KFailure) {
	    throw KFailure("KVPToolWBundle::GetCurveName: "
		"failed on `%s'.\n", mWBundle->GetName());
	  }
	}
	return mCurveName;
}



//---------------------------------------------------------------

void
KVPToolWBundle::Update()
{

	if (NeedsUpdate()) {

		UpdateDone();
	}
}


//---------------------------------------------------------------

KTSlice&
KVPToolWBundle::GetValue()
{
static	char	routine[] = "KVPToolWBundle::GetValue";
	int	idx;
	double	weight;
	KTSlice *tmpSlice = NULL;

const	String& curveName = GetCurveName();

	// If the storage slice has not been allocated,
	// allocate it now.
	if (mValue == NULL) {
	    mValue = new KTSlice(*mVTree,
			GetName() + String("::mValue"),
			curveName);
	}

	*mValue = 0e0;

	//
	// Temprary slice to store the underlying value, so that the value of underlying
	// slice is NOT changed during the operation.
	// 
	tmpSlice = new KTSlice(*mVTree,
				"Temp slice to store the underlying value",
				curveName);

        TDate currDate = mVTree->TPDateCurrent();

	for(idx=0; idx<NumDep(); idx++) {
		KTSlice	&ts = Dep(idx)->GetValue();
		weight = mWBundle->mWeights[idx];

		// Copy to the tmp slice and apply the weighting to tmp slice,
		// NOT to the original underlying slice.
		//
		*tmpSlice = ts;
		*tmpSlice *= weight;

		*mValue += *tmpSlice;

	}

	delete tmpSlice;

	return(*mValue);
}



