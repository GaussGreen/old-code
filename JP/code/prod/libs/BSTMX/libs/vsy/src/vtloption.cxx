/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      C. Daher, David Liu
 * Revision:	$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/crxvprod/src/vtloption.cxx,v 1.1.1.1 2005/06/27 19:23:41 dliu Exp $
 ************************************************************************/
#include "vtloption.h"

#include "kutilios.h"



extern "C" {
#include "drltime.h"
};


//===============================================================
//
//===============================================================


//---------------------------------------------------------------

KVPToolOption::KVPToolOption(SharedPointer<KVPOption> ins, KVTree &vt)
	: KVPToolAtom(vt)
{
static	char	routine[] = "KVPToolOption::KVPToolOption";
	int	idx, n;

    try {
	// Check valid
	if (ins->mSettleDates.size() <= 0) {
		throw KFailure("%s: no dates found.\n", routine);
	}


	//
	//

	mOption = ins;

	n = mOption->mSettleDates.size();
	if (mOption->mAmerican) {
		throw KFailure("%s: amer option not yet implemented.\n",
			routine);
	} else {
		for (idx=0; idx<n; idx++) {
			vt.Insert(mOption->mNotifDates[idx]);
			vt.Insert(mOption->mSettleDates[idx]);
		}
	}

	mValue = NULL;
	mDefValue = NULL;
	mXT0 = NULL;
	mXT1 = NULL;
	mXT2 = NULL;

	// Run statistics
	mRunStat = (debugLevel > 0 ? TRUE : FALSE);

    }
    catch (KFailure) {
	throw KFailure("%s: failed on option schedule:\n", routine);
    }
}

//---------------------------------------------------------------

KVPToolOption::~KVPToolOption()
{
	delete mValue;
	delete mDefValue;
	delete mXT0;
	delete mXT1;
	delete mXT2;

	if(debugLevel)
		dppLog << GetName() << ": deleted." << endl;
}


//---------------------------------------------------------------


inline const String&
KVPToolOption::GetCurveName()
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
            discIdx = mVTree->GetCurveIdx(mOption->GetDiscName());
            if (discIdx > cIdxMax)
                cIdxMax = discIdx;

	    mCurveName = mVTree->GetCurveName(cIdxMax);

	}
	  catch (KFailure) {
	    throw KFailure("KVPToolOption::GetCurveName: "
		"failed on option `%s'.\n", mOption->GetName());
	  }
	}
	return mCurveName;
}




//---------------------------------------------------------------


void
KVPToolOption::Update()
{
static	char	routine[] = "KVPToolOption::Update";
	int	idx, n = mOption->mNotifDates.size();
	TDate	curDate = mVTree->TPDateCurrent(),
		settleDate;
	int	idxD;


const	String&	discCurveName = mOption->mDiscZcName;
const	String&	curveName = GetCurveName();


    try {
	if (!NeedsUpdate()) return;

	double	exTime = mVTree->TPTimeCurrent();

	//
	// Perform dev on the value and on the exercise bank.
	//
	if (mValue != NULL) {
		mValue->Dev(discCurveName);
		if (mRunStat) {
			mXT0->Ev();
			mXT1->Ev();
			mXT2->Ev();
		}
	}

	for(KMap(TDate,KTSlice*)::iterator p=mExer.begin();
	     p != mExer.end();
	     ++p) {
		mExer[(*p).first]->Dev(discCurveName);
	}



	//
	// Add settlements
	//
	for (idx=0; idx<n; idx++) {
	    settleDate = mOption->mSettleDates[idx];

	    if (settleDate == curDate) {

		if (debugLevel) 
			dppLog << GetName() << ": Processing setllement "
				<< idx <<endl;

		//
		// Create new settlement value slice
		//
		mExer[settleDate] = new KTSlice(*mVTree, "Exercise", curveName);
		*(mExer[settleDate]) = (-mOption->mStrikes[idx]);

		//
		// Calculate value of the underlying
		// by summimg all dependencies.
		//
		for (idxD=0; idxD<this->NumDep(); idxD++) {
			(*mExer[settleDate]) += this->Dep(idxD)->GetValue(); 
		}

	    }
	}


	//
	// Perform exercises
	//
	for (idx=0; idx<n; idx++) {
	    //
	    // If is notifiable
	    //
	    if (mOption->mNotifDates[idx] == curDate) {

		if (debugLevel) 
		    dppLog << GetName() << format(
			": Processing notification %d (settlement %s %lf)",
			idx, DrlTDatePrint(NULL, mOption->mSettleDates[idx]),
			exTime) << endl;

		//
		// Ensure value slice exists
		//
		if (mValue == NULL) {
			mValue = new KTSlice(*mVTree, "mValue", curveName);
			*mValue = 0e0;
			if (mRunStat) {
				mXT0 = new KTSlice(*mVTree, "mXT0", curveName);
				*mXT0 = 0e0;
				mXT1 = new KTSlice(*mVTree, "mXT0", curveName);
				*mXT1 =	 0e0;
				mXT2 = new KTSlice(*mVTree, "mXT0", curveName);
				*mXT2 = 0e0;
			}

			/**mXT0 = 1e0;
			*mXT1 = exTime;
			*mXT2 = exTime*exTime;*/
		}


		//
		// Retrieve settlement value time slice
		//
		settleDate = mOption->mSettleDates[idx];

		KMap(TDate,KTSlice*)::iterator p = mExer.find(settleDate);

		if (p == mExer.end()) {
			throw KFailure("%s: can't find settle date %s.\n",
				routine, DrlTDatePrint(NULL, settleDate));
		}

		KTSlice& exerTs = *(mExer[(*p).first]);


		//
		// Perform exersize
		//

		switch (mOption->mType) {
		case KVPOption::CALL:
	        	if (mRunStat) {
				mVTree->TSliceSpecialOper(exerTs,
				"EXTIMESET",
				(KTSlice*) &exerTs,
				(KTSlice*) mValue,
				(KTSlice*) mXT0,
				(KTSlice*) mXT1,
				(KTSlice*) mXT2,
				exTime);
			}

			mValue->max(exerTs);
			break;
		case KVPOption::PUT:
			exerTs *= (-1e0);
	        	if (mRunStat) {
				mVTree->TSliceSpecialOper(exerTs,
				"EXTIMESET",
				(KTSlice*) &exerTs,
				(KTSlice*) mValue,
				(KTSlice*) mXT0,
				(KTSlice*) mXT1,
				(KTSlice*) mXT2,
				exTime);
			}

			mValue->max(exerTs);
			break;
		case KVPOption::TIMING_CALL:
		case KVPOption::TIMING_PUT:
			// For forward, if it is the last date
			// we force exercise.
			THROW_NA;
			break;
		}

		//
		// Erase exerTs
		//
		delete mExer[(*p).first];
		mExer.erase(p);

	    }
	}



	//
	// Last TP: store the result.
	//
	if (mVTree->TPIdxCurrent() == 0) {
		if (mValue == NULL) {
			throw KFailure("%s: no exercise event "
				"after today date.\n", routine);
		}
		// Store PV
		mResults.insert(
			String_Double("PV", mValue->GetCenter())
		);

#ifdef	__SKIP__	//$$$ NEED TO BE FIXED (CD)

		// Calculate PV value of the underlying
		// by summimg all dependencies.
		//
		double	underlyingValue = 0e0;
		for (idxD=0; idxD<this->NumDep(); idxD++) {
			underlyingValue += Dep(idxD)->Results()["PV"];
		}
		mResults.insert(
			String_Double("UNDER", underlyingValue)
		);
#endif
		mResults.insert(
			String_Double("UNDER", 0e0)
		);

		// Expected time to exercise
		//
		if (mRunStat) {
			double	xt0 = mXT0->GetCenter(),
				xt1 = mXT1->GetCenter(),
				xt2 = mXT2->GetCenter(),
				pe,			// prob exercise
				te,			// expected exer time
				se,			// var exer time
				tZ;			// last exer time
			TDate	lastDate =
			    mOption->mNotifDates[mOption->mNotifDates.size()-1];

			tZ = (lastDate - mVTree->TPToday()) / 365e0;
			pe = xt0;
			te = xt1 + (1-pe)*tZ;

			mResults["XTPROB"] = pe;
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
KVPToolOption::GetValue()
{
	// Before any exercise event occurs, mValue=NULL,
	// set the option to default value 0.
	if (mDefValue == NULL) {
		mDefValue = new KTSlice(*mVTree, "mDefValue", GetCurveName());
		ASSERT_OR_THROW(mDefValue != NULL);
		*mDefValue = 0e0;
	}

	if (mValue == NULL) {
		mDefValue->SetTpIdx(mVTree->TPIdxCurrent());	
		return(*mDefValue);
	}
	else	
		return(*mValue);
}



