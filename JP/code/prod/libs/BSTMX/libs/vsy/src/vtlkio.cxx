/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      D. Liu
 ************************************************************************/
#include "vtlkio.h"

#include "kutilios.h"
#include "kmrntree.h"



extern "C" {
#include "drltime.h"
#include "drlmath.h"
};


//===============================================================
//
//===============================================================


//---------------------------------------------------------------

KVPToolKnockIO::KVPToolKnockIO(SharedPointer<KVPKnockIO> ins, KVTree &vt)
	: KVPToolAtom(vt)
{
static	char	routine[] = "KVPToolKnockIO::KVPToolKnockIO";
	int	idx, n;

	TDate	today = vt.TPToday();

    try {
	// Check valid
	if (ins->mSettleDates.size() <= 0) {
		throw KFailure("%s: no dates found.\n", routine);
	}

	/* Trigger rate index curve name */
	mRateIdxName = ins->mRateIndex[0]->Rate().CurveName(); 

	// Only include events with observation date >= today 
	//
	ins->ValidEvents(today);

	mKnockIO = ins;

	// Add events in the tree.
	//
	n = mKnockIO->mObsDates.size();

	mModelObsDates.resize(n);

	for (idx=0; idx<n; idx++) {
		vt.Insert(mKnockIO->mObsDates[idx]);
		vt.Insert(mKnockIO->mSettleDates[idx]);
		
		mModelObsDates[idx] = vt.Insert(*(mKnockIO->mRateIndex[idx]));

	}

	mValue = NULL;
	mValueUS  = NULL;
	mGetValue = NULL;
        mDefValue = NULL;

	mRebateValue   = NULL;
	mRebateValueUS = NULL;

	mXT0 = NULL;
	mXT1 = NULL;
	mXT2 = NULL;

	// Run statistics
	mRunStat = (debugLevel > 0 ? TRUE : FALSE);


    }
    catch (KFailure) {

	throw KFailure("%s: failed on knock IO schedule:\n", routine);
    }
}

//---------------------------------------------------------------

KVPToolKnockIO::~KVPToolKnockIO()
{

	delete mValue;
	delete mValueUS;
	delete mGetValue;

	delete mRebateValue;
	delete mRebateValueUS;
    delete mDefValue;

	delete mXT0;
	delete mXT1;
	delete mXT2;

	for (KMap(TDate,KTSlice*)::iterator itEx = mExer.begin();
			itEx != mExer.end(); ++itEx)
	{
		delete (*itEx).second;
	}
	mExer.clear();

	for (KMap(TDate,KTSlice*)::iterator itRebate = mRebate.begin();
			itRebate != mRebate.end(); ++itRebate)
	{
		delete (*itRebate).second;
	}
	mRebate.clear();

	if(debugLevel)
		dppLog << GetName() << ": deleted." << endl;
}


//---------------------------------------------------------------


const String&
KVPToolKnockIO::GetCurveName()
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
            discIdx = mVTree->GetCurveIdx(mKnockIO->GetDiscName());
            if (discIdx > cIdxMax)
                cIdxMax = discIdx;

	    mCurveName = mVTree->GetCurveName(cIdxMax);
	  }
	  catch (KFailure) {
	    throw KFailure("KVPToolKnockIO::GetCurveName: "
		"failed on knock IO `%s'.\n", mKnockIO->GetName());
	  }
	}
	return mCurveName;
}




//---------------------------------------------------------------


void
KVPToolKnockIO::Update()
{
static	char	routine[] = "KVPToolKnockIO::Update";
	int	idx, n = mKnockIO->mObsDates.size();
	TDate	curDate = mVTree->TPDateCurrent(),
		settleDate, obsDate, modelObsDate;
	int	idxD;

	int	depCurveIdx, rateCurveIdx;

	bool	isKnockOut = false;

	double	idxRate, idxStep;
	double	value;				// Smoothed node value.
	double	valueUS;			// Unsmoothed payoff value
	double	rebateUS;			// Unsmoothed rebate value
	double	exerNodeValue, rebate;
	double	lowerStep, upperStep, rangeInStep; // Smoothing steps

	KRateReset      *indexRate = NULL;

	KTSlice		*rateTS = NULL;		// Knock IO rate index slice
	KTSliceNode	tsNode;			// Slice node

	double 		underlyingValue = 0e0;

	double	xt0,		// Probability
		xt1,		// expected ko time
		xt2;		// expected ko time^2 

const	String&	discCurveName = mKnockIO->mDiscZcName;

	String	curveName = GetCurveName();	// max dim of underlying slice
	String	rateIndexCV = mRateIdxName; 


 try {
	if (!NeedsUpdate()) return;

	double	exTime = mVTree->TPTimeCurrent();


	// 
	// The slice dimensions of mValue and mExer should be
	// the maximum of rate slice and underlying slice
	//
	depCurveIdx  = mVTree->GetCurveIdx(curveName);
	rateCurveIdx = mVTree->GetCurveIdx(mRateIdxName);

	if (rateCurveIdx > depCurveIdx)
		curveName = mRateIdxName;


	// Set Knock out flag
	//
	if (mKnockIO->mIOType == CRX_KNOCK_OUT)
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
	    settleDate = mKnockIO->mSettleDates[idx];

	    if (settleDate == curDate) {

		if(debugLevel)
		    dppLog << GetName() << format(": Processing settlement %d "
			      " (settlement %s)",
			      idx,
			      DrlTDatePrint(NULL, mKnockIO->mSettleDates[idx]))
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
			(*mRebate[settleDate]) = mKnockIO->mRebates[idx];
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

	    } 	// if settleDate
	}	// for idx


	//------------------------------------------------------
	//
	// (3) Perform knock-in exercises.
	//
	//------------------------------------------------------

	for (idx=0; idx<n; idx++) {
	    obsDate      = mKnockIO->mObsDates[idx];
	    modelObsDate = mModelObsDates[idx];

	    if (modelObsDate == curDate) {

		if(debugLevel)
		    dppLog << GetName() << format(": Processing observation "
			  "%d on date %s (settlement %s)", 
			  idx,
			  DrlTDatePrint(NULL, curDate), 
			  DrlTDatePrint(NULL, mKnockIO->mSettleDates[idx]))
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

		if (mKnockIO->mSmooth == DOUBLE_SMOOTH && mValueUS == NULL) {
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

		    if (mRebateValueUS == NULL && mKnockIO->mSmooth == DOUBLE_SMOOTH)
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

			if (mKnockIO->mSmooth == DOUBLE_SMOOTH) {
				dppLog << *mValueUS << endl;
			}

			if (isKnockOut)
				dppLog << *mRebateValue << endl;
		}


		//
		// Retrieve settlement value time slice
		//
		settleDate = mKnockIO->mSettleDates[idx];

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
		// Perform exersize
		//

		//
		// Create new observation rate slice
		//
		rateIndexCV = (mKnockIO->mRateIndex[idx])->Rate().CurveName(); 

		rateTS = new KTSlice(*mVTree, "Observ", rateIndexCV);

		mVTree->Get(*rateTS, *(mKnockIO->mRateIndex[idx]));
		
			
		
		if (debugLevel == DEBUG_LEVEL_KIO)
		{
			dppLog << "Knock IO rate index slice" << endl;
			dppLog << *rateTS << endl;

			dppLog << "Underlying slice" << endl;
			dppLog << exerTS << endl;
		}




		// Check each node for knock IO
		//
		if (mKnockIO->mSmooth == NO_SMOOTH)
		{
		    if(mKnockIO->mIOWindow == CRX_KNOCK_IN)
		    {
		    	for(tsNode.begin(*mValue); !tsNode.end(); ++tsNode)
		    	{
		    	   idxRate = (*rateTS)[tsNode];
			
/*
			   if (debugLevel == DEBUG_LEVEL_KIO)
			   {
  			   	dppLog << tsNode << endl;

			   	dppLog << "KI: Index rate: " << idxRate << endl;
			   	dppLog << "Low Barrier: " 
				       << mKnockIO->mBarrierLo[idx] << endl;
			   	dppLog << "High Barrier: " 
				       << mKnockIO->mBarrierHi[idx] << endl;
			   }
*/

			   // If within the range, = current underlying value,
			   // else, unchanged.
			   //
		    	   if(idxRate > mKnockIO->mBarrierLo[idx]*
					(1. - BARRIER_TOL)         &&
		       	      idxRate < mKnockIO->mBarrierHi[idx]*
					(1. + BARRIER_TOL))
		    	   {
			    	exerNodeValue = exerTS[tsNode];

			    	(*mValue)[tsNode] = exerNodeValue;
				if (mRunStat) {
					(*mXT0)[tsNode] = 1e0;
					(*mXT1)[tsNode] = exTime;
					(*mXT2)[tsNode] = exTime*exTime;
				}

				if (isKnockOut)
				    (*mRebateValue)[tsNode] 
						= (*rebateTS)[tsNode];

			   	if (debugLevel == DEBUG_LEVEL_KIO)
			   	{
					dppLog << "KIO: underlying value: " 
					       << exerNodeValue << endl;
					dppLog << *mValue << endl;
				}

		    	   }	
			}
		    }
		    else if(mKnockIO->mIOWindow == CRX_KNOCK_OUT)
		    {
		    	for(tsNode.begin(*mValue); !tsNode.end(); ++tsNode)
		    	{
		    	   idxRate = (*rateTS)[tsNode];

/*
			   if (debugLevel == DEBUG_LEVEL_KIO)
			   {
  			   	dppLog << tsNode << endl;

			   	dppLog << "KO: Index rate: " << idxRate << endl;
			   	dppLog << "Low Barrier: " 
				       << mKnockIO->mBarrierLo[idx] << endl;
			   	dppLog << "High Barrier: " 
				       << mKnockIO->mBarrierHi[idx] << endl;
			   }
*/


			   // If outside the range, = current underlying value,
			   // else, unchanged.
			   //
		           if(idxRate < mKnockIO->mBarrierLo[idx]*
					(1. + BARRIER_TOL)         ||
		              idxRate > mKnockIO->mBarrierHi[idx]*
					(1. - BARRIER_TOL))
		    	   {
			    	exerNodeValue = exerTS[tsNode];

				(*mValue)[tsNode] = exerNodeValue;
				if (mRunStat) {
					(*mXT0)[tsNode] = 1e0;
					(*mXT1)[tsNode] = exTime;
					(*mXT2)[tsNode] = exTime*exTime;
				}

				if (isKnockOut)
				    (*mRebateValue)[tsNode] 
						= (*rebateTS)[tsNode];

/*
			   	if (debugLevel == DEBUG_LEVEL_KIO)
			   	{
					dppLog << "KIO: underlying value: " 
					       << exerNodeValue << endl;
					dppLog << *mValue << endl;
				}
*/

		    	   }	
			}
		    }
		    else
		    	throw KFailure("%s: invalid knock window type (%d).\n",
				       routine,
				       mKnockIO->mIOWindow);	



		}
		else if (mKnockIO->mSmooth == SINGLE_SMOOTH)	// single smoothing
		{

		    if(mKnockIO->mIOWindow == CRX_KNOCK_IN)
		    {
		    	for(tsNode.begin(*mValue); !tsNode.end(); ++tsNode)
		    	{
			   idxStep = rateTS->GetNodeStepMax(tsNode);
		    	   idxRate = (*rateTS)[tsNode];

			   exerNodeValue = exerTS[tsNode];
			   value         = (*mValue)[tsNode];
			   if (mRunStat) {
			   	xt0           = (*mXT0)[tsNode];
			   	xt1		 = (*mXT1)[tsNode];
			   	xt2 		 = (*mXT2)[tsNode];
			   }

			   ////////////////////////////////////////////
			   // Given two step functions:
			   // 
			   //               f1
			   // 
			   //            ------------------------
			   //            |         
			   //            |         
			   //    --------         
			   //
			   //
			   //               f2
			   // 
			   //                      --------------
			   //                      |         
			   //                      |         
			   //    ------------------         
			   //
			   //
			   // Knock-in range payoff function
			   //
			   //             (1-f2)*f1
			   //
			   //            -----------
			   //            |         |
			   //            |         |
			   //    --------          -------------- 
			   //
			   //////////////////////////////////////////////

			   lowerStep = DrlSmoothStepFcn(
					idxRate - mKnockIO->mBarrierLo[idx],
					idxStep);
			   upperStep = DrlSmoothStepFcn(
					idxRate - mKnockIO->mBarrierHi[idx],
					idxStep);

			   rangeInStep = (1. - upperStep) * lowerStep;

			   // 
			   // Knock-in window:
			   // If < lower barrier, unchanged,
			   // if > higher barrier, unchanged.
			   // if lo barrier < && < hi barrier, underlying.
			   //

			   value = (exerNodeValue - value) * rangeInStep
				   + value;

			   (*mValue)[tsNode] = value;
			   if (mRunStat) {
			   (*mXT0)[tsNode] = (1e0 - xt0) * rangeInStep + xt0;
			   (*mXT1)[tsNode] = (exTime - xt1) * rangeInStep + xt1;
			   (*mXT2)[tsNode] = (exTime*exTime - xt2)*rangeInStep
					     + xt2;
			   }


			   if (isKnockOut)
			   	(*mRebateValue)[tsNode] 
				    = ( (*rebateTS)[tsNode] 
				      - (*mRebateValue)[tsNode]
				      ) * rangeInStep
				    + (*mRebateValue)[tsNode];

		    	}	
		    }
		    else if(mKnockIO->mIOWindow == CRX_KNOCK_OUT)
		    {
		    	for(tsNode.begin(*mValue); !tsNode.end(); ++tsNode)
		    	{
			   idxStep = rateTS->GetNodeStepMax(tsNode);
		    	   idxRate = (*rateTS)[tsNode];

/*
			   if (debugLevel == DEBUG_LEVEL_KIO)
			   {
  			   	dppLog << tsNode << endl;

			   	dppLog << "KO: Index rate: " << idxRate << endl;
			   	dppLog << "KO: Index MaxDiff: " << idxStep 
				       << endl;
			   }
*/

			   exerNodeValue = exerTS[tsNode];
			   value         = (*mValue)[tsNode];
			   if (mRunStat) {
			   xt0           = (*mXT0)[tsNode];
			   xt1		 = (*mXT1)[tsNode];
			   xt2 		 = (*mXT2)[tsNode];
			   }


			   ////////////////////////////////////////////
			   // Given two step functions:
			   // 
			   //               f1
			   // 
			   //            ------------------------
			   //            |         
			   //            |         
			   //    --------         
			   //
			   //
			   //               f2
			   // 
			   //                      --------------
			   //                      |         
			   //                      |         
			   //    ------------------         
			   //
			   //
			   // Knock-in range payoff function
			   //
			   //             (1-f2)*f1
			   //
			   //            -----------
			   //            |         |
			   //            |         |
			   //    --------          -------------- 
			   //
			   //
			   // Knock-out range payoff function
			   //
			   //            1 - (1-f2)*f1
			   //
			   //    --------          -------------- 
			   //            |         |
			   //            |         |
			   //            -----------
			   //
			   //////////////////////////////////////////////

			   lowerStep = DrlSmoothStepFcn(
					idxRate - mKnockIO->mBarrierLo[idx],
					idxStep);
			   upperStep = DrlSmoothStepFcn(
					idxRate - mKnockIO->mBarrierHi[idx],
					idxStep);
			
			   rangeInStep = (1. - upperStep) * lowerStep; 

			   // 
			   // Knock-out window:
			   // If < lower barrier, underlying,
			   // if > higher barrier, underlying.
			   // if lower barrier < && < higher barrier, unchanged.
			   //
			   value = (value - exerNodeValue) * rangeInStep
				   + exerNodeValue;

			   (*mValue)[tsNode] = value;
			   if (mRunStat) {
			   (*mXT0)[tsNode] = (xt0 - 1e0) * rangeInStep + 1e0;
			   (*mXT1)[tsNode] = (xt1 - exTime) * rangeInStep 
					     + exTime;
			   (*mXT2)[tsNode] = (xt2 - exTime*exTime)*rangeInStep
					     + exTime*exTime;
			   }

			   if (isKnockOut)
			   	(*mRebateValue)[tsNode] 
				    = ( (*mRebateValue)[tsNode]
				       -(*rebateTS)[tsNode] 
				      ) * rangeInStep
				    + (*rebateTS)[tsNode];

/*
			   if (debugLevel == DEBUG_LEVEL_KIO)
			   {
			   	dppLog << "Range Out Step: " << rangeInStep 
				       << endl;

			   	dppLog << "KO: underlying value: " 
				       << exerNodeValue << endl;

			   	dppLog << *mValue << endl;
			   }
*/

		    	}	
		    }
		    else
		    	throw KFailure("%s: invalid knocking type (%d).\n",
				       routine,
				       mKnockIO->mIOWindow);	



		}    // if single smoothing
		else if (mKnockIO->mSmooth == DOUBLE_SMOOTH)	// double smoothing
		{

		    if(mKnockIO->mIOWindow == CRX_KNOCK_IN)
		    {
		    	for(tsNode.begin(*mValue); !tsNode.end(); ++tsNode)
		    	{
			   idxStep = rateTS->GetNodeStepMax(tsNode);
		    	   idxRate = (*rateTS)[tsNode];

			   exerNodeValue = exerTS[tsNode];
			   value         = (*mValue)[tsNode];
			   valueUS       = (*mValueUS)[tsNode];

			   if (mRunStat) {
			   	xt0           = (*mXT0)[tsNode];
			   	xt1		 = (*mXT1)[tsNode];
			   	xt2 		 = (*mXT2)[tsNode];
			   }


			   if (isKnockOut)
			   {
			   	rebate   = (*mRebateValue)[tsNode];
			   	rebateUS = (*mRebateValueUS)[tsNode];
			   }

			   lowerStep = DrlSmoothStepFcn(
					idxRate - mKnockIO->mBarrierLo[idx],
					idxStep);
			   upperStep = DrlSmoothStepFcn(
					idxRate - mKnockIO->mBarrierHi[idx],
					idxStep);

			   rangeInStep = (1. - upperStep) * lowerStep;

			   // 
			   // Knock-in window:
			   // If < lower barrier, unchanged,
			   // if > higher barrier, unchanged.
			   // if lo barrier < && < hi barrier, underlying.
			   //

			   // Use unsmoothed value within smoothed region
			   // and smoothed value outside smoothed region
			   //
			   if (rangeInStep > 0.0 && rangeInStep < 1.0) 
			   {
			   	value = (exerNodeValue - valueUS) * rangeInStep
				   	+ valueUS;

			   	if (isKnockOut)
					rebate = ((*rebateTS)[tsNode]-rebateUS)
				       		* rangeInStep + rebateUS;
			   }
			   else
			   {
			   	value = (exerNodeValue - value) * rangeInStep
				   	+ value;

			   	if (isKnockOut)
					rebate = ((*rebateTS)[tsNode] - rebate)
				       		* rangeInStep + rebate;
			   }

			   (*mValue)[tsNode]  = value;
			   if (mRunStat) {
			   (*mXT0)[tsNode] = (1e0 - xt0) * rangeInStep + xt0;
			   (*mXT1)[tsNode] = (exTime - xt1) * rangeInStep + xt1;
			   (*mXT2)[tsNode] = (exTime*exTime - xt2)*rangeInStep
					     + xt2;
			   }

			   if (isKnockOut)
			   	(*mRebateValue)[tsNode] = rebate;


			   // Update unsmoothed value slice
			   //
		    	   if(idxRate > mKnockIO->mBarrierLo[idx]*
					(1. - BARRIER_TOL)         &&
		       	      idxRate < mKnockIO->mBarrierHi[idx]*
					(1. + BARRIER_TOL))
		    	   {
			       	(*mValueUS)[tsNode]  = exerNodeValue;	
				if (mRunStat) {
					(*mXT0)[tsNode] = 1e0;
					(*mXT1)[tsNode] = exTime;
					(*mXT2)[tsNode] = exTime*exTime;
				}

			   	if (isKnockOut)
			       	    	(*mRebateValueUS)[tsNode] 
						= (*rebateTS)[tsNode];
			   }


		    	}	
		    }
		    else if(mKnockIO->mIOWindow == CRX_KNOCK_OUT)
		    {
		    	for(tsNode.begin(*mValue); !tsNode.end(); ++tsNode)
		    	{
			   idxStep = rateTS->GetNodeStepMax(tsNode);
		    	   idxRate = (*rateTS)[tsNode];

/*
			   if (debugLevel == DEBUG_LEVEL_KIO)
			   {
  			   	dppLog << tsNode << endl;

			   	dppLog << "KO: Index rate: " << idxRate << endl;
			   	dppLog << "KO: Index MaxDiff: " << idxStep 
				       << endl;
			   }
*/

			   exerNodeValue = exerTS[tsNode];
			   value         = (*mValue)[tsNode];
			   valueUS       = (*mValueUS)[tsNode];

			   if (mRunStat) {
			  	xt0           = (*mXT0)[tsNode];
			   	xt1		 = (*mXT1)[tsNode];
			   	xt2 		 = (*mXT2)[tsNode];
			   }


			   if (isKnockOut)
			   {
			   	rebate   = (*mRebateValue)[tsNode];
			   	rebateUS = (*mRebateValueUS)[tsNode];
			   }


			   lowerStep = DrlSmoothStepFcn(
					idxRate - mKnockIO->mBarrierLo[idx],
					idxStep);
			   upperStep = DrlSmoothStepFcn(
					idxRate - mKnockIO->mBarrierHi[idx],
					idxStep);
			
			   rangeInStep = (1. - upperStep) * lowerStep; 

			   // 
			   // Knock-out window:
			   // If < lower barrier, underlying,
			   // if > higher barrier, underlying.
			   // if lower barrier < && < higher barrier, unchanged.
			   //

			   // Use unsmoothed value within smoothed region
			   // and smoothed value outside smoothed region
			   //
			   if (rangeInStep > 0.0 && rangeInStep < 1.0) 
			   {
			   	value = (valueUS - exerNodeValue) * rangeInStep
				       + exerNodeValue;

			   	if (isKnockOut)
					rebate = (rebateUS-(*rebateTS)[tsNode])
				       	    * rangeInStep + (*rebateTS)[tsNode];
			   }
			   else
			   {
			   	value = (value - exerNodeValue) * rangeInStep
				   	+ exerNodeValue;

			   	if (isKnockOut)
					rebate = (rebate - (*rebateTS)[tsNode])
				       	    * rangeInStep + (*rebateTS)[tsNode];
			   }

			   (*mValue)[tsNode]  = value;
			   if (mRunStat) {
			   (*mXT0)[tsNode] = (xt0 - 1e0) * rangeInStep + 1e0;
			   (*mXT1)[tsNode] = (xt1 - exTime) * rangeInStep 
					     + exTime;
			   (*mXT2)[tsNode] = (xt2 - exTime*exTime)*rangeInStep
					     + exTime*exTime;
			   }

			   if (isKnockOut)
			   	(*mRebateValue)[tsNode] = rebate;

/*
			   if (debugLevel == DEBUG_LEVEL_KIO)
			   {
			   	dppLog << "Range Out Step: " << rangeInStep 
				       << endl;
			   	dppLog << "KO: underlying value: " 
				       << exerNodeValue << endl;

			   	dppLog << *mValue << endl;
			   }
*/


			   // Update unsmoothed value slice
			   //
		           if(idxRate < mKnockIO->mBarrierLo[idx]*
					(1. + BARRIER_TOL)         ||
		              idxRate > mKnockIO->mBarrierHi[idx]*
					(1. - BARRIER_TOL))
		    	   {
			       	(*mValueUS)[tsNode]  = exerNodeValue;	

				if (isKnockOut)
			       	    (*mRebateValueUS)[tsNode]
						= (*rebateTS)[tsNode];
			   }

/*
			   if (debugLevel == DEBUG_LEVEL_KIO)
			   {
			   	dppLog << *mValueUS << endl;
			   }
*/

		    	}	
		    }
		    else
		    	throw KFailure("%s: invalid knocking type (%d).\n",
				       routine,
				       mKnockIO->mIOWindow);	
		}    // if double smoothing




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



	    } 	// if modelObsDate

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

			double	pe,			// prob exercise
				te,			// expected exer time
				se,			// var exer time
				tZ;			// last exer time
			TDate	lastDate = mKnockIO->mObsDates.back();

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
			double	koPV   = mResults["UNDER"]
				       - mResults["PV"]
				       + mRebateValue->GetCenter();

			mResults["PV"] = koPV;
		}

	}


	//
	// Done updating, set flag
	//
	UpdateDone();


	delete indexRate;
	delete rateTS;

    }
    catch (KFailure) {

	delete indexRate;
	delete rateTS;

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


//---------------------------------------------------------------

KTSlice&
KVPToolKnockIO::GetValue()
{
	bool	isKnockOut = false;
	if (mKnockIO->mIOType == CRX_KNOCK_OUT)
		isKnockOut = true;


	// Calculate value of the underlying
	// by summing all dependencies.
	// Correct for knock-out using parity, plus the rebate.
	//
	if (mGetValue == NULL) {
		mGetValue = new KTSlice(*mVTree, "value", mCurveName);
		ASSERT_OR_THROW(mGetValue != NULL);
	}

	if (isKnockOut) {
		*mGetValue = 0e0;
		for (int idxD=0; idxD<this->NumDep(); idxD++) {
		    if (!(Dep(idxD)->IsType("KVPToolRate")))
			*mGetValue += Dep(idxD)->GetValue();
		}
		if (mValue!=NULL)
			*mGetValue -= *mValue;
			
		*mGetValue += *mRebateValue;
	} else {
		if (mValue!=NULL)
			*mGetValue = *mValue;
		else
			*mGetValue = 0e0;
	}

	return(*mGetValue);
}



